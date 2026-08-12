#!/usr/bin/env python
"""Does COMMOT independently reproduce ALARMIST's motif 1 (the mGAM <-> MES-like loop)?

ALARMIST motif 1 is a *bidirectional* loop (CLAUDE.md): GRN->SORT1 carries mGAM -> MES-like,
ANXA1->FPR1 carries MES-like -> mGAM. COMMOT is a completely different model (collective optimal
transport, no factorization, no patches) but returns a per-cell amount SENT and RECEIVED for each
of those two pairs. That yields a falsifiable prediction rather than a vague "do they agree":

    GRN->SORT1 :  sender side enriched in mGAM      , receiver side enriched in MES-like
    ANXA1->FPR1:  sender side enriched in MES-like  , receiver side enriched in mGAM

and, per cell, motif-1 loading should track the mGAM side of both (s-GRN-SORT1 + r-ANXA1-FPR1).

Three controls keep this from being circular or trivial:
  1. ALL 20 motifs are correlated against the COMMOT motif-1 quantities. If motif 1 is not the
     argmax, the agreement is not specific.
  2. An unrelated but strong pair (FGF1-FGFR2, rank 1 in 7 of 13 cores) is correlated against
     motif 1. If motif 1 correlates with everything, the agreement is trivial.
  3. GRN and FPR1 are mGAM-specific *by expression* (CLAUDE.md), so the mGAM side of each
     direction is partly guaranteed. The non-trivial half is the COUNTERPART side: whether
     COMMOT's spatially-constrained transport lands SORT1-receipt on MES-like and ANXA1-sending
     on MES-like. That half is reported separately.

Replicate unit is the TMA core (7 high-grade + 6 low-grade), per spatial-workflow: per-cell
statistics are descriptive and are summarised ACROSS the 13 cores; grade is tested at core level.
COMMOT magnitudes are not comparable across cores (per-run OT normalisation, see DEVIATIONS.md),
so anything pooled is rank-normalised within core first.

Usage: python compare_motif1_commot.py [--commot-dir DIR] [--out-dir DIR] [--motif 1]
"""
import argparse, json, os, sys
from pathlib import Path
import numpy as np, pandas as pd, h5py
from scipy.stats import spearmanr, mannwhitneyu
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import plotting as _plot

ROOT = "/Users/jiayifan/tansey_lab/alarmist"
p = argparse.ArgumentParser()
p.add_argument("--h5ad", default=f"{ROOT}/data/xenium_mm_final_cell_id.h5ad")
p.add_argument("--loadings", default=f"{ROOT}/results/GBM/single_cell/cell_loadings.npy")
p.add_argument("--commot-dir", default=f"{ROOT}/results/comparators/commot/GBM/cellchatdb2")
p.add_argument("--out-dir", default=f"{ROOT}/results/comparators/commot/GBM/vs_alarmist")
p.add_argument("--motif", type=int, default=1, help="0-based column of cell_loadings.npy")
p.add_argument("--control-lr", default="FGF1-FGFR2")
p.add_argument("--sender-type", default="mGAM")
p.add_argument("--receiver-type", default="MES-like")
a = p.parse_args()
OUT = a.out_dir; os.makedirs(OUT, exist_ok=True)

# ---- plotting conventions: scripts/comparators/_common/plotting.py (CLAUDE.md house rules)
_plot.apply_publication_style(**{"figure.dpi": 150})


def save_all_formats(fig, path_no_ext):
    """Three formats at the figure's own dpi (150), closing as we go."""
    _plot.save_all_formats(fig, path_no_ext, dpi=None, close=True, verbose=True)

# ------------------------------------------------------------------ ALARMIST side
U = np.load(a.loadings)
f = h5py.File(a.h5ad, "r")
def col(n):
    g = f["obs"][n]
    if isinstance(g, h5py.Group):
        k = [x.decode() if isinstance(x, bytes) else str(x) for x in g["categories"][:]]
        return np.array([k[i] if i >= 0 else None for i in g["codes"][:]], dtype=object)
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in g[:]], dtype=object)
# the obs index is NOT always stored under '_index' -- anndata records its name in the attr
_ix = f["obs"].attrs.get("_index", "_index")
_ix = _ix.decode() if isinstance(_ix, bytes) else str(_ix)
obs_names = [x.decode() if isinstance(x, bytes) else str(x) for x in f["obs"][_ix][:]]
ct_all, grade_all, tma_all = col("cell_type"), col("grade"), col("tma_id")
f.close()
assert U.shape[0] == len(obs_names), "loadings / h5ad row mismatch"
row_of = {c: i for i, c in enumerate(obs_names)}

# Alignment guard, not decoration: motif `--motif` must peak on --sender-type. CLAUDE.md pins
# motif 1 = mGAM; if the row order or the column index were wrong this would not hold.
mean_by_ct = pd.Series(U[:, a.motif]).groupby(pd.Series(ct_all)).mean().sort_values(ascending=False)
if mean_by_ct.index[0] != a.sender_type:
    raise SystemExit(f"ALIGNMENT FAIL: motif {a.motif} peaks on '{mean_by_ct.index[0]}', "
                     f"expected '{a.sender_type}'. Refusing to continue.")
print(f"alignment OK: motif {a.motif} peaks on {a.sender_type} "
      f"({mean_by_ct.iloc[0]:.2e} vs next {mean_by_ct.index[1]} {mean_by_ct.iloc[1]:.2e})")

LR1, LR2 = "GRN-SORT1", "ANXA1-FPR1"
QUANT = [(f"s-{LR1}", "GRN->SORT1 sent"), (f"r-{LR1}", "GRN->SORT1 received"),
         (f"s-{LR2}", "ANXA1->FPR1 sent"), (f"r-{LR2}", "ANXA1->FPR1 received"),
         (f"s-{a.control_lr}", f"{a.control_lr} sent (control)"),
         (f"r-{a.control_lr}", f"{a.control_lr} received (control)"),
         # COMMOT's own grand total over all 671 pairs: the "this cell type is simply in dense
         # tissue and receives/sends everything" baseline. A per-pair cell-type profile that
         # merely reproduces this profile carries no pair-specific information.
         ("s-total-total", "all signalling sent (baseline)"),
         ("r-total-total", "all signalling received (baseline)")]

cores = sorted([d for d in os.listdir(a.commot_dir)
                if os.path.isdir(os.path.join(a.commot_dir, d))], key=lambda x: int(x))
print("cores:", cores)

rows_corr, rows_ct, rows_core, spatial_cache = [], [], [], {}
for c in cores:
    D = os.path.join(a.commot_dir, c, "data")
    meta = pd.read_csv(os.path.join(D, "cell_meta.csv"), dtype={"cell": str})
    snd = pd.read_csv(os.path.join(D, "sum_sender.csv.gz"), index_col=0)
    rcv = pd.read_csv(os.path.join(D, "sum_receiver.csv.gz"), index_col=0)
    idx = np.array([row_of[x] for x in meta.cell])
    u = U[idx]                                   # cells_in_core x 20
    ct = ct_all[idx]
    sig = {}
    for key, _ in QUANT:
        src = snd if key.startswith("s-") else rcv
        sig[key] = src[key].to_numpy() if key in src.columns else None
    # the mGAM-side composite: motif 1's loading peaks on mGAM, and mGAM is the sender of
    # GRN->SORT1 and the receiver of ANXA1->FPR1
    z = lambda v: pd.Series(v).rank(pct=True).to_numpy()
    if sig[f"s-{LR1}"] is not None and sig[f"r-{LR2}"] is not None:
        sig["composite-mGAM-side"] = z(sig[f"s-{LR1}"]) + z(sig[f"r-{LR2}"])
    # the composite's own baseline, built the SAME way from the grand totals, so the two are on
    # one scale and the residual below is meaningful
    if sig["s-total-total"] is not None and sig["r-total-total"] is not None:
        sig["composite-baseline"] = z(sig["s-total-total"]) + z(sig["r-total-total"])

    # (1) every motif vs every COMMOT quantity -- specificity control
    for key in list(sig):
        if sig[key] is None: continue
        for k in range(U.shape[1]):
            rho, _ = spearmanr(u[:, k], sig[key])
            rows_corr.append(dict(core=c, quantity=key, motif=k, rho=rho))

    # (2) directional test: within-core percentile of each signal, averaged per cell type
    for key in list(sig):
        if sig[key] is None: continue
        pct = pd.Series(sig[key]).rank(pct=True).to_numpy()
        for t, m in pd.Series(pct).groupby(pd.Series(ct)).mean().items():
            rows_ct.append(dict(core=c, quantity=key, cell_type=t, mean_pct=m,
                                n=int((ct == t).sum())))

    # (3) core-level summary for the grade contrast
    lr_tot = pd.read_csv(os.path.join(D, "lr_total_received.csv"), index_col=0).iloc[:, 0]
    lr_tot.index = [i[2:] if str(i).startswith("r-") else str(i) for i in lr_tot.index]
    order = list(lr_tot.sort_values(ascending=False).index)
    rows_core.append(dict(
        core=c, grade=grade_all[idx][0], n_cells=len(idx),
        motif_mean=float(u[:, a.motif].mean()),
        motif_mGAM_mean=float(u[ct == a.sender_type, a.motif].mean()) if (ct == a.sender_type).any() else np.nan,
        frac_mGAM=float((ct == a.sender_type).mean()),
        rank_GRN_SORT1=order.index(LR1) + 1 if LR1 in order else np.nan,
        rank_ANXA1_FPR1=order.index(LR2) + 1 if LR2 in order else np.nan,
        rho_composite=spearmanr(u[:, a.motif], sig["composite-mGAM-side"])[0]
                      if "composite-mGAM-side" in sig else np.nan))
    if c == "13":   # a mid-sized, high-grade core for the spatial panel
        spatial_cache = dict(x=meta.x.to_numpy(), y=meta.y.to_numpy(), ct=ct,
                             motif=u[:, a.motif], **{k: v for k, v in sig.items() if v is not None})
    print(f"  core {c}: {len(idx)} cells joined")

corr = pd.DataFrame(rows_corr); ctdf = pd.DataFrame(rows_ct); coredf = pd.DataFrame(rows_core)
corr.to_csv(f"{OUT}/motif_vs_commot_spearman_percore.csv", index=False)
ctdf.to_csv(f"{OUT}/celltype_signal_percentile_percore.csv", index=False)
coredf.to_csv(f"{OUT}/core_summary.csv", index=False)

# ---------------------------------------------------------------- summarise
med = corr.groupby(["quantity", "motif"]).rho.median().unstack()
med.to_csv(f"{OUT}/motif_vs_commot_spearman_median.csv")
print("\n=== median Spearman across 13 cores (COMMOT quantity x motif) ===")
print(med.round(3).to_string())
print("\n=== which motif does each COMMOT quantity match best? ===")
best = med.idxmax(axis=1)
for q in med.index:
    print(f"  {q:28s} -> motif {best[q]:2d} (rho={med.loc[q, best[q]]:.3f}); "
          f"motif {a.motif} rho={med.loc[q, a.motif]:.3f}")

print(f"\n=== directional test: mean within-core percentile by cell type ===")
piv = ctdf.groupby(["quantity", "cell_type"]).mean_pct.mean().unstack()
piv.to_csv(f"{OUT}/celltype_signal_percentile_mean.csv")
print(piv.round(3).to_string())
# --- BASELINE CORRECTION (do not skip) -------------------------------------------------------
# mGAM turns out to be COMMOT's most communicative cell type overall: it tops BOTH total-total
# baselines. So a raw "mGAM is top for s-GRN-SORT1" is not evidence of anything pair-specific --
# it is what mGAM does for every pair. Subtract the matching grand-total profile and judge the
# prediction on the residual, i.e. on signal that is SPECIFIC to this LR pair.
base = {"s": "s-total-total", "r": "r-total-total"}
BASELINES = list(base.values()) + ["composite-baseline"]
piv_corr = piv.copy()
for q in piv.index:
    if q in BASELINES: continue
    b = "composite-baseline" if q.startswith("composite") else base.get(q.split("-")[0])
    if b in piv.index:
        piv_corr.loc[q] = piv.loc[q] - piv.loc[b]
piv_corr = piv_corr.drop(index=[v for v in BASELINES if v in piv_corr.index])
piv_corr.to_csv(f"{OUT}/celltype_signal_percentile_baseline_corrected.csv")
print("\n=== SAME, baseline-corrected (minus the matching total-total profile) ===")
print(piv_corr.round(3).to_string())

pred = [(f"s-{LR1}", a.sender_type), (f"r-{LR1}", a.receiver_type),
        (f"s-{LR2}", a.receiver_type), (f"r-{LR2}", a.sender_type)]
print("\n  prediction -> rank of the predicted cell type   [raw]  ->  [baseline-corrected]")
verdict = {}
for q, want in pred:
    if q not in piv.index: continue
    o_raw = piv.loc[q].sort_values(ascending=False)
    o_cor = piv_corr.loc[q].sort_values(ascending=False)
    r_raw = list(o_raw.index).index(want) + 1
    r_cor = list(o_cor.index).index(want) + 1
    verdict[q] = dict(predicted=want, rank_raw=r_raw, rank_corrected=r_cor, of=len(o_raw),
                      top_raw=o_raw.index[0], top_corrected=o_cor.index[0],
                      value_corrected=float(piv_corr.loc[q, want]))
    print(f"    {q:16s} predict {want:9s}: rank {r_raw}/{len(o_raw)} (top {o_raw.index[0]:14s})"
          f" -> rank {r_cor}/{len(o_cor)} (top {o_cor.index[0]:14s}, "
          f"residual {piv_corr.loc[q, want]:+.3f})")

# grade contrast at the CORE level (7 high vs 6 low)
hi = coredf[coredf.grade == "high"]; lo = coredf[coredf.grade == "low"]
grade_tests = {}
for m in ["motif_mean", "motif_mGAM_mean", "rank_GRN_SORT1", "rho_composite", "frac_mGAM"]:
    x, y = hi[m].dropna(), lo[m].dropna()
    if len(x) > 1 and len(y) > 1:
        u_, pv = mannwhitneyu(x, y, alternative="two-sided")
        grade_tests[m] = dict(high_median=float(x.median()), low_median=float(y.median()),
                              U=float(u_), p=float(pv), n_high=len(x), n_low=len(y))
print(f"\n=== grade contrast, core = replicate (n={len(hi)} high vs {len(lo)} low) ===")
for m, d in grade_tests.items():
    print(f"  {m:18s} high {d['high_median']:.4g} vs low {d['low_median']:.4g}  p={d['p']:.3f}")

# ---------------------------------------------------------------- figures
q_order = [k for k, _ in QUANT if k in med.index] + (
    ["composite-mGAM-side"] if "composite-mGAM-side" in med.index else [])
fig, ax = plt.subplots(figsize=(9, 0.45 * len(q_order) + 2))
im = ax.imshow(med.loc[q_order].values, cmap="RdBu_r", aspect="auto",
               vmin=-np.abs(med.values).max(), vmax=np.abs(med.values).max())
ax.set_xticks(range(med.shape[1])); ax.set_xticklabels(med.columns, fontsize=7)
ax.set_yticks(range(len(q_order))); ax.set_yticklabels(q_order, fontsize=8)
ax.set_xlabel("ALARMIST motif (0-based column of cell_loadings.npy)", fontsize=9)
ax.set_title("Median Spearman across 13 TMA cores: COMMOT per-cell signal vs motif loading",
             fontsize=9)
for j, k in enumerate(med.columns):
    if k == a.motif:
        ax.add_patch(plt.Rectangle((j - .5, -.5), 1, len(q_order), fill=False,
                                   edgecolor="black", lw=1.6))
fig.colorbar(im, ax=ax, shrink=0.7, label="Spearman rho")
save_all_formats(fig, f"{OUT}/motif_vs_commot_heatmap")

sub = piv_corr.loc[[k for k in [f"s-{LR1}", f"r-{LR1}", f"s-{LR2}", f"r-{LR2}"] if k in piv.index]]
fig, ax = plt.subplots(figsize=(10, 4))
w = 0.8 / len(sub); types = list(piv.columns)
for i, q in enumerate(sub.index):
    hl = [a.sender_type if q.startswith("s-") and LR1 in q else
          a.receiver_type if q.startswith("r-") and LR1 in q else
          a.receiver_type if q.startswith("s-") and LR2 in q else a.sender_type][0]
    vals = sub.loc[q, [t for t in types if t in sub.columns]].values
    bars = ax.bar(np.arange(len(types)) + i * w, vals, w, label=q)
    for b, t in zip(bars, types):
        if t == hl: b.set_edgecolor("black"); b.set_linewidth(2)
ax.axhline(0.0, ls="--", c="grey", lw=0.8)
ax.set_xticks(np.arange(len(types)) + 0.4 - w / 2); ax.set_xticklabels(types, rotation=45, ha="right")
ax.set_ylabel("baseline-corrected percentile (minus total-total)")
ax.set_title("Directional test — black outline = the cell type ALARMIST motif "
             f"{a.motif} predicts", fontsize=9)
ax.legend(fontsize=7, ncol=2)
save_all_formats(fig, f"{OUT}/directional_celltype_test")

if spatial_cache:
    panels = [("motif", f"ALARMIST motif {a.motif} loading"),
              (f"s-{LR1}", "COMMOT GRN->SORT1 sent"),
              (f"r-{LR2}", "COMMOT ANXA1->FPR1 received")]
    panels = [(k, t) for k, t in panels if k in spatial_cache]
    fig, axes = plt.subplots(1, len(panels) + 1, figsize=(5 * (len(panels) + 1), 4.6))
    for ax_, (k, t) in zip(axes, panels):
        v = pd.Series(spatial_cache[k]).rank(pct=True)
        ax_.scatter(spatial_cache["x"], spatial_cache["y"], c=v, s=1.5, cmap="magma")
        ax_.set_title(t, fontsize=9); ax_.set_aspect("equal"); ax_.axis("off")
    m = spatial_cache["ct"] == a.sender_type
    axes[-1].scatter(spatial_cache["x"][~m], spatial_cache["y"][~m], c="lightgrey", s=1)
    axes[-1].scatter(spatial_cache["x"][m], spatial_cache["y"][m], c="crimson", s=3)
    axes[-1].set_title(f"{a.sender_type} cells (core 13)", fontsize=9)
    axes[-1].set_aspect("equal"); axes[-1].axis("off")
    fig.suptitle("Core 13 — within-core percentile, so panels are directly comparable", fontsize=9)
    save_all_formats(fig, f"{OUT}/spatial_concordance_core13")

json.dump({"script": "compare_motif1_commot.py", "motif": a.motif,
           "alignment_check": {"peaks_on": mean_by_ct.index[0], "expected": a.sender_type},
           "best_matching_motif_per_quantity": {q: int(best[q]) for q in med.index},
           "median_rho_motif_of_interest": {q: float(med.loc[q, a.motif]) for q in med.index},
           "directional_verdict": verdict, "grade_tests": grade_tests,
           "n_cores": len(cores), "n_cells": int(coredf.n_cells.sum())},
          open(f"{OUT}/summary.json", "w"), indent=2)
print(f"\nwrote -> {OUT}")
