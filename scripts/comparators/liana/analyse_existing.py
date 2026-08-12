#!/usr/bin/env python
"""Replicate-aware (per-TMA-punch) analysis of the LIANA+ GBM results that already exist.

WHY THIS EXISTS
---------------
`comparator-benchmark/SKILL.md` requires each method to be exercised in its NATIVE multi-sample /
differential mode, and for GBM specifies: split by `obs['grade']` with the 13 `obs['tma_id']` cores
as the UNITS. What the LIANA run actually did was
`li.mt.compute_global_specificity(groupby='grade')`, which permutes labels across CELLS -- a
per-group specificity test, not a contrast, and pseudoreplicated by ~4 orders of magnitude
relative to 13 cores. LIANA has no native paired/differential spatial mode to substitute, so the
contrast is done here, on the units the study design actually has: **7 high-grade vs 6 low-grade
punches**, one number per punch.

Everything is READ-ONLY on the existing run. Nothing under an existing `data/` or `plots/`
directory is touched; all output goes to --out-dir.

INPUTS (all already on disk)
  <nmf-dir>/data/nmf_WH.npz                    W (locations x factors), H, feature names
  <inflow-dir>/data/cell_meta.csv              cell, x, y, cell_type, core
  <inflow-dir>/data/inflow_scores.npz          per-cell inflow, for the required-LR test
  <inflow-dir>/data/global_interactions.csv    source x target x LR, for circle/tile plots
  <h5ad>                                       ONLY for --group-col; ALARMIST outputs are refused

Usage:
  python analyse_existing.py --nmf-dir DIR --inflow-dir DIR --h5ad PATH --out-dir DIR
"""
import argparse, json, math, os, subprocess, time
import numpy as np, pandas as pd
import scanpy as sc, anndata as ad
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

p = argparse.ArgumentParser()
p.add_argument("--nmf-dir", required=True)
p.add_argument("--inflow-dir", required=True)
p.add_argument("--h5ad", required=True)
p.add_argument("--out-dir", required=True)
p.add_argument("--group-col", default="grade", help="the CONTRAST, e.g. grade (high vs low)")
p.add_argument("--punch-col", default="tma_id", help="the REPLICATE UNIT")
p.add_argument("--cell-type-col", default="cell_type")
p.add_argument("--required-lrs", default="GRN^SORT1,ANXA1^FPR1")
p.add_argument("--alpha", type=float, default=0.05)
p.add_argument("--seed", type=int, default=0)
p.add_argument("--no-liana-plots", action="store_true", default=False,
               help="skip li.pl.circle_plot / li.pl.tileplot (they need the comp-liana env)")
a = p.parse_args()
np.random.seed(a.seed)

# ALARMIST outputs must never enter a comparator run -- refuse them by name, as the sibling
# downstream script does.
for banned in ("motif", "patch_id"):
    if banned in (a.group_col, a.punch_col, a.cell_type_col):
        raise SystemExit(f"STOP: '{banned}' is an ALARMIST output and must not be used here.")

OUT = a.out_dir
DATA, PLOTS = os.path.join(OUT, "data"), os.path.join(OUT, "plots")
for d in (DATA, PLOTS): os.makedirs(d, exist_ok=True)
t0 = time.time()
log = lambda *m: print(f"[{time.strftime('%H:%M:%S')}]", *m, flush=True)
SAVED, SKIPPED = [], []


def save(fig, name):
    """Refuse blank figures (the plotnine/plt.gcf() trap that produced silent white PNGs)."""
    if fig is None or not fig.axes or not any(
            ax.has_data() or ax.get_images() or ax.collections or ax.patches or ax.texts
            for ax in fig.axes):
        plt.close("all"); SKIPPED.append((name, "blank")); log(f"  BLANK, not written: {name}")
        return
    fig.savefig(os.path.join(PLOTS, name), dpi=180, bbox_inches="tight")
    plt.close("all"); SAVED.append(name)


def bh(p):
    """Benjamini-Hochberg. NaN-safe; returns q in the input order."""
    p = np.asarray(p, float)
    ok = ~np.isnan(p)
    q = np.full_like(p, np.nan)
    if not ok.any():
        return q
    v = p[ok]; n = v.size
    o = np.argsort(v)
    qs = np.minimum.accumulate((v[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    out = np.empty(n); out[o] = np.clip(qs, 0, 1)
    q[ok] = out
    return q


# ------------------------------------------------------------------ inputs
log("loading")
meta = pd.read_csv(os.path.join(a.inflow_dir, "data", "cell_meta.csv"))
meta["cell"] = meta["cell"].astype(str)
if "core" not in meta.columns:
    raise SystemExit("STOP: cell_meta.csv has no 'core' column -- rerun the upstream script.")

z = np.load(os.path.join(a.nmf_dir, "data", "nmf_WH.npz"), allow_pickle=True)
W = np.asarray(z["W"], float)
feats = z["features"].astype(str)
cells = z["cells"].astype(str)
if W.shape[0] != len(meta):
    raise SystemExit(f"STOP: W has {W.shape[0]} rows, cell_meta.csv has {len(meta)}")
if not (cells == meta["cell"].values).all():
    raise SystemExit("STOP: nmf_WH.npz cell order does not match cell_meta.csv -- join by ID.")
factors = [f"Factor{i+1}" for i in range(W.shape[1])]
log(f"  W {W.shape} | {len(feats)} features | {len(factors)} factors")

# --group-col comes from the source h5ad, joined BY CELL ID. Only the two columns we need are
# read across; nothing else from that object enters this analysis.
src = sc.read_h5ad(a.h5ad)
for c in (a.group_col, a.punch_col):
    if c not in src.obs:
        raise SystemExit(f"STOP: '{c}' not in the h5ad obs; have {list(src.obs.columns)}")
grp = (src.obs[[a.group_col, a.punch_col]].astype(str)
          .set_axis(src.obs_names.astype(str)))
del src
meta = meta.join(grp, on="cell")
missing = int(meta[a.group_col].isna().sum())
if missing:
    log(f"  WARNING: {missing} cells have no {a.group_col}; dropped from the punch table")

# ------------------------------------------------- punch table (the replicate unit)
key = meta[a.punch_col].fillna(meta["core"]).astype(str)
ok = meta[a.group_col].notna().values
punch_group = (pd.DataFrame({"punch": key[ok], "group": meta[a.group_col][ok]})
                 .drop_duplicates())
bad = punch_group.groupby("punch")["group"].nunique()
if (bad > 1).any():
    raise SystemExit(f"STOP: punches spanning >1 {a.group_col}: {list(bad[bad>1].index)} -- "
                     f"the punch is not a valid replicate unit for this contrast.")
punch_group = punch_group.set_index("punch")["group"]
log(f"  {punch_group.nunique()} groups over {len(punch_group)} punches: "
    f"{punch_group.value_counts().to_dict()}")
if len(punch_group) < 4:
    raise SystemExit("STOP: fewer than 4 punches -- a rank test on this is not worth running.")


def punch_means(mat, colnames, tag):
    """One number per punch per column: the mean over that punch's cells."""
    df = pd.DataFrame(mat[ok], columns=colnames)
    df["punch"] = key[ok].values
    out = df.groupby("punch").mean()
    out = out.loc[punch_group.index]
    out.to_csv(os.path.join(DATA, f"punch_{tag}_means.csv"))
    return out


def contrast(pm, tag):
    """Mann-Whitney U on the PUNCH means -- n = number of punches, not number of cells."""
    lvls = sorted(punch_group.unique())
    if len(lvls) != 2:
        log(f"  {a.group_col} has {len(lvls)} levels; only a 2-level contrast is implemented")
        return None
    A, B = [pm.index[punch_group.loc[pm.index] == L] for L in lvls]
    rows = []
    for c in pm.columns:
        x, y = pm.loc[A, c].values, pm.loc[B, c].values
        try:
            u, pv = mannwhitneyu(x, y, alternative="two-sided")
        except ValueError:
            u, pv = np.nan, np.nan
        rows.append({"feature": c, f"mean_{lvls[0]}": float(np.mean(x)),
                     f"mean_{lvls[1]}": float(np.mean(y)),
                     "log2fc": float(np.log2((np.mean(x) + 1e-12) / (np.mean(y) + 1e-12))),
                     "U": float(u) if u == u else np.nan, "pval": float(pv) if pv == pv else np.nan,
                     f"n_{lvls[0]}": len(x), f"n_{lvls[1]}": len(y)})
    res = pd.DataFrame(rows)
    res["qval_bh"] = bh(res["pval"].values)
    res = res.sort_values("pval")
    res.to_csv(os.path.join(DATA, f"punch_{tag}_tests.csv"), index=False)
    n_sig = int((res["qval_bh"] < a.alpha).sum())
    log(f"  {tag}: {n_sig}/{len(res)} at BH q<{a.alpha} "
        f"(min p {res['pval'].min():.4g}); n = {len(A)} vs {len(B)} PUNCHES")
    # The honest caveat, printed every run: with 7 vs 6 the smallest attainable two-sided
    # Mann-Whitney p is 1/C(13,6) * 2 = 0.0012, and BH over many features will rarely survive.
    log(f"  (floor on a two-sided rank test with {len(A)} vs {len(B)} units: "
        f"p_min = {2.0 / math.comb(len(A) + len(B), min(len(A), len(B))):.4g})")
    return res


# ------------------------------------------------------- 1. factors vs grade
log(f"[1] NMF factor usage per punch, contrasted on {a.group_col}")
pm_fac = punch_means(W, factors, "factor")
res_fac = contrast(pm_fac, "factor")


def strip(pm, res, tag, title):
    lvls = sorted(punch_group.unique())
    cols = list(res["feature"]) if res is not None else list(pm.columns)
    ncol = min(4, len(cols)); nrow = int(np.ceil(len(cols) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.1 * ncol, 3.0 * nrow), squeeze=False)
    axes = axes.ravel()
    for i, c in enumerate(cols):
        d = pd.DataFrame({"v": pm[c].values, "g": punch_group.loc[pm.index].values,
                          "punch": pm.index})
        sns.boxplot(data=d, x="g", y="v", order=lvls, ax=axes[i], width=.55,
                    showfliers=False, boxprops=dict(alpha=.35))
        sns.stripplot(data=d, x="g", y="v", order=lvls, ax=axes[i], size=6,
                      edgecolor="black", linewidth=.6, jitter=.12)
        sub = "" if res is None else (
            f"\np={float(res.loc[res.feature == c, 'pval'].iloc[0]):.3g}  "
            f"q={float(res.loc[res.feature == c, 'qval_bh'].iloc[0]):.3g}")
        axes[i].set_title(f"{c}{sub}", fontsize=8)
        axes[i].set_xlabel(""); axes[i].set_ylabel("punch mean", fontsize=7)
        axes[i].tick_params(labelsize=7)
    for j in range(len(cols), len(axes)): axes[j].axis("off")
    fig.suptitle(f"{title}  —  one point per TMA punch (n={len(pm)}), "
                 f"Mann-Whitney on punch means", fontsize=10)
    fig.tight_layout()
    save(fig, f"punch_{tag}_by_{a.group_col}.png")


strip(pm_fac, res_fac, "factor", "NMF factor usage")

# ------------------------------------------- 2. the required LRs, same units
log(f"[2] required LRs per punch, contrasted on {a.group_col}")
required = [s.strip() for s in a.required_lrs.split(",") if s.strip()]
inf_path = os.path.join(a.inflow_dir, "data", "inflow_scores.npz")
req_report = {}
if os.path.exists(inf_path):
    iz = np.load(inf_path, allow_pickle=True)
    ifeat = iz["features"].astype(str)
    icells = iz["cells"].astype(str)
    if not (icells == meta["cell"].values).all():
        raise SystemExit("STOP: inflow_scores.npz cell order != cell_meta.csv")
    lr_of = np.array(["^".join(f.split("^")[1:]) for f in ifeat])
    cols, names = [], []
    IX = iz["values"]
    for lr in required:
        hit = np.where(lr_of == lr)[0]
        req_report[lr] = {"n_features": int(hit.size),
                          "features": sorted(ifeat[hit].tolist())} if hit.size else None
        if not hit.size:
            log(f"  REQUIRED {lr}: absent from the inflow feature space -- that is the result")
            continue
        # per-sender features, plus the summed "all senders" aggregate
        for i in hit:
            cols.append(IX[:, i]); names.append(ifeat[i])
        cols.append(IX[:, hit].sum(1)); names.append(f"ALL^{lr}")
        log(f"  REQUIRED {lr}: {hit.size} sender-typed features + 1 aggregate")
    if cols:
        pm_lr = punch_means(np.column_stack(cols), names, "requiredLR")
        res_lr = contrast(pm_lr, "requiredLR")
        agg = [n for n in names if n.startswith("ALL^")]
        if agg:
            strip(pm_lr[agg], res_lr[res_lr.feature.isin(agg)] if res_lr is not None else None,
                  "requiredLR", "required-LR inflow (summed over senders)")
else:
    log(f"  inflow_scores.npz absent at {inf_path} -- required-LR punch test skipped")
    SKIPPED.append(("required-LR punch test", "inflow_scores.npz absent"))

# ------------------------------ 3. li.pl plots never drawn, from existing tables
gi_path = os.path.join(a.inflow_dir, "data", "global_interactions.csv")
if not a.no_liana_plots and os.path.exists(gi_path):
    log("[3] li.pl.circle_plot / li.pl.tileplot from global_interactions.csv")
    import liana as li
    gi = pd.read_csv(gi_path)
    gi["lr"] = gi["ligand_complex"] + "^" + gi["receptor_complex"]
    sig = gi[gi["pval"] < a.alpha].copy()
    log(f"  {len(sig):,} of {len(gi):,} rows at p<{a.alpha}")
    # liana/plotting/_common.py:16 asserts `uns_key in adata.uns` WHENEVER an adata is passed,
    # even if liana_res was supplied explicitly. So: tileplot/dotplot take liana_res alone
    # (no adata); circle_plot requires adata positionally, so the frame is parked in .uns instead.
    dummy = ad.AnnData(X=np.zeros((len(meta), 1), np.float32),
                       obs=pd.DataFrame({a.cell_type_col: meta["cell_type"].astype("category").values},
                                        index=meta["cell"].values))
    dummy.uns["liana_res"] = sig

    def _gg(fn, name, w=9, h=7, use_adata=False, **kw):
        # liana's plotting entry points disagree on return type: tileplot/dotplot return a
        # PLOTNINE ggplot (own .save), circle_plot returns a matplotlib Axes. Handle both rather
        # than assuming -- assuming is what produced blank PNGs elsewhere in this suite.
        try:
            g = (fn(dummy, uns_key="liana_res", **kw) if use_adata
                 else fn(liana_res=sig, **kw))
            out = os.path.join(PLOTS, name)
            if hasattr(g, "save"):
                g.save(out, dpi=180, width=w, height=h, verbose=False)
            else:
                fig = g.get_figure() if hasattr(g, "get_figure") else g
                fig.set_size_inches(w, h)
                fig.savefig(out, dpi=180, bbox_inches="tight"); plt.close(fig)
            SAVED.append(name); log(f"  ok: {name}")
        except Exception as e:
            plt.close("all"); log(f"  FAILED {name}: {type(e).__name__}: {e}")
            SKIPPED.append((name, f"{type(e).__name__}: {e}"))

    _gg(li.pl.circle_plot, "circle_plot_top50.png", w=9, h=9, use_adata=True,
        groupby=a.cell_type_col, score_key="lr_mean", top_n=50, orderby="lr_mean",
        orderby_ascending=False)
    _gg(li.pl.tileplot, "tileplot_top25.png", w=11, h=9, fill="lr_mean", label="pval",
        top_n=25, orderby="lr_mean", orderby_ascending=False)
    _gg(li.pl.dotplot, "dotplot_top25.png", w=9, h=8, colour="lr_mean", size="pval",
        inverse_size=True, top_n=25, orderby="lr_mean", orderby_ascending=False)
    for lr in required:
        lig, rec = lr.split("^", 1)
        _gg(li.pl.tileplot, f"tileplot_{lr.replace('^','-')}.png", w=10, h=6,
            fill="lr_mean", label="pval", ligand_complex=lig, receptor_complex=rec,
            top_n=25, orderby="lr_mean", orderby_ascending=False)
else:
    SKIPPED.append(("liana plots", "skipped or global_interactions.csv absent"))

# ------------------------------------------------------------------ manifest
try:
    _sha = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                          capture_output=True, text=True).stdout.strip() or None
except Exception:
    _sha = None
json.dump({"script": "analyse_existing.py", "dataset": "GBM", "tier": "cellchatdb2",
           "nmf_dir": a.nmf_dir, "inflow_dir": a.inflow_dir, "h5ad": a.h5ad,
           "group_col": a.group_col, "punch_col": a.punch_col,
           "replicate_unit": f"{a.punch_col} (n={len(punch_group)})",
           "group_sizes": punch_group.value_counts().to_dict(),
           "n_factors": len(factors), "alpha": a.alpha,
           "n_factors_signif_bh": (None if res_fac is None
                                   else int((res_fac["qval_bh"] < a.alpha).sum())),
           "required_lrs": required, "required_lr_report": req_report,
           "figures_saved": SAVED, "skipped": SKIPPED, "seed": a.seed, "git_sha": _sha,
           "note": "Mann-Whitney on PUNCH means. This replaces the cell-level p-values from "
                   "compute_global_specificity, which are pseudoreplicated w.r.t. the cores.",
           "wall_min": round((time.time() - t0) / 60, 2)},
          open(os.path.join(OUT, "manifest.json"), "w"), indent=2, default=str)
log(f"DONE in {(time.time()-t0)/60:.2f} min -> {OUT} "
    f"({len(SAVED)} figures, {len(SKIPPED)} skipped)")
