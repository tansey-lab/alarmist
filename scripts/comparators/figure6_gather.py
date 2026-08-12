#!/usr/bin/env python
"""Gather every number Figure 6 needs into tidy CSVs.

Figure 6 = the CCC comparator benchmark panel: seven competing cell-cell
communication methods vs ALARMIST on the GBM/LGG TMA (13 cores, 100,197 cells).

Reads only files that already exist on disk. Runs no inference, refits nothing.
Writes results/comparators/_benchmark/figure6/panel_{a,b,c,d,e,f}.csv.

Env: /Users/jiayifan/anaconda3/envs/bptf/bin/python
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/Users/jiayifan/tansey_lab/alarmist")
CMP = ROOT / "results/comparators"
OUT = ROOT / "results/comparators/_benchmark/figure6"
OUT.mkdir(parents=True, exist_ok=True)

# The two arms of ALARMIST motif 1 (mGAM <-> MES-like), as each method spells them.
ARM1 = "GRN->SORT1"
ARM2 = "ANXA1->FPR1"

rows_a: list[dict] = []


# Where the RANKING comes from. Only stLearn ships a ranked list; for everyone else
# we sorted something. For COMMOT and NICHES the method emits no per-interaction
# scalar at all, so the quantity being sorted is ours too -- that must be disclosed,
# because it is why those two rank so high in panel a.
#   native  = the method's own ranked output, read as-is
#   sorted  = the method's own per-interaction statistic, ranked by us
#   derived = no per-interaction statistic exists; we defined the ranked quantity
def add(method: str, ranked_object: str, arm: str, rank: float, n: float,
        prov: str = "sorted", note: str = "") -> None:
    """A whole-slide method: one run, one rank, one denominator."""
    rows_a.append(
        dict(method=method, ranked_object=ranked_object, arm=arm, rank_provenance=prov,
             rank=rank, n_tested=n, percentile=100.0 * rank / n,
             per_core=False, n_cores=np.nan, pct_q25=np.nan, pct_q75=np.nan,
             pct_min=np.nan, pct_max=np.nan, note=note)
    )


def add_per_core(method: str, ranked_object: str, arm: str,
                 ranks: np.ndarray, denoms: np.ndarray,
                 prov: str = "sorted", note: str = "") -> None:
    """A per-core method: percentile is computed WITHIN each core, then summarised.

    Ranks are never pooled across cores -- each core re-runs its own expression
    filter, so the denominators differ (SpatialDM: 686-1,661) and a pooled rank
    would mix incomparable universes. The plotted point is the median of the
    per-core percentiles and the whisker is their interquartile range.
    """
    pct = 100.0 * np.asarray(ranks, float) / np.asarray(denoms, float)
    rows_a.append(
        dict(method=method, ranked_object=ranked_object, arm=arm, rank_provenance=prov,
             rank=float(np.median(ranks)), n_tested=float(np.median(denoms)),
             percentile=float(np.median(pct)),
             per_core=True, n_cores=int(len(pct)),
             pct_q25=float(np.percentile(pct, 25)), pct_q75=float(np.percentile(pct, 75)),
             pct_min=float(pct.min()), pct_max=float(pct.max()), note=note)
    )


# ---------------------------------------------------------------- panel a
# Rank of each arm within the set of interactions THAT METHOD itself tested.
# Percentile, never a raw count: the denominators differ ~30-fold (panel b).

# CytoSignal -- row position in signif_summary (file is pre-sorted by n_hq desc)
cs = pd.read_csv(CMP / "cytosignal/GBM/cellchatdb2/run_full/quant/signif_summary_diffusion_Raw_smooth.csv")
n_cs = len(cs)
for arm, name in ((ARM1, "GRN - SORT1"), (ARM2, "ANXA1 - FPR1")):
    add("CytoSignal", "significant-cell count (diffusion)", arm,
        int(cs.index[cs["name"] == name][0]) + 1, n_cs, prov="sorted")

# stLearn -- rank by number of BH-significant spots
sl = pd.read_csv(CMP / "stlearn/GBM/cellchatdb2/data/lr_summary.csv", index_col=0)
sl_rank = sl["n_spots_sig"].rank(ascending=False, method="min")
n_sl = len(sl)
for arm, name in ((ARM1, "GRN_SORT1"), (ARM2, "ANXA1_FPR1")):
    # the one native ranking in the benchmark: stLearn writes uns['lr_summary'] itself
    add("stLearn", "significant-spot count", arm, float(sl_rank.loc[name]), n_sl, prov="native")

# SpatialDM -- per-core Moran's R rank, median of per-core percentiles
sd = pd.read_csv(CMP / "spatialdm/GBM/cellchatdb2/per_split_summary.csv")
for arm, col in ((ARM1, "GRN_SORT1_rank"), (ARM2, "ANXA1_FPR1_rank")):
    ok = sd[col].notna()
    add_per_core("SpatialDM", "bivariate Moran's R (per core)", arm,
                 sd.loc[ok, col].to_numpy(), sd.loc[ok, "n_pairs_valid"].to_numpy(),
                 prov="sorted", note=f"{int(ok.sum())}/13 cores testable; denominator 686-1,661")

# COMMOT -- per-core transported-signal rank (13/13 cores after the rerun)
cm = pd.read_csv(CMP / "commot/GBM/cellchatdb2/per_split_summary.csv")
for arm, col in ((ARM1, "GRN_SORT1_rank"), (ARM2, "ANXA1_FPR1_rank")):
    ok = cm[col].notna()
    add_per_core("COMMOT", "transported signal mass (per core)", arm,
                 cm.loc[ok, col].to_numpy(), cm.loc[ok, "n_pairs_used"].to_numpy(),
                 prov="derived", note=f"{int(ok.sum())}/13 cores; denominator fixed at 671")

# NICHES -- ALRA tier (the favourable configuration), median per-core rank
nc = pd.read_csv(CMP / "niches/GBM/cellchatdb2/summary_requested_lr.csv")
nc = nc[nc["tier"] == "alra"]
N_NICHES = 1088
for arm, name in ((ARM1, "GRN—SORT1"), (ARM2, "ANXA1—FPR1")):
    sub = nc[nc["lr"] == name]
    add_per_core("NICHES", "neighbourhood detection rate (per core)", arm,
                 sub["rank_by_frac"].to_numpy(), np.full(len(sub), N_NICHES),
                 prov="derived", note=f"ALRA tier; {len(sub)}/13 cores; denominator fixed at 1,088")

# CellChat -- summed communication probability, high-grade object
cc = pd.read_csv(CMP / "cellchat/GBM/default/quant/high_lr_ranked.csv")
n_cc = len(cc)
for arm, name in ((ARM1, "GRN_SORT1"), (ARM2, "ANXA1_FPR1")):
    add("CellChat", "communication probability (high grade)", arm,
        int(cc.index[cc["interaction_name"] == name][0]) + 1, n_cc, prov="sorted")

# LIANA+ -- inflow, max lr_mean over all sender/receiver pairs (file pre-sorted)
ln = pd.read_csv(CMP / "liana/GBM/cellchatdb2_inflow/data/lr_ranking_all_pairs.csv")
n_ln = len(ln)
for arm, name in ((ARM1, "GRN^SORT1"), (ARM2, "ANXA1^FPR1")):
    add("LIANA+", "inflow specificity", arm,
        int(ln.index[ln["lr"] == name][0]) + 1, n_ln, prov="sorted")

# ALARMIST -- rank of each arm inside motif 1, among the detectable LR pairs.
# lri_motifs.csv is 122 MB; read only the four columns needed.
lm = pd.read_csv(ROOT / "results/GBM/bptf/lri_motifs.csv",
                 usecols=["motif_idx", "ligand", "receptor", "factor_lrnorm"])
m1 = lm[lm["motif_idx"] == 1]
lr_w = m1.groupby(["ligand", "receptor"])["factor_lrnorm"].sum().sort_values(ascending=False)
n_al = len(lr_w)
for arm, key in ((ARM1, ("GRN", "SORT1")), (ARM2, ("ANXA1", "FPR1"))):
    add("ALARMIST", "motif-1 LRI weight", arm,
        int(lr_w.index.get_loc(key)) + 1, n_al, prov="derived", note="rank within motif 1")

pa = pd.DataFrame(rows_a)
pa.to_csv(OUT / "panel_a_recovery.csv", index=False)

# ---------------------------------------------------------------- panel b
# The denominators, and WHY each is what it is.
pb = pd.DataFrame([
    ("CytoSignal", n_cs, "1,088 scored; 895 with >=1 significant cell"),
    ("stLearn", n_sl, "526 of 3,233 DB rows: cannot encode heteromeric complexes"),
    ("SpatialDM", float(sd["n_pairs_valid"].median()), "per core; expression filter re-run each core"),
    ("COMMOT", float(cm["n_pairs_used"].median()), "per core; min_cell=100 filter"),
    ("NICHES", N_NICHES, "fixed mechanism set on the 5,119-gene panel"),
    ("CellChat", n_cc, "over-expression pre-filter, high-grade object"),
    ("LIANA+", n_ln, "inflow feature set, 633 distinct LR pairs"),
    ("ALARMIST", n_al, f"{n_al} LR pairs -> 25,271 sender|receiver|L|R|mode columns"),
], columns=["method", "n_tested", "reason"])
pb.to_csv(OUT / "panel_b_denominators.csv", index=False)

# ---------------------------------------------------------------- panel c
# Which OBJECT of each method's output does each arm live in, and is it the same one?
# obj_type: flat = no grouping object exists at all
#           curated = a grouping exists but its membership is a database column
#           learned = a grouping was estimated from the data
pc = pd.DataFrame([
    ("CytoSignal", "flat", "row 66 of one ranked table", "row 255 of the same table", False,
     "no grouping object of any kind"),
    ("stLearn", "flat", "per_lr_cci_cell_type/GRN_SORT1.csv", "per_lr_cci_cell_type/ANXA1_FPR1.csv", False,
     "one file per LR; most aggregated output pools all 526"),
    ("SpatialDM", "flat", "row of global_res.csv", "row of global_res.csv", False,
     "SparseAEH pattern clustering is an optional package, not installed"),
    ("COMMOT", "curated", "CellChatDB pathway GRN", "CellChatDB pathway ANNEXIN", False,
     "grouping is a database column, not estimated"),
    ("NICHES", "flat", "receiver block Vascular", "receiver block mGAM", False,
     "one-vs-rest marker ranking, not a co-activity model"),
    ("CellChat", "curated", "NMF outgoing Pattern 3 (GRN 0.604)", "NMF incoming Pattern 2 (ANNEXIN 1.000)", False,
     "members are pathways; outgoing/incoming are two separate factorizations"),
    ("LIANA+", "learned", "NMF inflow Factor 1 (Glial-Neuronal 92.9%)", "NMF inflow Factor 3 (MES-like 82.5%)", False,
     "factors organise by sender identity; 6/7 are >=75% one sender"),
    ("ALARMIST", "learned", "motif 1", "motif 1", True,
     "single factor over a patch x LRI matrix"),
], columns=["method", "obj_type", "arm1_object", "arm2_object", "same_object", "note"])
pc.to_csv(OUT / "panel_c_objects.csv", index=False)

# ---------------------------------------------------------------- panel d
# The mechanism: co-occurrence of the two arms at the cell vs at the 50 um patch.
d = json.loads((CMP / "liana/GBM/vs_alarmist/why_no_mgam_motif.json").read_text())
pd_rows = []
for key, label in (("cell_level", "cell"), ("patch_level", "50 um patch")):
    s = d[key]
    pd_rows.append(dict(
        unit=label, n=s["n"], arm1_pct=s["arm1_pct"], arm2_pct=s["arm2_pct"],
        both_pct=s["both_pct"], pearson=s["pearson"], spearman=s["spearman"],
        enrichment=d["enrichment_cell"] if key == "cell_level" else d["enrichment_patch"],
        source=s["unit"],
    ))
pdf = pd.DataFrame(pd_rows)
pdf.to_csv(OUT / "panel_d_unit.csv", index=False)

# ---------------------------------------------------------------- panel e
# Grade test, annotated by the replicate unit each method can actually test at.
pe = pd.DataFrame([
    ("ALARMIST", "TMA punch", 13, 0.013986, 0.016454, "MWU on per-punch motif-1 ON fraction"),
    ("SpatialDM", "TMA punch", 13, 0.6196, np.nan, "native differential_test; 0/1,662 at FDR<0.1"),
    ("LIANA+", "TMA punch", 13, 0.945, 1.0, "hand-rolled MWU on punch means (no native mode)"),
    ("CytoSignal", "cell (punch as random intercept)", 100197, 0.7054, 0.8043, "NEBULA NB-GLMM"),
    ("CellChat", "2 pooled objects", 2, np.nan, np.nan, "arm appears in neither net_up nor net_down"),
    ("NICHES", "cell (pooled)", 100197, np.nan, np.nan, "test.use='roc' emits no p-value at all"),
    ("stLearn", "none", 0, np.nan, np.nan, "no multi-sample or differential mode exists"),
    ("COMMOT", "none", 0, np.nan, np.nan, "no native differential mode"),
], columns=["method", "replicate_unit", "n", "p_arm1", "q_arm1", "note"])
pe.to_csv(OUT / "panel_e_grade.csv", index=False)

# honesty strip: all 20 ALARMIST motifs' grade p-values
mf = pd.read_csv(CMP / "cytosignal/GBM/cellchatdb2/nebula_grade/alarmist_motif_fraction_grade.csv")
mf.to_csv(OUT / "panel_e_all_motifs.csv", index=False)

# ---------------------------------------------------------------- panel f
# Density control. High-grade cores hold 3.4x more cells than low-grade ones, and
# ALARMIST's loading is a projection over a 50 um neighbourhood -- so grade
# separation is expected under a pure density effect. Test it with the continuous
# per-cell loadings (no GMM refit, no ON/OFF call).
import h5py
from scipy.stats import mannwhitneyu, spearmanr

U = np.load(ROOT / "results/GBM/single_cell/cell_loadings.npy")  # (100197, 20)

with h5py.File(ROOT / "data/xenium_mm_final_cell_id.h5ad", "r") as f:
    def cat(name: str) -> np.ndarray:
        g = f["obs"][name]
        if isinstance(g, h5py.Group):  # modern categorical
            codes = g["codes"][:]
            cats = np.asarray([c.decode() if isinstance(c, bytes) else str(c) for c in g["categories"][:]])
            return cats[codes]
        arr = g[:]
        return np.asarray([a.decode() if isinstance(a, bytes) else str(a) for a in arr])

    tma = cat("tma_id")
    grade = cat("grade")
    x = f["obs"]["centroid_x"][:].astype(float)
    y = f["obs"]["centroid_y"][:].astype(float)

assert len(tma) == U.shape[0], (len(tma), U.shape)

core_rows = []
for c in np.unique(tma):
    m = tma == c
    xc, yc = x[m], y[m]
    # convex-hull-free area proxy: the bounding box of the punch, in mm^2
    area_mm2 = ((xc.max() - xc.min()) * (yc.max() - yc.min())) / 1e6
    core_rows.append(dict(
        tma_id=c, grade=grade[m][0], n_cells=int(m.sum()),
        area_mm2=area_mm2, density_per_mm2=m.sum() / area_mm2,
        **{f"mean_loading_m{k}": float(U[m, k].mean()) for k in range(U.shape[1])},
    ))
cores = pd.DataFrame(core_rows).sort_values("tma_id")
cores.to_csv(OUT / "panel_f_cores.csv", index=False)

hi = cores["grade"] == "high"
dens = np.log10(cores["density_per_mm2"].to_numpy())
f_rows = []
for k in range(U.shape[1]):
    v = cores[f"mean_loading_m{k}"].to_numpy()
    p_grade = mannwhitneyu(v[hi.to_numpy()], v[~hi.to_numpy()], alternative="two-sided").pvalue
    rho_d, p_d = spearmanr(v, dens)
    # grade effect after removing the density trend (residualise loading on log density)
    b = np.polyfit(dens, v, 1)
    res = v - np.polyval(b, dens)
    p_adj = mannwhitneyu(res[hi.to_numpy()], res[~hi.to_numpy()], alternative="two-sided").pvalue
    f_rows.append(dict(motif=k, p_grade_raw=p_grade, rho_density=rho_d, p_density=p_d,
                       p_grade_density_adjusted=p_adj,
                       mean_high=float(v[hi.to_numpy()].mean()), mean_low=float(v[~hi.to_numpy()].mean())))
pf = pd.DataFrame(f_rows)
pf.to_csv(OUT / "panel_f_density.csv", index=False)

rho_gd, p_gd = spearmanr(hi.to_numpy().astype(float), dens)
(OUT / "panel_f_meta.json").write_text(json.dumps(dict(
    n_cores=int(len(cores)), n_high=int(hi.sum()), n_low=int((~hi).sum()),
    mean_cells_high=float(cores.loc[hi, "n_cells"].mean()),
    mean_cells_low=float(cores.loc[~hi, "n_cells"].mean()),
    density_grade_spearman=float(rho_gd), density_grade_p=float(p_gd),
    mwu_floor_7v6=0.0011655,
), indent=2))

print(f"wrote {OUT}")
for f_ in sorted(OUT.iterdir()):
    print(" ", f_.name)
