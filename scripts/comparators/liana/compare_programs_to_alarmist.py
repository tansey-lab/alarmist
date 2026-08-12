#!/usr/bin/env python
"""Head-to-head: LIANA+ communication programs vs ALARMIST motifs, per cell.

This is the only fully like-for-like comparison available between the two methods. Both produce
a per-cell x per-program activity matrix over the SAME cells:

    ALARMIST   results/GBM/single_cell/cell_loadings.npy        (100,197 x 20)  motifs, non-negative
    MOFA-Flex  mofaflex_inflow/data/factor_scores.csv.gz        (100,190 x 17)  factors, SIGNED
    NMF        nmf_inflow/data/NMF_W_factor_scores.csv          (100,190 x  7)  factors, non-negative

so they can be aligned by cell name and correlated directly. Everything else in the benchmark
compares objects of different kinds.

Spearman is used throughout: BPTF loadings and inflow scores are heavy-tailed and ~99% zero, so
Pearson is dominated by a handful of cells. Factor SIGN is arbitrary for MOFA-Flex (it is a signed
factorisation, unlike BPTF and NMF which are non-negative), so |rho| is the quantity of interest
and the sign is reported only for reference.

The two diagnostics that matter, and why:
  * DEGENERACY -- if many ALARMIST motifs pick the SAME LIANA factor as their best match, the
    mapping is not a correspondence, it is a hub.
  * GENERAL-ACTIVITY CONFOUND -- a factor that tracks "how much signalling is here" will correlate
    with everything. Tested by correlating each factor against (a) the summed ALARMIST loading over
    all motifs and (b) the summed inflow score per cell.

Run:
    /Users/jiayifan/anaconda3/envs/bptf/bin/python \\
        scripts/comparators/liana/compare_programs_to_alarmist.py
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr

ROOT = Path("/Users/jiayifan/tansey_lab/alarmist")

p = argparse.ArgumentParser()
p.add_argument("--h5ad", default=str(ROOT / "data/xenium_mm_final_cell_id.h5ad"))
p.add_argument("--alarmist-loadings", default=str(ROOT / "results/GBM/single_cell/cell_loadings.npy"))
p.add_argument("--liana-dir", default=str(ROOT / "results/comparators/liana/GBM"))
p.add_argument("--out-dir", default=str(ROOT / "results/comparators/liana/GBM/vs_alarmist"))
p.add_argument("--motif-of-interest", type=int, default=1,
               help="0-based column of cell_loadings.npy. 1 = the mGAM<->MES-like loop (CLAUDE.md).")
a = p.parse_args()

OUT = Path(a.out_dir); (OUT / "data").mkdir(parents=True, exist_ok=True)
LI = Path(a.liana_dir)

ad = sc.read_h5ad(a.h5ad, backed="r")
names = ad.obs_names.astype(str).values
ct = pd.Series(ad.obs["cell_type"].astype(str).values, index=names)
U = pd.DataFrame(np.load(a.alarmist_loadings), index=names,
                 columns=[f"motif{i}" for i in range(np.load(a.alarmist_loadings).shape[1])])
MOI = f"motif{a.motif_of_interest}"

# total inflow per cell -- the "how much signalling is here" null
z = np.load(LI / "cellchatdb2_inflow/data/inflow_scores.npz", allow_pickle=True)
tot_inflow = pd.Series(z["values"].sum(1), index=z["cells"].astype(str))

SOURCES = [
    ("MOFA-Flex (tutorial QC, nzf>0.01)", LI / "mofaflex_inflow/data/factor_scores.csv.gz", "csv"),
    ("MOFA-Flex (sensitivity, nzf>0.001)",
     LI / "mofaflex_inflow/sensitivity_nzf0.001/data/factor_scores.csv.gz", "csv"),
    ("MOFA-Flex (reachability-normalised >0.01, all 9 views)",
     LI / "mofaflex_inflow_reachnorm/data/factor_scores.csv.gz", "csv"),
    ("NMF on inflow", LI / "nmf_inflow/data/NMF_W_factor_scores.csv", "csv_nmf"),
]

report = {}
for label, path, kind in SOURCES:
    if not path.exists():
        print(f"  SKIP (absent): {path}"); continue
    F = pd.read_csv(path, index_col=0)
    if kind == "csv_nmf":                       # NMF scores are written without a usable index
        F = pd.read_csv(path)
        F.index = z["cells"].astype(str)[: len(F)]
    F.index = F.index.astype(str)
    F = F.select_dtypes("number")
    common = U.index.intersection(F.index)
    Ua, Fa = U.loc[common], F.loc[common]

    rho = pd.DataFrame(
        np.vstack([spearmanr(Ua.values[:, i], Fa.values, axis=0).statistic[0, 1:]
                   for i in range(Ua.shape[1])]),
        index=Ua.columns, columns=Fa.columns)
    rho.to_csv(OUT / "data" / f"rho_{label.split()[0].lower()}_{len(Fa.columns)}f.csv")

    best = rho.abs().idxmax(1)
    hub, hub_n = best.value_counts().index[0], int(best.value_counts().iloc[0])
    tot_alz = Ua.sum(1)
    r_act_alz = float(spearmanr(Fa[hub], tot_alz).statistic)
    r_act_inf = float(spearmanr(Fa[hub], tot_inflow.loc[common]).statistic)
    j = rho.loc[MOI].abs().idxmax()

    report[label] = dict(
        n_cells=int(len(common)), n_factors=int(Fa.shape[1]),
        max_abs_rho=float(np.abs(rho.values).max()),
        n_motifs_above_0p3=int((np.abs(rho.values).max(1) > 0.3).sum()),
        hub_factor=hub, n_motifs_mapping_to_hub=hub_n,
        hub_vs_total_alarmist_loading=r_act_alz, hub_vs_total_inflow=r_act_inf,
        moi=MOI, moi_best_factor=j, moi_best_rho=float(rho.loc[MOI, j]),
        moi_vs_total_inflow=float(spearmanr(Ua[MOI], tot_inflow.loc[common]).statistic),
    )
    print(f"\n=== {label}: {len(common):,} cells x {Fa.shape[1]} factors ===")
    print(f"  max |rho| anywhere ............... {report[label]['max_abs_rho']:.3f}")
    print(f"  DEGENERACY: {hub_n} of {Ua.shape[1]} motifs best-match the same factor ({hub})")
    print(f"  ACTIVITY CONFOUND: rho({hub}, total ALARMIST loading) = {r_act_alz:+.3f}")
    print(f"                     rho({hub}, total inflow)          = {r_act_inf:+.3f}")
    print(f"                     rho({MOI}, total inflow)          = "
          f"{report[label]['moi_vs_total_inflow']:+.3f}")
    print(f"  {MOI} best match: {j}  rho={rho.loc[MOI, j]:+.3f}")

    # where each side places the signal -- rank-normalised mean by cell type
    cmp = pd.DataFrame({
        f"ALARMIST {MOI}": Ua[MOI].rank(pct=True).groupby(ct.loc[common]).mean(),
        f"{label} {j}": Fa[j].rank(pct=True).groupby(ct.loc[common]).mean(),
    }).sort_values(f"ALARMIST {MOI}", ascending=False)
    cmp.to_csv(OUT / "data" / f"celltype_placement_{label.split()[0].lower()}.csv")
    print(cmp.round(3).to_string())

json.dump(report, open(OUT / "comparison_summary.json", "w"), indent=2)
print(f"\nwrote {OUT}/comparison_summary.json")
