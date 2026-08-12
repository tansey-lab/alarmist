#!/usr/bin/env python
"""Why LIANA's inflow + MOFA-Flex cannot recover an ALARMIST-style mGAM motif.

The claim tested here is structural, not a matter of tuning: ALARMIST motif 1 is a BIDIRECTIONAL
loop between two DIFFERENT cell types --- GRN->SORT1 (mGAM to MES-like) and ANXA1->FPR1 (MES-like
to mGAM). LIANA's inflow assigns each score to the RECEIVING cell, so

    mGAM^GRN^SORT1        lands on the cells that RECEIVE GRN   -> tumour states (SORT1+)
    MES-like^ANXA1^FPR1   lands on the cells that RECEIVE ANXA1 -> mGAM (FPR1+)

i.e. the two arms are scored on largely DISJOINT cell populations. A MOFA-Flex factor is a
per-cell score, so for both arms to load on one factor the factor would have to be high on
mGAM and on MES-like cells simultaneously --- but the two feature vectors are near-orthogonal
across cells, so there is no signal to bind them.

ALARMIST factorises a PATCH x LRI matrix instead. A 50 um patch contains both cell types, so
both arms are non-zero in the SAME ROW and the factorisation can bind them.

This script measures exactly that, on the same tissue and the same two interactions, changing
only the unit of analysis. It is read-only.

Run:
    /Users/jiayifan/anaconda3/envs/bptf/bin/python \\
        scripts/comparators/liana/why_no_mgam_motif.py
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import spearmanr

ROOT = Path("/Users/jiayifan/tansey_lab/alarmist")

p = argparse.ArgumentParser()
p.add_argument("--inflow-npz", default=str(ROOT / "results/comparators/liana/GBM/cellchatdb2_inflow/data/inflow_scores.npz"))
p.add_argument("--cell-meta", default=str(ROOT / "results/comparators/liana/GBM/cellchatdb2_inflow/data/cell_meta.csv"))
p.add_argument("--patch-npz", default=str(ROOT / "results/GBM/patch_lri_matrix.npz"))
p.add_argument("--patch-cols", default=str(ROOT / "results/GBM/patch_lri_columns.csv"))
p.add_argument("--arm1-cell", default="mGAM^GRN^SORT1")
p.add_argument("--arm2-cell", default="MES-like^ANXA1^FPR1")
p.add_argument("--arm1-patch", default="mGAM|MES-like|GRN|SORT1")
p.add_argument("--arm2-patch", default="MES-like|mGAM|ANXA1|FPR1")
p.add_argument("--out-dir", default=str(ROOT / "results/comparators/liana/GBM/vs_alarmist"))
a = p.parse_args()
OUT = Path(a.out_dir); OUT.mkdir(parents=True, exist_ok=True)


def stats(x, y, unit):
    nz1, nz2 = x > 0, y > 0
    both = nz1 & nz2
    return dict(
        unit=unit, n=int(len(x)),
        arm1_nonzero=int(nz1.sum()), arm1_pct=float(100 * nz1.mean()),
        arm2_nonzero=int(nz2.sum()), arm2_pct=float(100 * nz2.mean()),
        both_nonzero=int(both.sum()), both_pct=float(100 * both.mean()),
        pearson=float(np.corrcoef(x, y)[0, 1]),
        spearman=float(spearmanr(x, y).statistic),
        p_arm2_given_arm1=float(both.sum() / max(nz1.sum(), 1)),
        p_arm2_marginal=float(nz2.mean()),
    )


# ------------------------------------------------------------------ CELL level (LIANA)
z = np.load(a.inflow_npz, allow_pickle=True)
X, feats, cells = z["values"], z["features"].astype(str), z["cells"].astype(str)
fi = {f: i for i, f in enumerate(feats)}
meta = pd.read_csv(a.cell_meta); meta["cell"] = meta["cell"].astype(str)
ct = meta.set_index("cell").loc[cells, "cell_type"].values
c1, c2 = X[:, fi[a.arm1_cell]], X[:, fi[a.arm2_cell]]
cell = stats(c1, c2, "cell (LIANA inflow)")

print("=== WHICH CELLS receive each arm? (inflow scores the RECEIVER) ===")
comp = {}
for nm, v in [(a.arm1_cell, c1), (a.arm2_cell, c2)]:
    top = pd.Series(ct[v > 0]).value_counts(normalize=True).head(4)
    comp[nm] = {k: round(float(p_), 3) for k, p_ in top.items()}
    print(f"  {nm:24s} -> " + ", ".join(f"{k} {100*p_:.0f}%" for k, p_ in top.items()))

# ----------------------------------------------------------------- PATCH level (ALARMIST)
cols = pd.read_csv(a.patch_cols)["column_name"].values
M = sp.load_npz(a.patch_npz).tocsc()
pi = {c: i for i, c in enumerate(cols)}
g1 = [c for c in cols if c.startswith(a.arm1_patch)]
g2 = [c for c in cols if c.startswith(a.arm2_patch)]
p1 = np.asarray(M[:, [pi[c] for c in g1]].todense()).sum(1)
p2 = np.asarray(M[:, [pi[c] for c in g2]].todense()).sum(1)
patch = stats(p1, p2, "50 um patch (ALARMIST)")

print("\n=== SAME two interactions, only the UNIT of analysis differs ===")
hdr = f"{'':26s} {'cells (LIANA)':>16s} {'patches (ALARMIST)':>20s}"
print(hdr); print("-" * len(hdr))
for k, lab in [("n", "rows in the matrix"), ("both_pct", "rows with BOTH arms (%)"),
               ("pearson", "Pearson r"), ("spearman", "Spearman rho"),
               ("p_arm2_given_arm1", "P(arm2>0 | arm1>0)"), ("p_arm2_marginal", "P(arm2>0) marginal")]:
    f = ",.0f" if k == "n" else ".4f"
    print(f"{lab:26s} {format(cell[k], f):>16s} {format(patch[k], f):>20s}")
enrich_c = cell["p_arm2_given_arm1"] / max(cell["p_arm2_marginal"], 1e-12)
enrich_p = patch["p_arm2_given_arm1"] / max(patch["p_arm2_marginal"], 1e-12)
print(f"{'co-occurrence enrichment':26s} {enrich_c:>15.1f}x {enrich_p:>19.1f}x")
print(f"\n  Pearson r rises {patch['pearson']/max(cell['pearson'],1e-12):.0f}x on the same data, "
      f"purely by aggregating cells into 50 um patches.")
print("  The loop is a property of a NEIGHBOURHOOD, and a neighbourhood is not a row in")
print("  LIANA's matrix. No amount of factorisation can recover a structure the input lacks.")

json.dump({"cell_level": cell, "patch_level": patch,
           "receiver_composition": comp,
           "enrichment_cell": enrich_c, "enrichment_patch": enrich_p},
          open(OUT / "why_no_mgam_motif.json", "w"), indent=2)
print(f"\nwrote {OUT}/why_no_mgam_motif.json")
