#!/usr/bin/env bash
# LIANA+ `default` tier for GBM -- the authors' OWN LR resource (`consensus`, 4,624 pairs),
# run with every other parameter identical to the `cellchatdb2` tier.
#
# comparator-benchmark/SKILL.md:51-54 requires `default` FIRST and `cellchatdb2` second:
# `default` proves the method was used as its authors recommend, `cellchatdb2` removes the LR
# database as a confounder. Only `cellchatdb2` had been run; this closes that gap.
#
# consensus overlaps CellChatDB v2 by 1,663 of 4,624 pairs (36%), so the two tiers are a real
# comparison, not a formality. GRN^SORT1 and ANXA1^FPR1 are present in BOTH resources.
#
# Directory names mirror the existing tier's convention (cellchatdb2 / cellchatdb2_inflow /
# nmf_bivariate / nmf_inflow) rather than SKILL.md's <dataset>/<tier>/, which this method's
# tree already diverges from -- see DEVIATIONS.md.
#
# Usage: bash scripts/comparators/liana/run_default_tier.sh
set -euo pipefail

PY=/Users/jiayifan/anaconda3/envs/comp-liana/bin/python
S=scripts/comparators/liana
R=results/comparators/liana/GBM
H5AD=data/xenium_mm_final_cell_id.h5ad
L=logs/comparators
mkdir -p "$L"

step () { echo; echo "=== $* ==="; }

step "1/5 bivariate, consensus resource"
$PY $S/run_liana.py --h5ad "$H5AD" --db consensus \
    --out-dir "$R/default" 2>&1 | tee "$L/liana-GBM-default-bivariate.log" | tail -3

step "2/5 inflow, consensus resource"
$PY $S/run_inflow.py --h5ad "$H5AD" --db consensus \
    --out-dir "$R/default_inflow" 2>&1 | tee "$L/liana-GBM-default-inflow.log" | tail -3

step "3/5 inflow downstream (global specificity, proximity, rank_aggregate)"
$PY $S/run_inflow_downstream.py --in-dir "$R/default_inflow" --h5ad "$H5AD" \
    --db consensus 2>&1 | tee "$L/liana-GBM-default-downstream.log" | tail -3

step "4/5 NMF on bivariate"
# --k-max 21 MUST match the cellchatdb2 tier's window. run_nmf.py defaults it to 11, and the
# first pass of this script omitted it -- so the default tier was ranked over k in 1..10 while
# cellchatdb2 was ranked over 1..20. Kneedle normalises over the window it is given, so ranks
# from different windows are NOT comparable and the resulting 6/7-vs-4/4 gap confounded the LR
# database with the search window. Never omit this.
$PY $S/run_nmf.py --k-max 21 --input bivariate --in-dir "$R/default" \
    --out-dir "$R/nmf_bivariate_default" 2>&1 | tee "$L/liana-GBM-default-nmf-biv.log" | tail -3

step "5/5 NMF on inflow"
# --k-max 21 MUST match the cellchatdb2 tier's window. run_nmf.py defaults it to 11, and the
# first pass of this script omitted it -- so the default tier was ranked over k in 1..10 while
# cellchatdb2 was ranked over 1..20. Kneedle normalises over the window it is given, so ranks
# from different windows are NOT comparable and the resulting 6/7-vs-4/4 gap confounded the LR
# database with the search window. Never omit this.
$PY $S/run_nmf.py --k-max 21 --input inflow --in-dir "$R/default_inflow" \
    --out-dir "$R/nmf_inflow_default" 2>&1 | tee "$L/liana-GBM-default-nmf-inf.log" | tail -3

echo; echo "=== default tier complete ==="
