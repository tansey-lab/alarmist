#!/usr/bin/env bash
# ======================================================================================
# CellChat — TWO objects (one per condition, two sections each), then its NATIVE
# cross-condition comparison.
#
#   bash scripts/comparators/run_luad/04_cellchat.sh [--dry-run] [--workers N]
#
# THE ONLY METHOD THAT LEGITIMATELY PUTS SEVERAL SECTIONS IN ONE OBJECT.
# computeRegionDistance (CellChat/R/modeling.R:1194-1228) loops over meta$samples and calls
# BiocNeighbors::queryKNN strictly WITHIN each sample; nothing ever queries across samples.
# So the four overlapping coordinate frames are inert here -- two sections whose x/y happen
# to coincide are still never treated as neighbours. That is why this one reads the
# CONCATENATED file while stLearn and LIANA must not.
#
#   AIS  object: meta$samples = {P17_AIS,  P21_AIS}    475,240 cells
#   LUAD object: meta$samples = {P17_LUAD, P21_LUAD} 1,200,922 cells
#   -> mergeCellChat -> compareInteractions / rankNet(do.stat=TRUE) /
#      netVisual_diffInteraction / manifold learning / presto-backed cross-condition DEA
#
# All 19 annotation_coarse categories are present in all four sections, so levels() is
# identical across the two objects: liftCellChat is NOT needed and the functional-similarity
# comparison is applicable.
#
# ROW ORDER OF spatial_factors.csv IS LOAD-BEARING. computeRegionDistance indexes ratio[k]
# and tol[k] POSITIONALLY over levels(meta$samples) (modeling.R:1165,1212). prepare_gbm_input.py
# writes one row per sample sorted the way R sorts character factor levels. Do not hand-edit it.
#
# TIER: for CellChat "cellchatdb2" is not a different database -- the BUNDLED CellChatDB.human
# already IS CellChatDB v2 (audit_db_equivalence.R found Jaccard 1.0000 against our export).
# The tier flag selects which annotation categories are used: `default` = subsetDB(search=
# "Secreted Signaling") = 1,280 interactions, `cellchatdb2` = subsetDB(CellChatDB) = 2,239.
# So skipping `default` here drops the Secreted-only VIEW, not a database. Say it that way.
#
# MEMORY: the dominant term is modeling.R:207-215 -- 100 bootstraps of aggregate(t(data.use))
# over a dense nCells x ~733 signaling-gene matrix, run `workers` at a time. One copy is
# 2.6 GB (AIS) / 6.6 GB (LUAD); at workers=4 budget ~21 GB (AIS) and ~52 GB (LUAD). Drop to
# --workers 1 for ~5 / ~13 GB at roughly 4x the wall clock.
# The other risk is Matrix::readMM on the LUAD exchange matrix (~302M nnz, ~6.9 GB of ASCII)
# at run_cellchat.R:85 -- allow ~15 GB transient for it.
# ======================================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../_common/luad_config.sh
source "$HERE/../_common/luad_config.sh"
source "$HERE/_lib.sh"
require_prep

WORKERS=4
for i in $(seq 1 $#); do
    if [ "${!i}" = "--workers" ]; then j=$((i + 1)); WORKERS="${!j}"; fi
done

banner "CellChat — 2 condition objects + native cross-condition comparison"

IN="$LUAD_RESULTS_DIR/cellchat/LUAD/input"
OUT="$LUAD_RESULTS_DIR/cellchat/LUAD/$LUAD_TIER"
run mkdir -p "$IN" "$OUT"

if [ ! -f "$SCRIPTS/cellchat/run_cellchat_luad.R" ]; then
    echo "ERROR: run_cellchat_luad.R missing. Run 00_prep.sh (it calls make_luad_variants.sh)."
    exit 1
fi

step "export AIS/ and LUAD/ exchange directories"
# prepare_gbm_input.py is already fully parameterised (--sample-column/--condition-column)
# and asserts that X is log-normalized, which prep guarantees. It is used UNCHANGED; only
# its flags differ from the GBM invocation. Xenium spatial factors per the CellChat FAQ:
# ratio = 1 (coordinates already in microns), tol = spot.size/2 = 5.
run "$PY_PREP" "$SCRIPTS/cellchat/prepare_gbm_input.py" \
    --h5ad    "$LUAD_PREPPED_DIR/AIS_LUAD_4sections.h5ad" \
    --out-dir "$IN" \
    --sample-column    sample \
    --cell-type-column cell_type \
    --condition-column stage \
    --ratio 1.0 --spot-size 10.0

step "run_cellchat_luad.R  (workers=$WORKERS)"
# future's default globals ceiling is 500 MB; the per-worker payload here is gigabytes.
export R_FUTURE_GLOBALS_MAXSIZE=${R_FUTURE_GLOBALS_MAXSIZE:-16000000000}
echo "    R_FUTURE_GLOBALS_MAXSIZE=$R_FUTURE_GLOBALS_MAXSIZE"
run "$RS_CELLCHAT" "$SCRIPTS/cellchat/run_cellchat_luad.R" \
    --input-dir "$IN" \
    --out-dir   "$OUT" \
    --tier      "$LUAD_TIER" \
    --conditions AIS,LUAD \
    --interaction-range 250 --contact-range 10 --distance-use FALSE \
    --type truncatedMean --trim 0.1 --nboot 100 \
    --seed 1 --workers "$WORKERS" --min-cells 10

step "plot_cellchat_luad.R  (includes the mergeCellChat comparison stage)"
# A FRESH process on purpose: selectK calls NMF::nmfEstimateRank with NMF's parallel foreach
# backend and dies with "All the runs produced an error" inside a session that has already
# done Seurat/ComplexHeatmap work -- while the identical call succeeds in a clean session.
run "$RS_CELLCHAT" "$SCRIPTS/cellchat/plot_cellchat_luad.R" \
    --out-dir "$OUT" \
    --conditions AIS,LUAD \
    --top-n-lr 10

done_banner "CellChat"
cat <<'EOF'

DISK: the GBM CellChat tree is 10 GB, most of it plots/ (9,256 files). LUAD sections carry
8-30x more points per scatter and there are 19 cell types instead of 9. Check the size of
plots/ before starting the next method, and prune the oversized per-pathway spatial SVGs if
needed -- the GBM run had to.

Stage F (the mergeCellChat comparison) only fires when exactly 2 conditions are given
(plot_cellchat.R:356-357). --conditions AIS,LUAD satisfies that. AIS is passed FIRST so
"increased in the second dataset" reads as increased in LUAD.
EOF
