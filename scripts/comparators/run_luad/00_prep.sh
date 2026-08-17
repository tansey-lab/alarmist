#!/usr/bin/env bash
# ======================================================================================
# STEP 0 — build the shared inputs. Everything else depends on this.
#
#   bash scripts/comparators/run_luad/00_prep.sh [--dry-run]
#
# Produces:
#   $LUAD_PREPPED_DIR/{P17_AIS,P17_LUAD,P21_AIS,P21_LUAD}.prepped.h5ad
#   $LUAD_PREPPED_DIR/AIS_LUAD_4sections.h5ad
#   $LUAD_PREPPED_DIR/prep_manifest.json      <- carries the stLearn grid geometry
#   scripts/comparators/cellchat/{run,plot}_cellchat_luad.R
#
# Cost: ~1 h, peak ~8 GB, ~15 GB of output.
# ======================================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../_common/luad_config.sh
source "$HERE/../_common/luad_config.sh"
source "$HERE/_lib.sh"

banner "STEP 0 — prep"

run mkdir -p "$LUAD_PREPPED_DIR"

step "0a  prepare_luad_input.py — four prepped sections + one concatenated"
run "$PY_PREP" "$SCRIPTS/_common/prepare_luad_input.py" \
    --in-dir  "$LUAD_SRC_DIR" \
    --out-dir "$LUAD_PREPPED_DIR" \
    --cell-type-column annotation_coarse

step "0b  make_luad_variants.sh — generate the two CellChat LUAD variants"
run bash "$SCRIPTS/_common/make_luad_variants.sh"

step "0c  CytoSignal LR database (reuse if already built)"
CS_DB="$LUAD_RESULTS_DIR/cytosignal/cellchat_db_human.rds"
if [ -f "$CS_DB" ] && [ "${DRY_RUN:-0}" != "1" ]; then
    echo "    reusing $CS_DB"
else
    run mkdir -p "$(dirname "$CS_DB")"
    run "$RS_CYTOSIGNAL" "$SCRIPTS/cytosignal/build_cellchat_db.R" "$LUAD_DB" "$CS_DB"
fi

step "0d  LUAD gene panel -> stLearn LR list"
# build_cellchat_lrs.py takes POSITIONAL args: <cellchat_csv> <out_txt> [panel_genes.tsv].
# The panel file matters: the existing GBM list is 1,371 of CellChatDB's 3,218 pairs because
# it was filtered to the GBM 5,119-gene panel. LUAD's panel is a DIFFERENT 5,101 genes, so
# it needs its OWN list -- reusing the GBM one would test the wrong pair set.
PANEL="$LUAD_PREPPED_DIR/luad_panel_genes.tsv"
ST_LRS="$LUAD_RESULTS_DIR/stlearn/LUAD/cellchatdb2_lrs.txt"
run mkdir -p "$(dirname "$ST_LRS")"
run "$PY_PREP" "$SCRIPTS/_common/dump_panel_genes.py" \
    "$LUAD_PREPPED_DIR/P17_AIS.prepped.h5ad" "$PANEL"
run "$PY_STLEARN" "$SCRIPTS/stlearn/build_cellchat_lrs.py" "$LUAD_DB" "$ST_LRS" "$PANEL"

done_banner "STEP 0"
echo
echo "Next: bash scripts/comparators/run_luad/01_stlearn.sh"
