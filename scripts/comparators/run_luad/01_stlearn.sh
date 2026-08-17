#!/usr/bin/env bash
# ======================================================================================
# stLearn — FOUR independent per-section runs.
#
#   bash scripts/comparators/run_luad/01_stlearn.sh [--dry-run] [SECTION ...]
#
# NO NATIVE MULTI-SAMPLE MODE. stlearn.tl.cci exports exactly seven functions and every one
# takes a single AnnData; there is no condition/sample/batch argument anywhere in the CCI
# API. So the AIS-vs-LUAD comparison for stLearn is made by US, reading four lr_summary
# tables side by side -- stLearn emits no between-condition statistic. That is a capability
# gap and is reported as one.
#
# It must NOT see the concatenated h5ad: st.tl.cci.grid() bins the GLOBAL bounding box, so
# four overlapping coordinate frames would be gridded into one tissue.
#
# --n-col / --n-row come from prep_manifest.json, which derives them from each section's own
# ANNOTATED extent under the tutorial's spot-AREA rule (2,637 um^2 -> 51.35 um square edge;
# see stlearn/DEVIATIONS.md row 11). run_stlearn.py drops unannotated cells before it reads
# the coordinates, so the annotated extent is the one it actually grids.
#
# Cost: ~1.5-4 h per section, ~6-8 GB peak, ~15 GB total output.
# ======================================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../_common/luad_config.sh
source "$HERE/../_common/luad_config.sh"
source "$HERE/_lib.sh"
require_prep

banner "stLearn — 4 per-section runs"

SECS="${*:-}"
SECS="$(echo "$SECS" | tr ' ' '\n' | grep -v '^--' | tr '\n' ' ')"
[ -z "${SECS// /}" ] && SECS="$LUAD_SECTIONS"

ST_LRS="$LUAD_RESULTS_DIR/stlearn/LUAD/cellchatdb2_lrs.txt"

for S in $SECS; do
    NC=$(prep_json "[s for s in m['sections'] if s['section']=='$S'][0]['stlearn_grid']['n_col']")
    NR=$(prep_json "[s for s in m['sections'] if s['section']=='$S'][0]['stlearn_grid']['n_row']")
    step "stLearn $S   (grid ${NC} x ${NR})"
    run "$PY_STLEARN" "$SCRIPTS/stlearn/run_stlearn.py" \
        --h5ad    "$LUAD_PREPPED_DIR/$S.prepped.h5ad" \
        --out-dir "$LUAD_RESULTS_DIR/stlearn/LUAD/$LUAD_TIER/$S" \
        --cell-type-col cell_type \
        --count-layer   counts \
        --lrs           "$ST_LRS" \
        --n-col "$NC" --n-row "$NR" \
        --distance 250 --n-pairs 10000 --n-perms 1000 --min-spots 20 \
        --seed "$LUAD_SEED" \
        --requested-lrs ""
done

# The heavy plotting/export passes replay the saved grid.h5ad; they are cheap and read-only
# with respect to the run.
for S in $SECS; do
    D="$LUAD_RESULTS_DIR/stlearn/LUAD/$LUAD_TIER/$S"
    step "stLearn $S — quant export + full plot pass"
    run "$PY_STLEARN" "$SCRIPTS/stlearn/export_stlearn_quant.py" \
        --run-dir "$D" --label cell_type || true
    run "$PY_STLEARN" "$SCRIPTS/stlearn/plot_stlearn_full.py" \
        --run-dir "$D" --out-dir "$D/plots_full" --requested-lrs "" || true
done

done_banner "stLearn"
cat <<'EOF'

NOTE FOR THE WRITE-UP: run_stlearn.py creates plots/requested/ and leaves it EMPTY,
because --requested-lrs "" was passed deliberately (no LRs of interest were named for this
dataset, and ANXA1 -- half of the GBM pair -- is not on the LUAD panel at all). An empty
requested/ directory is expected output here, not lost output.
EOF
