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

# Optional positional SECTION args; any flag (--dry-run, -n) is not a section.
# NB: the obvious `echo "$*" | grep -v '^--'` cannot be used here. Under `set -o pipefail`
# (line 23) a grep that filters EVERY line exits 1, the whole pipeline inherits that, and
# `set -e` then kills the script silently right after the banner -- which is exactly what
# `01_stlearn.sh --dry-run` did (fixed 2026-08-17). A plain loop has no such hazard, and it
# also drops `-n`, which the grep form would have passed through as a section name.
# --skip-run does only the export/plot tail, for sections whose run_stlearn.py already
# finished. Added 2026-08-18: the ArrowStringArray crash killed the script inside the FIRST
# loop, so the second loop never ran for the section that had already completed, and re-running
# it wholesale would redo a ~3.8 h fit and overwrite good output.
SKIP_RUN=0
SECS=""
for _arg in "$@"; do
    case "$_arg" in
        --skip-run) SKIP_RUN=1 ;;
        -*) ;;                          # a flag, already handled by _lib.sh
        *)  SECS="$SECS $_arg" ;;
    esac
done
[ -z "${SECS// /}" ] && SECS="$LUAD_SECTIONS"

ST_LRS="$LUAD_RESULTS_DIR/stlearn/LUAD/cellchatdb2_lrs.txt"

# --n-cpus MUST be passed. stLearn's grid(), run() and run_cci() all do
#     if n_cpus is not None: numba.set_num_threads(n_cpus)
#     else:                  numba.set_num_threads(os.cpu_count())
# (analysis.py:112-115, :280-283, :619-622). os.cpu_count() reports the NODE's core count --
# 56 on isca* -- while numba's ceiling is the cgroup affinity, i.e. --cpus-per-task. So the
# default path raises
#     ValueError: The number of threads must be between 1 and <cpus-per-task>
# the moment gridding starts. This CANNOT be fixed by asking for more CPUs: componc_cpu's
# MaxCPUsPerNode is 52 and the nodes have 56 cores, so os.cpu_count() is unreachable by
# construction. Passing n_cpus is the authors' own documented parameter, not a workaround.
#
# The default is 1, deliberately, for two reasons:
#   1. It is what the GBM runs used. run_stlearn.py:40-46 pins the thread caps to 1 only on
#      Darwin, and those runs were on the Mac with no --n-cpus, so they were single-threaded.
#      The runbook's "1.5-4 h per section" is that single-thread measurement.
#   2. run_cci's permutation p-values are THREAD-COUNT DEPENDENT. het.py:157-200 draws
#      np.random.choice inside a `prange(n_perms)` loop after np.random.seed(seed); under
#      numba each thread carries its own RNG state, so the realised permutations change with
#      the thread count. random_state alone does not pin the result -- n_cpus is part of it.
# Raise it with LUAD_STLEARN_NCPUS=<n> if you want the speed and accept that the run_cci
# p-values will not match a 1-thread run. Keep it the SAME for all four sections either way.
STLEARN_NCPUS="${LUAD_STLEARN_NCPUS:-1}"
_avail=$(nproc)          # affinity-aware, i.e. what the cgroup actually grants
if [ "$STLEARN_NCPUS" -gt "$_avail" ]; then
    echo "ERROR: LUAD_STLEARN_NCPUS=$STLEARN_NCPUS but only $_avail CPUs are visible here."
    echo "       numba would raise 'number of threads must be between 1 and $_avail'."
    echo "       Either lower it or request --cpus-per-task $STLEARN_NCPUS."
    exit 1
fi
echo "  n_cpus   $STLEARN_NCPUS  (of $_avail visible; recorded in each run_manifest.json)"

if [ "$SKIP_RUN" = "1" ]; then
    for S in $SECS; do
        M="$LUAD_RESULTS_DIR/stlearn/LUAD/$LUAD_TIER/$S/run_manifest.json"
        [ -f "$M" ] || { echo "ERROR: --skip-run given but $M is missing -- $S never completed."; exit 1; }
    done
    step "skipping run_stlearn.py (--skip-run); run_manifest.json verified for:$SECS"
fi

for S in $SECS; do
    [ "$SKIP_RUN" = "1" ] && continue
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
        --n-cpus "$STLEARN_NCPUS" \
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
