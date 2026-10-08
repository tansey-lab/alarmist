#!/usr/bin/env bash
# ======================================================================================
# CytoSignal — four per-section runs, then its NATIVE differential test (mergeCytoSignal +
# runNEBULA) contrasting AIS vs LUAD.
#
#   bash scripts/comparators/run_luad/03_cytosignal.sh [--dry-run] [--skip-significance]
#                                                      [--pairs]
#
# NATIVE MULTI-SAMPLE: yes, and the contrast is a fitted model coefficient, not us comparing
# two output tables. Documented in the authors' own
# "Differential signaling interaction analysis across multiple datasets" tutorial
# (vendored at results/comparators/cytosignal/reference_notebook/).
#
# TWO INDEPENDENT STAGES, and they have very different costs:
#
#   (a) per-section SIGNIFICANCE (run_cytosignal.R). This is the one that is expensive:
#       cytosignal:::permuteLR floors perm.size at ncol(dge.raw), so the null scales with the
#       cell count. Measured 57 GB at 498,422 cells; projected 17 / 28 / 53 / 61 GB for
#       P17_AIS / P21_AIS / P21_LUAD / P17_LUAD. --skip-significance drops this stage if you
#       only want the differential.
#
#   (b) the DIFFERENTIAL (run_nebula_stage.R). mergeCytoSignal never calls inferIntrScore, so
#       it does NOT pay (a)'s permutation cost. It does hold all sections' LRscore matrices at
#       once -- see the dgCMatrix int32 warning in run_nebula_stage.R, and --pairs for the
#       within-patient fallback.
#
# NEEDS A SECOND R: `nebula` does not build inside comp-cytosignal. Set SYS_RSCRIPT and
# SYS_RLIB in _common/luad_config.sh to an R that has it (on iris: the comp-nebula env that
# _common/install_envs_iris.sh creates). run_nebula_stage.R re-invokes itself there.
#
# PLOTS: run_cytosignal.R:69 gates all figures at <= 200,000 cells. Only P17_AIS (182,378)
# is under it, so THREE of the four sections ship quantitative output but no images. That is
# the tracked script's own behaviour and is left alone; it is declared, not silently absorbed.
# ======================================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../_common/luad_config.sh
source "$HERE/../_common/luad_config.sh"
source "$HERE/_lib.sh"
require_prep

SKIP_SIG=0; SKIP_DIFF=0; SKIP_EXPORT=0; PAIRS=""; SECS=""
for a in "$@"; do
    case "$a" in
        --skip-significance) SKIP_SIG=1 ;;
        # Added 2026-08-18, while resuming after the P17_LUAD OOM. The differential is a single
        # call over ALL FOUR sections (run_nebula_stage.R:59 hardcodes the list), so it cannot
        # run until every section has completed -- which is exactly when you are resuming a subset.
        --skip-differential) SKIP_DIFF=1 ;;
        # The four input dirs already exist and check out against their own provenance.json;
        # re-exporting costs ~95 s and rewrites 3.4 GB of counts.mtx for nothing.
        --skip-export)       SKIP_EXPORT=1 ;;
        --pairs)             PAIRS="--pairs" ;;
        -*) ;;                          # any other flag, e.g. --dry-run (handled by _lib.sh)
        *)  SECS="$SECS $a" ;;          # positional SECTION arguments, as 01_stlearn.sh takes
    esac
done
[ -z "${SECS// /}" ] && SECS="$LUAD_SECTIONS"

banner "CytoSignal — 4 per-section runs + native NEBULA differential"

IN="$LUAD_RESULTS_DIR/cytosignal/LUAD/input"
OUT="$LUAD_RESULTS_DIR/cytosignal/LUAD/$LUAD_TIER"
CS_DB="$LUAD_RESULTS_DIR/cytosignal/cellchat_db_human.rds"
run mkdir -p "$IN" "$OUT"

# BUGFIX 2026-08-21. These two lines used to gate only the BANNER: the loop below ran
# unconditionally, so --skip-export printed "skipping" and then exported anyway. It cost no
# correctness -- export_cs_input.py is deterministic, so the bytes were identical -- but it
# rewrote ~5.7 GB of input/*.mtx on every resume and left every input file with an mtime LATER
# than the runs that consumed it, which misreads as "the inputs changed after the analysis".
# (Visible on disk: input/* carry 2026-08-20 12:51 while the section runs finished before that.)
# The sibling flags --skip-significance and --skip-differential were always proper if/else;
# only this one was written with the `[ ] && / [ ] ||` idiom. Now it matches them.
if [ "$SKIP_EXPORT" = "1" ]; then
    step "skipping the CytoSignal input export (--skip-export)"
    # Refuse rather than fail later inside R, after the expensive stages have started.
    for S in $SECS; do
        [ -d "$IN/$S" ] || { echo "ERROR: --skip-export but $IN/$S is missing."; exit 1; }
    done
    echo "    reusing $IN ($(du -sh "$IN" 2>/dev/null | cut -f1))"
else
    step "export the selected sections to the CytoSignal input contract"
    # export_cs_input.py reads layers['counts'], NEVER X. After prep, X is log1p(CP10K) and the
    # old bundle_bignode/export_p21_full.py -- which reads X and casts to int32 -- would silently
    # truncate every count to 0 or 1.
    for S in $SECS; do
        run "$PY_PREP" "$SCRIPTS/cytosignal/export_cs_input.py" \
            "$LUAD_PREPPED_DIR/$S.prepped.h5ad" "$IN/$S"
    done
fi

if [ "$SKIP_SIG" = "0" ]; then
    # Explicit at the call site as well as in luad_config.sh: no LR pair is privileged on
    # LUAD/AIS. GRN_SORT1 / ANXA1_FPR1 belong to the LGG motif-1 story, not to this dataset.
    export CS_EXTRA=""
    for S in $SECS; do
        step "CytoSignal significance — $S   (heavy: perm.size is floored at n_cells)"
        # `nosave` skips the multi-GB cs_checkpoint.rds / cs_result.rds; the quantitative
        # outputs under quant/ are written either way.
        run "$RS_CYTOSIGNAL" "$SCRIPTS/cytosignal/run_cytosignal.R" \
            "$IN/$S" "$OUT/$S" nosave "$CS_DB"
    done
else
    step "skipping the per-section significance stage (--skip-significance)"
fi

if [ "$SKIP_DIFF" = "1" ]; then
    step "skipping the NEBULA differential (--skip-differential) -- it needs ALL FOUR sections"
else
step "NATIVE differential — mergeCytoSignal + runNEBULA on stage (AIS = reference)"
echo "    SYS_RSCRIPT = $SYS_RSCRIPT"
echo "    SYS_RLIB    = $SYS_RLIB"
run "$RS_CYTOSIGNAL" "$SCRIPTS/cytosignal/run_nebula_stage.R" \
    --input-root "$IN" \
    --db         "$CS_DB" \
    --out-dir    "$OUT/differential_stage" \
    $PAIRS
fi

done_banner "CytoSignal"
cat <<'EOF'

IF THE MERGE DIED WITH A LENGTH / ALLOCATION ERROR: that is the dgCMatrix int32 nnz cap
(1,676,162 cells x ~2,683 interactions = ~4.5e9 elements vs a 2,147,483,647 limit). Re-run
this script with --pairs, which merges WITHIN patient instead (P17_AIS+P17_LUAD, then
P21_AIS+P21_LUAD) and reports two matched within-patient contrasts. Given the 2x2 design that
is arguably the better analysis anyway.
EOF
