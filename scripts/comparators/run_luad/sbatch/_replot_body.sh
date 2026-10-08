#!/usr/bin/env bash
# Shared body for the four cytosignal_replot_<SECTION>.sbatch files. Not submitted directly.
#   bash _replot_body.sh <SECTION>
# Everything host-specific comes from _common/luad_config.sh -- the one file to edit.
set -euo pipefail

SECTION="${1:?usage: _replot_body.sh <SECTION>}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DEFAULT="$(cd "$HERE/../../../.." && pwd)"
: "${REPO:=$REPO_DEFAULT}"

# activate_env.sh sources luad_config.sh itself and pins PATH + R_LIBS_USER to
# comp-cytosignal. Plain `conda activate` is NOT used: it falls through to the system R.
# shellcheck source=../../cytosignal/activate_env.sh
source "$REPO/scripts/comparators/cytosignal/activate_env.sh"

export TMPDIR="${TMPDIR:-/data1/tanseyw/projects/fanj2/tmp}"
mkdir -p "$TMPDIR"
# R's BLAS will otherwise grab every core on the node and inflate RSS; the package's own
# parallelism is numCores=4 inside inferIntrScore, which is unaffected by these.
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"

IN="$LUAD_RESULTS_DIR/cytosignal/LUAD/input/$SECTION"
OUT="$LUAD_RESULTS_DIR/cytosignal/LUAD/$LUAD_TIER/$SECTION"
CS_DB="$LUAD_RESULTS_DIR/cytosignal/cellchat_db_human.rds"

echo "================================================================="
echo " CytoSignal replot -- $SECTION"
echo "   host      $(hostname)   job ${SLURM_JOB_ID:-<none>}"
echo "   cpus      ${SLURM_CPUS_PER_TASK:-?}   mem ${SLURM_MEM_PER_NODE:-?} MB"
echo "   Rscript   $(command -v Rscript)"
echo "   input     $IN"
echo "   out       $OUT"
echo "   db        $CS_DB"
echo "   mode      ${CS_MODE:-inject}"
echo "================================================================="

# Fail before burning an allocation, not an hour in.
[ -d "$IN" ]      || { echo "FATAL: input dir missing: $IN"; exit 1; }
[ -f "$CS_DB" ]   || { echo "FATAL: db missing: $CS_DB"; exit 1; }
[ -d "$OUT/quant" ] || { echo "FATAL: no quant/ under $OUT -- run 03_cytosignal.sh for this section first,"; \
                         echo "       or set CS_MODE=full to compute significance from scratch."; exit 1; }
Rscript -e 'q(status = !requireNamespace("cytosignal", quietly = TRUE))' \
  || { echo "FATAL: comp-cytosignal has no cytosignal package"; exit 1; }
# scattermore backs raster=TRUE; without it plotIntrValue stops on every large section.
Rscript -e 'q(status = !requireNamespace("scattermore", quietly = TRUE))' \
  || { echo "FATAL: scattermore missing -- raster=TRUE plotting will fail. Install it into comp-cytosignal first."; exit 1; }
# circlize is a cytosignal *Suggests* and gates plotCircosNIntr ONLY. Warn, never fail:
# the other six plot types are unaffected. Fix with cytosignal/install_plot_deps.sh.
Rscript -e 'q(status = !requireNamespace("circlize", quietly = TRUE))' \
  || echo "WARNING: circlize missing -> plotCircosNIntr will be skipped. Run: bash $REPO/scripts/comparators/cytosignal/install_plot_deps.sh"

# Validate everything above without spending the allocation on the run itself.
if [ "${CS_PREFLIGHT_ONLY:-0}" = "1" ]; then echo ">> preflight OK for $SECTION (CS_PREFLIGHT_ONLY=1); not running"; exit 0; fi

exec Rscript "$REPO/scripts/comparators/cytosignal/replot_from_quant.R" "$IN" "$OUT" "$CS_DB"
