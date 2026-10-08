#!/usr/bin/env bash
# ======================================================================================
# CellChat — TWO objects (one per condition, two sections each), then its NATIVE
# cross-condition comparison.
#
#   bash scripts/comparators/run_luad/04_cellchat.sh [--dry-run] [--workers N] [--plan P]
#        [--conditions AIS,LUAD] [--reuse-existing] [--skip-export] [--skip-plots]
#
#   RESUME (2026-08-20): AIS is done and LUAD is not. Re-run just LUAD, keep AIS, and still
#   get a two-condition manifest and the mergeCellChat comparison:
#       bash scripts/comparators/run_luad/04_cellchat.sh --skip-export --reuse-existing ...
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
# MEMORY -- TWO ESTIMATES HAVE NOW BEEN BLOWN HERE. The numbers below are the third pass and
# the first one built on a measured ledger rather than on the size of `data.use`. 2026-08-20.
#
#   run 1  --mem=192G, multisession, workers=4  -> OOM at MaxRSS 192.0 GiB, in the AIS condition
#   run 2  --mem=192G, multicore,    workers=4  -> AIS finished in 2166.7 s; OOM at MaxRSS
#                                                  201,326,424K = 191.99989 GiB in LUAD
# Both are cgroup OOM kills (ConstrainRAMSpace=yes, no OverMemoryKill), so the process is
# capped and shot -- MaxRSS reports the CAP, never the requirement. Do not size off it.
#
# WHERE IT ACTUALLY GOES. computeCommunProb calls my.sapply TWICE, and the earlier call is the
# expensive one, which the previous note missed:
#   deparse line 160 (ONCE per condition):  data.use.avg.boot <- my.sapply(1:nboot, function(nE)
#                                             aggregate(t(data.use), list(group[permutation[,nE]]),
#                                                       FUN = trimmed mean))
#   deparse line 214 (once per LR pair):    Pboot <- my.sapply(1:nboot, ...) on 733 x 19 matrices
# Line 214 is cheap and independent of nCells. Line 160 runs `nboot`=100 full aggregations of
# the whole dense matrix, and it is where run 2 died -- the ">>> Run CellChat on spatial <<<"
# print is deparse line 100, i.e. immediately before it.
#
# MEASURED (this repo, 2026-08-20): one line-160 iteration costs 5.9e-5 s per cell and peaks at
# 4.6-4.8x the dense matrix, of which ~3.6x is allocated by the worker (t() plus the
# as.data.frame/split copies inside stats:::aggregate.data.frame). Dense `data.use` is
# nCells x 733 x 8 B = 2.60 GiB (AIS) and 6.56 GiB (LUAD).
#
# LEDGER for the LUAD condition, at deparse line 160 (GiB, ESTIMATED unless marked):
#   parent   sparse @data ~4 + data 6.56 (MEASURED, from the coercion warning)
#            + data.use 6.56 + permutation 0.45 + coords/meta ~0.3      =  ~18
#   per fork t(data.use) + aggregate internals = 3.6 x 6.56 = 23.6 (from the MEASURED ratio)
#            + parent pages privatised when the child's GC mark phase dirties them ~15
#                                                                       =  ~40 each
#   workers=4 -> 18 + 4 x ~43 = ~190, bracketing the observed 192.0.
#
# READ THE ANCHOR CORRECTLY. 192.0 GiB is LEFT-CENSORED: the cgroup killer fired while the four
# children were still allocating, so the observation says "> 192", never "= 192". It bounds the
# true per-child private demand from BELOW at (192 - 18)/4 = 43.5 GiB and not at all from above.
# Every band here is therefore a floor with a soft ceiling. Three independent passes measured
# the deparse-160 transient at 3.96-4.09x, 4.19x and 4.6-4.8x of the dense matrix, and the
# per-cell cost at 5.9e-5, 6.9e-5 and 7.3e-5 s (all still sub-linear at 50-100k, i.e. the small
# end over-states it). Quote the RANGE; none of the three is settled.
#
# BUDGET (ESTIMATED peak RSS, LUAD condition):
#   --workers 1                   46-70 GiB     no fork at all; my.sapply is plain sapply
#   --workers 2 --plan multicore  >=106 GiB, no upper bound  -- 256G would NOT be safe here
#   --workers 4 --plan multicore  >=192 GiB     -- MEASURED to exceed it; budget 448G
#
# WALL CLOCK (LUAD condition; MEASURED stage times + the 5.9e-5..7.3e-5 s/cell range):
#   read mtx 2 min (M) | KNN 1.6 min (M) | identifyOverExpressed* ~5 min | deparse 16-100
#   45.8 min (M, workers-INDEPENDENT: computeRegionDistance is called unconditionally at
#   deparse 76) | deparse 160 = 100 x (5.9e-5..7.3e-5) x 1,200,922 = 118-146 min at workers=1,
#   ~30-37 min at workers=4 | LR loop 1-6 min (FASTER at workers=1: it drops 2,856 forks)
#   | tail ~6 min
#   -> --workers 1  3.0-3.5 h     --workers 4  ~1.4 h
#   AIS at workers=1 is ~70-85 min; with --reuse-existing it is 41 s (MEASURED).
#
# RECOMMENDED: --workers 1 with --mem=256G.
#   workers=1 makes CellChat take the path it takes on any single-core machine --
#   my.sapply <- sapply, `future` is never entered -- so there is no fork, no serialisation,
#   and no UNRELIABLE VALUE warning. It removes the entire failure class that has now cost two
#   runs, for about 1.5 extra hours.
#   256G, not 128G, because the error bar on the peak is ONE-SIDED in the dangerous direction
#   (see the censoring note above) and this quantity has already been mis-sized twice.
#   MaxMemPerNode is 975 GiB, so the margin costs queue time and nothing else.
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
# multicore, not the tutorial's multisession -- see the MEMORY block above and the C2 comment
# in cellchat/run_cellchat.R. Override with `--plan multisession` to get the tutorial's exact
# backend back, or pair `--workers 1` with anything to bypass `future` altogether.
PLAN=multicore
# Resume switches, added 2026-08-20 after the LUAD condition was OOM-killed with AIS already
# complete on disk. All default to the previous behaviour, so a bare invocation is unchanged.
CONDITIONS=AIS,LUAD      # --conditions: run a subset. run_cellchat.R has no skip logic.
REUSE=0                  # --reuse-existing: load objects/<cond>.rds for finished conditions
SKIP_EXPORT=0            # --skip-export: input/ is 8.7 GB and already checked out
SKIP_PLOTS=0             # --skip-plots: plot_cellchat_luad.R reads the RDS, so it runs later
for i in $(seq 1 $#); do
    case "${!i}" in
        --workers)        j=$((i + 1)); WORKERS="${!j}"    ;;
        --plan)           j=$((i + 1)); PLAN="${!j}"       ;;
        --conditions)     j=$((i + 1)); CONDITIONS="${!j}" ;;
        --reuse-existing) REUSE=1       ;;
        --skip-export)    SKIP_EXPORT=1 ;;
        --skip-plots)     SKIP_PLOTS=1  ;;
    esac
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
if [ "$SKIP_EXPORT" = "1" ]; then
    # Refuse rather than silently continue: a missing exchange dir would fail much later,
    # inside R, after the expensive stages had already started.
    for _c in ${CONDITIONS//,/ }; do
        [ -f "$IN/$_c/data.mtx" ] || { echo "ERROR: --skip-export but $IN/$_c/data.mtx is missing."; exit 1; }
    done
    echo "    --skip-export: reusing $IN ($(du -sh "$IN" 2>/dev/null | cut -f1))"
else
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
fi

step "run_cellchat_luad.R  (conditions=$CONDITIONS, workers=$WORKERS, plan=$PLAN, reuse=$REUSE)"
# future's default globals ceiling is 500 MB; the per-worker payload here is gigabytes.
export R_FUTURE_GLOBALS_MAXSIZE=${R_FUTURE_GLOBALS_MAXSIZE:-16000000000}
# Forking a process that already holds live OpenMP/BLAS threads can deadlock. CellChat's inner
# loop is Matrix::crossprod, i.e. BLAS. Pin to one thread before R starts -- it also stops the
# `workers` forks from each spawning a full BLAS thread pool inside an 8-CPU cgroup.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
echo "    R_FUTURE_GLOBALS_MAXSIZE=$R_FUTURE_GLOBALS_MAXSIZE  OMP_NUM_THREADS=1"
run "$RS_CELLCHAT" "$SCRIPTS/cellchat/run_cellchat_luad.R" \
    --input-dir "$IN" \
    --out-dir   "$OUT" \
    --tier      "$LUAD_TIER" \
    --conditions "$CONDITIONS" \
    --interaction-range 250 --contact-range 10 --distance-use FALSE \
    --type truncatedMean --trim 0.1 --nboot 100 \
    --seed 1 --workers "$WORKERS" --min-cells 10 --plan "$PLAN" \
    --reuse-existing "$([ "$REUSE" = 1 ] && echo TRUE || echo FALSE)"

step "plot_cellchat_luad.R  (includes the mergeCellChat comparison stage)"
# A FRESH process on purpose: selectK calls NMF::nmfEstimateRank with NMF's parallel foreach
# backend and dies with "All the runs produced an error" inside a session that has already
# done Seurat/ComplexHeatmap work -- while the identical call succeeds in a clean session.
if [ "$SKIP_PLOTS" = "1" ]; then
    echo "    --skip-plots: skipping plot_cellchat_luad.R (it only reads $OUT/objects/*.rds,"
    echo "                  so it can be run on its own once every condition is on disk)"
else
run "$RS_CELLCHAT" "$SCRIPTS/cellchat/plot_cellchat_luad.R" \
    --out-dir "$OUT" \
    --conditions AIS,LUAD \
    --top-n-lr 10
fi

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
