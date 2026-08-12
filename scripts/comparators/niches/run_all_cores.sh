#!/usr/bin/env bash
# Drive run_niches.R over every TMA core of the GBM slide, both imputation sub-runs.
#
#   bash scripts/comparators/niches/run_all_cores.sh [--force]
#
# Cores run SMALLEST FIRST (cores.csv is sorted by n_cells) because NICHES' ComputeEdgelist
# is dense O(N^2): core1 at 26,456 cells needs ~4 x N^2 x 8 B ~= 22 GB of R heap on a 36 GB
# box, so it is the one most likely to die -- running it last means a failure there still
# leaves the other 12 cores complete.
#
# Each core writes a DONE sentinel, so re-running this script skips finished work; nothing
# is ever recomputed unless you pass --force. Peak RSS per core is captured with
# /usr/bin/time -l and appended to run_timings.csv.
# NB: no `set -u` -- the conda-forge compiler activation scripts sourced by activate_env.sh
# reference unset variables and abort the shell under it (rc=127).
set -o pipefail

REPO=/Users/jiayifan/tansey_lab/alarmist
cd "$REPO"
source scripts/comparators/niches/activate_env.sh

IN=results/comparators/niches/GBM/input
OUT=results/comparators/niches/GBM/cellchatdb2
LOGS=$OUT/logs
FORCE=""
[[ "${1:-}" == "--force" ]] && FORCE="--force"

mkdir -p "$LOGS"
TIMINGS=$OUT/run_timings.csv
[[ -f $TIMINGS ]] || echo "core,n_cells,grade,impute,status,wall_sec,peak_rss_gb" > "$TIMINGS"

# cores.csv: core,tma_id,n_cells,grade,... (already ascending by n_cells)
tail -n +2 "$IN/cores.csv" | while IFS=, read -r core tma n_cells grade rest; do
  for impute in none alra; do
    tag=$([[ $impute == none ]] && echo noimpute || echo alra)
    odir=$OUT/$tag/$core
    log=$LOGS/${core}_${tag}.log

    if [[ -f $odir/DONE && -z $FORCE ]]; then
      echo "SKIP  $core ($n_cells cells, $grade) [$tag] -- already done"
      continue
    fi

    echo "RUN   $core ($n_cells cells, $grade) [$tag] -> $log"
    start=$(date +%s)
    /usr/bin/time -l Rscript scripts/comparators/niches/run_niches.R \
        --input "$IN/$core" --output "$odir" \
        --lr-db "$IN/custom_LR_database.csv" --impute "$impute" $FORCE \
        > "$log" 2>&1
    rc=$?
    end=$(date +%s)
    wall=$((end - start))

    # /usr/bin/time -l on macOS reports "maximum resident set size" in BYTES
    rss=$(grep "maximum resident set size" "$log" | awk '{printf "%.2f", $1/1024/1024/1024}')
    [[ -z $rss ]] && rss="NA"
    status=$([[ $rc -eq 0 ]] && echo ok || echo "FAIL(rc=$rc)")
    echo "$core,$n_cells,$grade,$tag,$status,$wall,$rss" >> "$TIMINGS"
    echo "      -> $status  ${wall}s  peak ${rss} GB"

    if [[ $rc -ne 0 ]]; then
      echo "      !! failed; see $log (continuing with remaining cores)"
      tail -5 "$log" | sed 's/^/      | /'
    fi
  done
done

echo
echo "===== timings ====="
cat "$TIMINGS"
