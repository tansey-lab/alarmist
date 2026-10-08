#!/usr/bin/env bash
# ======================================================================================
# Submit the CytoSignal replot jobs.
#
#   bash scripts/comparators/run_luad/sbatch/submit_all.sh            # the 3 with NO figures
#   bash scripts/comparators/run_luad/sbatch/submit_all.sh --all      # + P17_AIS
#   bash scripts/comparators/run_luad/sbatch/submit_all.sh --check    # preflight only, no submit
#   bash scripts/comparators/run_luad/sbatch/submit_all.sh P21_AIS    # named sections only
#
# The four .sbatch files are independent -- `sbatch <file>` works on any of them on its
# own. This just submits them together and prints the job ids.
#
# Default set is P17_LUAD, P21_AIS, P21_LUAD: the three sections that finished
# run_cytosignal.R but produced no images because ncol(cs@counts) > 200,000.
# P17_AIS (182,378 cells) already has cluster_map.png + signif_plots/ from the inline
# block; --all adds the wider surface (edge / intrValue / circos / sigCluster) it never got.
# ======================================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

MISSING=(P17_LUAD P21_AIS P21_LUAD)
CHECK=0; SECS=()
for a in "$@"; do
    case "$a" in
        --all)   SECS=(P17_LUAD P21_AIS P21_LUAD P17_AIS) ;;
        --check) CHECK=1 ;;
        -*)      echo "unknown flag: $a" >&2; exit 2 ;;
        *)       SECS+=("$a") ;;
    esac
done
[ ${#SECS[@]} -eq 0 ] && SECS=("${MISSING[@]}")

if [ "$CHECK" = "1" ]; then
    for S in "${SECS[@]}"; do
        echo "== preflight $S =="
        CS_PREFLIGHT_ONLY=1 bash "$HERE/_replot_body.sh" "$S"
    done
    exit 0
fi

for S in "${SECS[@]}"; do
    f="$HERE/cytosignal_replot_${S}.sbatch"
    [ -f "$f" ] || { echo "no sbatch file for $S ($f)" >&2; exit 1; }
    id=$(sbatch --parsable "$f")
    printf "  submitted %-10s job %s   (%s)\n" "$S" "$id" "$(basename "$f")"
done

cat <<'EOF'

  squeue -u $USER -o "%.10i %.16j %.10M %.10L %.10m %.20R"
  tail -f /data1/tanseyw/projects/fanj2/results/_logs/csreplot_<SECTION>_<jobid>.out

  When a job ends, size the next one from what it actually used:
  sacct -j <jobid> --format=JobID,JobName%18,State,Elapsed,ReqMem,MaxRSS
  (the R script also logs its own peak RSS, and records it in plots/replot_manifest.json)
EOF
