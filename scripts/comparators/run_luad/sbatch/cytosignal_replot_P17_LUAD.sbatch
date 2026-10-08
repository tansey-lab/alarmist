#!/bin/bash
# ======================================================================================
# CytoSignal replot -- P17_LUAD   (536,388 cells post-QC)
#
#   sbatch scripts/comparators/run_luad/sbatch/cytosignal_replot_P17_LUAD.sbatch
#
# Produces the figures run_cytosignal.R skipped, by rebuilding only as far as
# inferIntrScore and injecting the already-banked quant/reslist_*.rds. See
# cytosignal/replot_from_quant.R for why inferSignif does not have to run again.
# Nothing under quant/ is read-modify-written; figures land in <out>/plots/.
#
# WHY --mem=128G AND NOT MORE. The relevant measurement is DEVIATIONS CS-1: P17_LUAD --
# the LARGEST section at 536,388 post-QC cells -- completed findNN + imputeLR +
# inferIntrScore + inferSignif inside a 96 GiB cgroup and was killed only later, inside
# rankIntrSpatialVar -> SPARK::sparkx. CS_MODE=inject stops after inferIntrScore, so it
# never reaches EITHER of the two stages that have ever OOMed. Its peak is therefore
# bounded well under 96 GiB even on the largest section; 128G is headroom, not a guess.
# Largest section, and the one CS-1 measured at 96 GiB through inferSignif.
# The script logs its own peak RSS and records it in plots/replot_manifest.json -- read it
# and cut this further. Check afterwards with:
#   sacct -j <jobid> --format=JobID,State,Elapsed,ReqMem,MaxRSS
# ======================================================================================
#SBATCH --job-name=csreplot_P17_LUAD
#SBATCH --partition=componc_cpu
#SBATCH --account=tanseyw
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=06:00:00
#SBATCH --output=/data1/tanseyw/projects/fanj2/results/_logs/%x_%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=fanj2@mskcc.org

set -euo pipefail
# Slurm COPIES this file to /var/spool/slurmd/job<id>/slurm_script and runs it FROM THERE,
# so $0 and ${BASH_SOURCE[0]} resolve into the spool dir, not the repo. The path below must
# be absolute. `sbatch --test-only` cannot catch a mistake here -- it never executes the
# script. Override with: sbatch --export=ALL,REPO=/other/checkout <this file>
REPO="${REPO:-/home/fanj2/alarmist}"
exec bash "$REPO/scripts/comparators/run_luad/sbatch/_replot_body.sh" P17_LUAD
