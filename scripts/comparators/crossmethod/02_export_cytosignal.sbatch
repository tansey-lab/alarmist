#!/bin/bash
#SBATCH --job-name=cs_export
#SBATCH --partition=componc_cpu
#SBATCH --account=tanseyw
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=/home/fanj2/alarmist/logs/%x_%j.out
set -euo pipefail
export TMPDIR=/data1/tanseyw/fanj2/tmp
export R_LIBS_USER=/data1/tanseyw/fanj2/envs/comp-cytosignal/lib/R/library
/data1/tanseyw/fanj2/envs/comp-cytosignal/bin/Rscript \
    /home/fanj2/alarmist/scripts/comparators/crossmethod/02_export_cytosignal.R \
    /data1/tanseyw/projects/fanj2/results/_crossmethod/label_auc/cytosignal_export
