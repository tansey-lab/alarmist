#!/bin/bash
#SBATCH --job-name=xm_auc
#SBATCH --partition=componc_cpu
#SBATCH --account=tanseyw
#SBATCH --cpus-per-task=4
#SBATCH --mem=160G
#SBATCH --time=08:00:00
#SBATCH --output=/home/fanj2/alarmist/logs/%x_%j.out
# --mem=160G: the peak is LIANA on P17_LUAD -- its inflow npz decompresses to a
# dense 640,739 x 9,313 float32 array (23.9 GB) before it is chunked into CSC.
# P17_AIS measured 10.4 GB peak at 182,378 cells; 160G is ~4x that scaled up.
set -euo pipefail
export TMPDIR=/data1/tanseyw/fanj2/tmp
export PYTHONNOUSERSITE=1
cd /home/fanj2/alarmist
/data1/tanseyw/fanj2/envs/comp-liana/bin/python scripts/comparators/crossmethod/03_run_auc.py
