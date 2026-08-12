#!/usr/bin/env bash
# End-to-end CellChat on GBM: DB audit -> tier `default` -> tier `cellchatdb2`.
# Run from the repo root. Each stage logs to results/comparators/cellchat/GBM/<tier>/{run,plot}.log
#
#   bash scripts/comparators/cellchat/run_all_gbm.sh [default|cellchatdb2|audit|all]
set -euo pipefail
cd "$(dirname "$0")/../../.."
source scripts/comparators/cellchat/activate_env.sh

STAGE="${1:-all}"
IN=results/comparators/cellchat/GBM/input
S=scripts/comparators/cellchat

run_tier () {
  local tier="$1"
  local out="results/comparators/cellchat/GBM/${tier}"
  echo "### tier ${tier} -> ${out}"
  Rscript "$S/run_cellchat.R"  --input-dir "$IN" --out-dir "$out" --tier "$tier" \
                               --conditions low,high
  Rscript "$S/plot_cellchat.R" --out-dir "$out" --conditions low,high
}

case "$STAGE" in
  audit) Rscript "$S/audit_db_equivalence.R" \
           --csv data/LRdatabase/CellChatDBv2.0.human.csv \
           --out-dir results/comparators/cellchat/db_audit ;;
  default|cellchatdb2) run_tier "$STAGE" ;;
  all)
    Rscript "$S/audit_db_equivalence.R" \
      --csv data/LRdatabase/CellChatDBv2.0.human.csv \
      --out-dir results/comparators/cellchat/db_audit
    run_tier default
    run_tier cellchatdb2 ;;
  *) echo "usage: $0 [default|cellchatdb2|audit|all]" >&2; exit 1 ;;
esac
echo "### done: $STAGE"
