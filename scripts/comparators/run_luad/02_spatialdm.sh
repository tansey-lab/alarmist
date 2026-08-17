#!/usr/bin/env bash
# ======================================================================================
# SpatialDM — four per-section fits, then its NATIVE cross-condition differential test.
#
#   bash scripts/comparators/run_luad/02_spatialdm.sh [--dry-run]
#
# NATIVE MULTI-SAMPLE: yes. spatialdm.diff_utils (concat_obj -> differential_test ->
# group_differential_pairs) is the authors' own documented workflow -- it is one of the two
# tutorials on the package's front page (tutorial/differential_test_intestine.ipynb, "Differential
# analyses of whole interactome among varying conditions"). It consumes the separately fitted
# per-section objects and RECOMPUTES NOTHING. run_diff_spatialdm.py drives exactly that chain.
#
# WHY FOUR SEPARATE RUNS AND NOT --split-col ON THE CONCATENATED FILE: spatialdm_global returns
# ONE Moran's R per pair per object, so the section must be its own object either way; and
# run_spatialdm.py has no resume and writes per_split_summary.csv only after ALL splits finish,
# so one process over four sections loses everything on a late failure. Four invocations into
# the SAME --out-dir give the layout run_diff_spatialdm.py discovers (one subdir per split),
# with a checkpoint between each. Each invocation OVERWRITES per_split_summary.csv, so we copy
# it aside and concatenate at the end.
#
# MEMORY, MEASURED FROM SOURCE: spatialdm/utils.py:181-183 materialises TWO dense
# 5,101 x N float64 copies (`csr_matrix([norm_max(X) for X in raw_norm.X.T])`), plus L_mat0 /
# R_mat0 and four N x n_pairs arrays. Peak per section:
#     P17_AIS  ~25 GB   P21_AIS  ~40 GB   P21_LUAD ~76 GB   P17_LUAD ~87 GB
# Give this a >=128 GB node. It is the second-hungriest thing in the whole plan.
#
# STATISTICAL HEALTH WARNING: at n_sub = 4 the differential test is ANTI-conservative, not
# merely underpowered. diff_utils.py:98-110 fits a 2-parameter OLS (2 residual df) but scores
# it with `chi2.sf(LR_statistic, 1)` -- df is hardcoded to 1. A clean 2-vs-2 separation like
# z = [0.5, 6.0, 0.6, 5.8] returns p = 1.1e-7. Read the output as exploratory ranking, never
# as a test.
# ======================================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../_common/luad_config.sh
source "$HERE/../_common/luad_config.sh"
source "$HERE/_lib.sh"
require_prep

banner "SpatialDM — 4 per-section fits + native differential"

OUT="$LUAD_RESULTS_DIR/spatialdm/LUAD/$LUAD_TIER"
run mkdir -p "$OUT"

for S in $LUAD_SECTIONS; do
    step "SpatialDM fit — $S"
    run "$PY_SPATIALDM" "$SCRIPTS/spatialdm/run_spatialdm.py" \
        --h5ad    "$LUAD_PREPPED_DIR/$S.prepped.h5ad" \
        --out-dir "$OUT" \
        --split-col sample \
        --cell-type-col cell_type \
        --count-layer   counts \
        --db "$LUAD_DB" \
        --l 75 --cutoff 0.2 --single-cell --min-cell 3 --n-perm 1000 --nproc 1 \
        --n-neighbors auto --n-nearest-neighbors 6 \
        --seed "$LUAD_SEED" \
        --requested-lrs ""
    # each invocation overwrites this file; keep a per-section copy before the next one
    if [ "${DRY_RUN:-0}" != "1" ] && [ -f "$OUT/per_split_summary.csv" ]; then
        cp "$OUT/per_split_summary.csv" "$OUT/per_split_summary_$S.csv"
    fi
done

step "concatenate the four per-section summaries"
run "$PY_PREP" "$SCRIPTS/_common/merge_split_summaries.py" "$OUT"

step "NATIVE differential: concat_obj -> differential_test -> group_differential_pairs"
# --splits is passed explicitly and sorted: run_diff_spatialdm.py:145-147 sorts discovered
# split dirs with `int(x) if x.isdigit() else 0`, which returns 0 for every non-numeric name
# and degenerates to os.listdir order. Correctness is unaffected but column order in
# p_df / zscore_df / tf_df would not be reproducible.
run "$PY_SPATIALDM" "$SCRIPTS/spatialdm/run_diff_spatialdm.py" \
    --run-dir "$OUT" \
    --out-dir "$OUT/differential_stage" \
    --condition-col stage --c1 LUAD --c2 AIS \
    --splits "P17_AIS,P17_LUAD,P21_AIS,P21_LUAD" \
    --requested-lrs "" \
    --seed "$LUAD_SEED"

step "stamp the dataset label onto the manifests"
# run_diff_spatialdm.py:357 hardcodes "dataset": "GBM". The tracked script is not edited;
# the label is corrected afterwards and the original value is preserved in the manifest.
run "$PY_PREP" "$SCRIPTS/_common/fix_manifest.py" --dataset LUAD --find "$OUT"

done_banner "SpatialDM"
cat <<'EOF'

NOT PRODUCED, DECLARE IT: spatialdm/compare_tiers.py needs both a `default` and a
`cellchatdb2` tree. Only cellchatdb2 is being run for LUAD (user decision), so there is no
LUAD tier_comparison/ and there will not be one.
EOF
