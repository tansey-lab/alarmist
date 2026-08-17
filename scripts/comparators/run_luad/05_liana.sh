#!/usr/bin/env bash
# ======================================================================================
# LIANA+ — inflow (per section) + ONE joint MOFA-Flex fit. Bivariate and NMF are NOT run.
#
#   bash scripts/comparators/run_luad/05_liana.sh [--dry-run] [--skip-mofaflex]
#
# SCOPE: the user restricted LIANA to inflow + MOFA-Flex, so the bivariate branch, the
# bivariate-morans variant, NMF-on-bivariate, NMF-on-inflow, LR-MISTy and LRIC/cross-PCF are
# deliberately not run here. That is a scope decision, not a failure -- record it as one.
#
# INFLOW IS PER SECTION AND MUST BE. run_inflow.py builds li.ut.spatial_neighbors over the
# WHOLE object; --punch-col is metadata only and does NOT split the graph. At bandwidth
# 13.1454 the kernel support is 28.2 um, and the four sections' coordinate frames coincide
# at distance zero, so a concatenated input would make cells from two different patients
# nearest neighbours. Four invocations on the per-section prepped files.
#
# NO NATIVE MULTI-SAMPLE MODE for the per-pair spatial branches. That is open contract
# deviation CD-2 in METHODS.md, and SKILL.md forbids hand-rolling a substitute. MOFA-Flex is
# the one branch that puts all four sections into a single comparable factor space.
#
# BANDWIDTH stays at 13.1454 um -- the GBM value. The strict ALARMIST analogue for this
# dataset would be 21.0326 (ALARMIST used 80 um patches on these sections; bw = s/3.803630).
# Keeping 13.1454 is the choice that makes LIANA-GBM and LIANA-LUAD comparable, and it keeps
# max_neighbours=100 non-binding (measured max 92 neighbours at r=28.21 um across the four
# sections; at r=45.14 the cap starts truncating 0.43-2.48% of cells). CD-1 -- that this
# number is ALARMIST-derived at all, which violates the do-not-harmonise-kernel-scale rule --
# remains OPEN either way and must be restated in the write-up.
#
# THE MOFA-FLEX COORDINATE TRAP: run_mofaflex.py:433-434 uses a GaussianProcess factor prior
# over obsm[spatial]. Concatenating four overlapping frames would have the GP model four
# tissues as one field. concat_inflow_lrdata.py applies li.ut.expand_coordinates first --
# LIANA's own utility, whose 1.8.1 CHANGELOG names this exact use case -- which translates
# each section onto its own grid cell. Within-section distances are preserved exactly.
# No liana tutorial calls it: log this in liana/DEVIATIONS.md.
# ======================================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../_common/luad_config.sh
source "$HERE/../_common/luad_config.sh"
source "$HERE/_lib.sh"
require_prep

SKIP_MOFA=0
for a in "$@"; do [ "$a" = "--skip-mofaflex" ] && SKIP_MOFA=1; done

BW=13.1454
ROOT="$LUAD_RESULTS_DIR/liana/LUAD"

banner "LIANA+ — inflow x4 + one joint MOFA-Flex"

for S in $LUAD_SECTIONS; do
    step "inflow — $S"
    run "$PY_LIANA" "$SCRIPTS/liana/run_inflow.py" \
        --h5ad    "$LUAD_PREPPED_DIR/$S.prepped.h5ad" \
        --out-dir "$ROOT/$S/${LUAD_TIER}_inflow" \
        --cell-type-col cell_type \
        --count-layer   counts \
        --db "$LUAD_DB" \
        --bandwidth "$BW" --cutoff 0.1 --nz-prop 0.001 --svg-filter \
        --punch-col sample \
        --seed "$LUAD_SEED"

done

# ---------------------------------------------------------------------------------------
# run_inflow_downstream.py is NOT run for LUAD. This is a structural incompatibility, not a
# choice to skip work, and it is verified rather than assumed:
#
#   run_inflow_downstream.py:114-118 does
#       CT, REG = a.cell_type_col, a.region_col
#       for c in (CT, REG): if c not in lrdata.obs: raise SystemExit(f"STOP: obs has no '{c}'")
#   so --region-col cannot be empty, and everything the script adds over run_inflow.py is
#   conditioned on that column being a WITHIN-object contrast:
#       :237-238  compute_global_specificity(groupby=REG) on the non-NaN subset
#       :245-246  group_labels = f"{CT}::{REG}"
#       :147      one spatial embedding per REG level
#
#   On GBM, REG was `grade` and a single object held both grades. Here inflow is per SECTION
#   and a section is exactly one stage -- `stage`, `patient` and `sample` are all constant
#   within the object, so the contrast is empty by construction. Feeding it a constant column
#   would produce a one-level "comparison" that looks like a result and is not one.
#
#   The cross-stage question is answered instead by the joint MOFA-Flex fit below, which is
#   the branch that legitimately spans all four sections.
# ---------------------------------------------------------------------------------------

if [ "$SKIP_MOFA" = "1" ]; then
    step "skipping MOFA-Flex (--skip-mofaflex)"
else
    step "concatenate the four inflow outputs + expand_coordinates"
    IN_DIRS=""
    for S in $LUAD_SECTIONS; do IN_DIRS="$IN_DIRS $ROOT/$S/${LUAD_TIER}_inflow"; done
    # shellcheck disable=SC2086
    run "$PY_LIANA" "$SCRIPTS/liana/concat_inflow_lrdata.py" \
        --in-dirs $IN_DIRS \
        --out "$ROOT/joint_inflow_lrdata.h5ad" \
        --join inner --sample-key sample --spatial-key spatial --n-cols 2 --margin 0.1

    step "MOFA-Flex — ONE joint fit over all four sections"
    # --sample-key sample : the replicate unit for the group-level summary (4 sections).
    # --region-key stage  : run_mofaflex.py asserts exactly 2 levels; AIS/LUAD satisfies it.
    # --drop-obs / --drop-uns default to GBM-specific keys; pass empty lists instead.
    run "$PY_LIANA" "$SCRIPTS/liana/run_mofaflex.py" \
        --lrdata "$ROOT/joint_inflow_lrdata.h5ad" \
        --outdir "$ROOT/mofaflex_inflow_joint" \
        --celltype-key cell_type \
        --region-key   stage \
        --sample-key   sample \
        --spatial-key  spatial \
        --drop-obs --drop-uns \
        --bandwidth "$BW" --cutoff 0.1 \
        --n-factors 20 --max-epochs 1000 --patience 50 \
        --seed "$LUAD_SEED" \
        --lr-of-interest \
        --tag mofaflex_inflow_joint
fi

# ---------------------------------------------------------------------------------------
# annotate_factors.py is NOT run for LUAD, for the same kind of reason. Verified:
#   annotate_factors.py:368-369
#       for branch in a.branches:
#           src = os.path.join(a.results_root, f"nmf_{branch}", "data", "NMF_H_loadings.csv")
# It annotates NMF factor loadings, and NMF is out of scope for this run (inflow +
# MOFA-Flex only). With no nmf_*/ directory there is no input; it would exit on the first
# branch. If factor annotation is wanted later it needs a MOFA-Flex-aware variant.
# ---------------------------------------------------------------------------------------

step "stamp the dataset label onto the manifests"
run "$PY_PREP" "$SCRIPTS/_common/fix_manifest.py" --dataset LUAD --tier "$LUAD_TIER" --find "$ROOT"

done_banner "LIANA+"
cat <<'EOF'

SCOPE NOTE FOR METHODS.md: on GBM, LIANA ran 8 branches. For LUAD only inflow and MOFA-Flex
were run, by user decision. Bivariate, bivariate-morans, NMF (both inputs), LR-MISTy and
LRIC/cross-PCF are NOT missing outputs -- they were not requested. Say so explicitly rather
than leaving a reader to infer a failure.

STATISTICAL HEALTH WARNING: run_mofaflex.py's step [9] aggregates factor scores to
--sample-key and tests them against --region-key. That is n = 2 vs 2. A two-sided rank test
on 2-vs-2 floors at p = 0.333 and can never reach 0.05. Read the factor-vs-stage output as a
direction and a concordance check between P17 and P21, never as a p-value.
EOF
