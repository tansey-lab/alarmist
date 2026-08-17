#!/usr/bin/env bash
# ======================================================================================
# Generate the LUAD-specific variants of the two CellChat R scripts.
#
# The user's instruction is: do NOT edit the tracked scripts that produced the signed-off
# GBM outputs. But hand-maintaining a copy of plot_cellchat.R (521 lines) guarantees drift.
# So the copies are GENERATED here by a short list of documented substitutions, and can be
# regenerated at any time. The tracked originals are never modified -- this script asserts
# that at the end.
#
# Exactly THREE substitutions are made, all verified against the current files:
#
#   cellchat/run_cellchat.R:76    REQUESTED_LR <- c("GRN_SORT1", "ANXA1_FPR1")  ->  character(0)
#   cellchat/run_cellchat.R:231   dataset = "GBM"                               ->  dataset = "LUAD"
#   cellchat/plot_cellchat.R:41   REQUESTED_LR <- c("GRN_SORT1", "ANXA1_FPR1")  ->  character(0)
#
# Why REQUESTED_LR is emptied rather than repointed: GRN_SORT1 / ANXA1_FPR1 are the GBM
# motif-1 mGAM<->MES-like loop and mean nothing here, and ANXA1 is not even on this panel
# (verified: absent from all four LUAD/AIS 5,101-gene var indices; GRN, SORT1 and FPR1 are
# present). Leaving it would ship a quant/requested_lr_status.csv reading
# "in_db=TRUE, tested=FALSE" for a gene that cannot be measured. SKILL.md is explicit that
# a dataset's LRs of interest are named by the user and must not be invented, and for this
# run the user named none.
#
# Usage:  bash scripts/comparators/_common/make_luad_variants.sh
# ======================================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CC="$(cd "$HERE/../cellchat" && pwd)"

RUN_SRC="$CC/run_cellchat.R"
RUN_DST="$CC/run_cellchat_luad.R"
PLOT_SRC="$CC/plot_cellchat.R"
PLOT_DST="$CC/plot_cellchat_luad.R"

for f in "$RUN_SRC" "$PLOT_SRC"; do
    [ -f "$f" ] || { echo "ERROR: missing $f"; exit 1; }
done

# Record the inputs so a stale variant is detectable.
SRC_SUM=$(shasum -a 256 "$RUN_SRC" "$PLOT_SRC" | awk '{print $1}' | tr '\n' ' ')

banner() {
    cat <<EOF
# ------------------------------------------------------------------------------------
# GENERATED FILE -- DO NOT EDIT BY HAND.
# Produced by scripts/comparators/_common/make_luad_variants.sh from
#   $1
# Regenerate with: bash scripts/comparators/_common/make_luad_variants.sh
# Substitutions applied: REQUESTED_LR -> character(0); dataset -> "LUAD".
# Source sha256 at generation time: $SRC_SUM
# ------------------------------------------------------------------------------------
EOF
}

# ---- run_cellchat_luad.R
{
    banner "run_cellchat.R"
    sed \
      -e 's/^REQUESTED_LR <- c("GRN_SORT1", "ANXA1_FPR1")$/REQUESTED_LR <- character(0)   # LUAD: no requested LRIs (ANXA1 is not on this panel)/' \
      -e 's/dataset = "GBM"/dataset = "LUAD"/' \
      "$RUN_SRC"
} > "$RUN_DST"

# ---- plot_cellchat_luad.R
{
    banner "plot_cellchat.R"
    sed \
      -e 's/^REQUESTED_LR <- c("GRN_SORT1", "ANXA1_FPR1")$/REQUESTED_LR <- character(0)   # LUAD: no requested LRIs (ANXA1 is not on this panel)/' \
      -e 's|^#   plots/requested_lr/  -- always GRN_SORT1 and ANXA1_FPR1, whatever their rank|#   plots/requested_lr/  -- EMPTY on LUAD: no requested LRIs were named for this dataset|' \
      "$PLOT_SRC"
} > "$PLOT_DST"

# ---- verify the substitutions actually landed (sed silently no-ops on a changed source)
# Guards look at CODE only: a leading-# line is a comment and cannot affect a run.
# NB: no pipelines here. `grep -v ... | grep -q ...` exits the upstream grep with SIGPIPE as
# soon as the downstream one matches, and under `set -o pipefail` that turns a SUCCESSFUL
# match into a non-zero pipeline status.
fail=0
code() { grep -v '^[[:space:]]*#' "$1" || true; }
want() {  # want <file> <regex> <message>
    local body; body="$(code "$1")"
    grep -q -- "$2" <<<"$body" || { echo "ERROR: $3 -- in $1"; fail=1; }
}
reject() {  # reject <file> <regex> <message>
    local body; body="$(code "$1")"
    if grep -q -- "$2" <<<"$body"; then echo "ERROR: $3 -- in $1"; fail=1; fi
}

want   "$RUN_DST"  'REQUESTED_LR <- character(0)' "REQUESTED_LR not substituted"
want   "$RUN_DST"  'dataset = "LUAD"'             "dataset not substituted"
want   "$PLOT_DST" 'REQUESTED_LR <- character(0)' "REQUESTED_LR not substituted"
reject "$RUN_DST"  'GRN_SORT1\|ANXA1_FPR1'        "a requested-LR literal survives in code"
reject "$PLOT_DST" 'GRN_SORT1\|ANXA1_FPR1'        "a requested-LR literal survives in code"
reject "$RUN_DST"  '"GBM"'                        "a GBM dataset label survives in code"

# ---- assert the tracked originals were not touched
NOW_SUM=$(shasum -a 256 "$RUN_SRC" "$PLOT_SRC" | awk '{print $1}' | tr '\n' ' ')
[ "$SRC_SUM" = "$NOW_SUM" ] || { echo "ERROR: a tracked source file changed during generation"; fail=1; }

if [ "$fail" -ne 0 ]; then
    echo "make_luad_variants.sh FAILED -- the source scripts have changed shape; update the sed rules."
    exit 1
fi

echo "wrote:"
echo "  $RUN_DST"
echo "  $PLOT_DST"
echo "tracked originals unchanged (sha256 verified)."
