#!/usr/bin/env bash
# Shared helpers for the run_luad/*.sh entry points. Sourced, never executed.
#
#   run <cmd...>   echo the command, then execute it (or only echo when --dry-run)
#   step "text"    a labelled sub-step
#   banner/done_banner
#
# Every script accepts --dry-run, which prints the exact command lines without running
# anything. Use it first, on iris, to confirm the paths resolve.

DRY_RUN=0
for _a in "$@"; do
    case "$_a" in
        --dry-run|-n) DRY_RUN=1 ;;
    esac
done
export DRY_RUN

_c_bold=$'\033[1m'; _c_dim=$'\033[2m'; _c_off=$'\033[0m'
[ -t 1 ] || { _c_bold=""; _c_dim=""; _c_off=""; }

banner() {
    echo
    echo "${_c_bold}==================================================================${_c_off}"
    echo "${_c_bold} $* ${_c_off}"
    [ "$DRY_RUN" = "1" ] && echo "${_c_bold} (DRY RUN — nothing will be executed)${_c_off}"
    echo "${_c_bold}==================================================================${_c_off}"
    echo "  repo     $REPO"
    echo "  prepped  $LUAD_PREPPED_DIR"
    echo "  results  $LUAD_RESULTS_DIR"
    echo "  db       $LUAD_DB"
    echo
}

done_banner() {
    echo
    echo "${_c_bold}--- $* complete ---${_c_off}"
    [ "$DRY_RUN" = "1" ] && echo "(dry run: nothing was actually executed)"
    return 0
}

step() {
    echo
    echo "${_c_bold}>> $*${_c_off}"
}

# Print the command with each argument quoted, then run it (unless dry).
run() {
    local q="" a
    for a in "$@"; do
        case "$a" in
            *[![:alnum:]/._=-]*) q+=" '${a//\'/\'\\\'\'}'" ;;
            *)                   q+=" $a" ;;
        esac
    done
    echo "${_c_dim}   \$${q}${_c_off}"
    if [ "$DRY_RUN" = "1" ]; then return 0; fi
    "$@"
}

# Read a value out of prep_manifest.json (needs prep to have run).
prep_json() {  # prep_json <python expression over `m`>
    "$PY_PREP" -c "
import json,sys
m=json.load(open('$LUAD_PREPPED_DIR/prep_manifest.json'))
print($1)
"
}

require_prep() {
    if [ ! -f "$LUAD_PREPPED_DIR/prep_manifest.json" ]; then
        echo "ERROR: $LUAD_PREPPED_DIR/prep_manifest.json not found."
        echo "       Run  bash scripts/comparators/run_luad/00_prep.sh  first."
        exit 1
    fi
}
