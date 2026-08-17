#!/usr/bin/env bash
# ======================================================================================
# LUAD/AIS comparator benchmark — THE ONLY FILE YOU EDIT WHEN MOVING MACHINES.
#
# Every run_luad_*.sh sources this. Change the paths below for iris (or any host) and
# nothing else needs to move.
#
#   source scripts/comparators/_common/luad_config.sh
#
# Edit the four HOST PATHS, then run `bash scripts/comparators/_common/luad_config.sh`
# on its own to self-check that every path and interpreter actually resolves.
# ======================================================================================

# ---------------------------------------------------------------- HOST PATHS -- EDIT
# Repo checkout root (contains scripts/, data/, results/).
: "${REPO:=/Users/jiayifan/tansey_lab/alarmist}"

# Where the FOUR SOURCE h5ads live. On iris this is different from the repo's data/.
# Must contain: P17_AIS_Xenium.h5ad P17_LUAD_Xenium.h5ad P21_AIS_Xenium.h5ad P21_LUAD_Xenium.h5ad
: "${LUAD_SRC_DIR:=$REPO/data/linghua}"

# Where prepare_luad_input.py writes its prepped h5ads (needs ~15 GB).
: "${LUAD_PREPPED_DIR:=$LUAD_SRC_DIR/prepped}"

# Where all comparator outputs go (needs ~250-300 GB).
: "${LUAD_RESULTS_DIR:=$REPO/results/comparators}"

# Conda envs root — the directory that CONTAINS comp-stlearn/, comp-spatialdm/, ...
: "${CONDA_ENVS:=/Users/jiayifan/anaconda3/envs}"

# ---------------------------------------------------------------- DERIVED -- rarely edit
: "${LUAD_DB:=$REPO/data/LRdatabase/CellChatDBv2.0.human.csv}"
: "${LUAD_TIER:=cellchatdb2}"
: "${LUAD_SEED:=0}"
: "${SCRIPTS:=$REPO/scripts/comparators}"

# The four sections, in the canonical order used by ALARMIST's own row ordering
# (results/AIS_LUAD/single_cell/sample_info.csv). Do not reorder: obs['alarmist_row']
# depends on it.
: "${LUAD_SECTIONS:=P17_AIS P17_LUAD P21_AIS P21_LUAD}"

# Per-method interpreters. One env per method is a hard rule (their numpy/scanpy/Seurat
# pins conflict) -- see scripts/comparators/METHODS.md.
: "${PY_PREP:=$CONDA_ENVS/bptf/bin/python}"
: "${PY_STLEARN:=$CONDA_ENVS/comp-stlearn/bin/python}"
: "${PY_SPATIALDM:=$CONDA_ENVS/comp-spatialdm/bin/python}"
: "${PY_LIANA:=$CONDA_ENVS/comp-liana/bin/python}"
: "${RS_CELLCHAT:=$CONDA_ENVS/comp-cellchat/bin/Rscript}"
: "${RS_CYTOSIGNAL:=$CONDA_ENVS/comp-cytosignal/bin/Rscript}"

# CytoSignal's differential stage needs a SECOND R that has the `nebula` CRAN package,
# because nebula does not build inside comp-cytosignal (see METHODS.md gotchas).
# run_nebula_stage.R reads these two env vars; set them to the comp-nebula env on iris.
: "${SYS_RSCRIPT:=$CONDA_ENVS/comp-nebula/bin/Rscript}"
: "${SYS_RLIB:=$CONDA_ENVS/comp-nebula/lib/R/library}"

export REPO LUAD_SRC_DIR LUAD_PREPPED_DIR LUAD_RESULTS_DIR CONDA_ENVS LUAD_DB LUAD_TIER \
       LUAD_SEED SCRIPTS LUAD_SECTIONS PY_PREP PY_STLEARN PY_SPATIALDM PY_LIANA \
       RS_CELLCHAT RS_CYTOSIGNAL SYS_RSCRIPT SYS_RLIB

# ---------------------------------------------------------------- self-check
luad_check() {
    local bad=0 p
    echo "== paths =="
    for v in REPO LUAD_SRC_DIR LUAD_PREPPED_DIR LUAD_RESULTS_DIR CONDA_ENVS SCRIPTS; do
        p="${!v}"
        if [ -d "$p" ]; then printf "  OK    %-18s %s\n" "$v" "$p"
        else printf "  MISS  %-18s %s\n" "$v" "$p"; bad=1; fi
    done
    printf "  %-5s %-18s %s\n" "$([ -f "$LUAD_DB" ] && echo OK || echo MISS)" LUAD_DB "$LUAD_DB"
    [ -f "$LUAD_DB" ] || bad=1

    echo "== source h5ads =="
    for s in $LUAD_SECTIONS; do
        p="$LUAD_SRC_DIR/${s}_Xenium.h5ad"
        if [ -f "$p" ]; then printf "  OK    %-10s %s\n" "$s" "$p"
        else printf "  MISS  %-10s %s\n" "$s" "$p"; bad=1; fi
    done

    echo "== interpreters =="
    for v in PY_PREP PY_STLEARN PY_SPATIALDM PY_LIANA RS_CELLCHAT RS_CYTOSIGNAL SYS_RSCRIPT; do
        p="${!v}"
        if [ -x "$p" ]; then printf "  OK    %-14s %s\n" "$v" "$p"
        else printf "  MISS  %-14s %s\n" "$v" "$p"; bad=1; fi
    done

    if [ "$bad" -eq 0 ]; then echo "ALL OK"; else echo "SOMETHING IS MISSING -- fix luad_config.sh before running anything"; fi
    return $bad
}

# Executed directly (not sourced) -> run the self-check.
if [ "${BASH_SOURCE[0]}" = "${0}" ]; then luad_check; fi
