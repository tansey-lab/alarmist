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
: "${REPO:=/home/fanj2/alarmist}"

# Where the FOUR SOURCE h5ads live. On iris this is different from the repo's data/.
# Must contain: P17_AIS_Xenium.h5ad P17_LUAD_Xenium.h5ad P21_AIS_Xenium.h5ad P21_LUAD_Xenium.h5ad
#
# NB (iris): the source h5ads are NOT named this way on disk. They are four identically
# named files in four per-section directories:
#   /data1/tanseyw/projects/spatial_data/linghua/j.ccell.2025.10.004/<SECTION>_Xenium/xenium_mm.h5ad
# prepare_luad_input.py:233 builds "<in-dir>/<SECTION>_Xenium.h5ad", which no single
# directory on iris satisfies, so LUAD_SRC_DIR points at a SYMLINK FARM of that shape.
# The symlinks are read-only pointers; the originals are never written to.
: "${LUAD_SRC_DIR:=/data1/tanseyw/projects/fanj2/luad_src}"

# Where prepare_luad_input.py writes its prepped h5ads (needs ~15 GB).
: "${LUAD_PREPPED_DIR:=/data1/tanseyw/projects/fanj2/prepped}"

# Where all comparator outputs go (needs ~250-300 GB).
# NOT under $HOME: home is a 100 GB quota; /data1/tanseyw is a 100 TB WekaFS volume with
# ~15 TB free (verified 2026-08-17 -- and note /data1 itself is node-local ext4, so a
# `df -h /data1` reading of ~67 GB is the node's root disk, NOT this storage).
: "${LUAD_RESULTS_DIR:=/data1/tanseyw/projects/fanj2/results}"

# ---------------------------------------------------------------------------------------
# CONDA_ENVS DELIBERATELY USES THE OLD ROOT. This is the one path that is not
# /data1/tanseyw/projects/fanj2/..., and it is not an oversight.
#
# On 2026-08-17 the project directory was moved:
#     mv /data1/tanseyw/fanj2 /data1/tanseyw/projects/fanj2
#     ln -s /data1/tanseyw/projects/fanj2 /data1/tanseyw/fanj2
# so /data1/tanseyw/fanj2 is now a SYMLINK to the real location.
#
# Conda environments are not relocatable: the absolute prefix is baked into entry-point
# shebangs in envs/*/bin/*, into R's etc/Makeconf CC/CXX/FC, into .pc files, and into
# conda-meta. The six envs here were expensive to build (comp-nebula took three attempts;
# comp-liana needed a manual mofaflex checkout at a specific SHA), so they were left in
# place and the symlink keeps their baked-in prefix valid.
#
# *** DO NOT DELETE /data1/tanseyw/fanj2. *** It is load-bearing. Removing it breaks all
# six environments. The data directories above have no such constraint -- nothing writes
# an absolute path into them -- which is why they point at the real location and this
# does not. Same for envs_dirs/pkgs_dirs in ~/.condarc.
# ---------------------------------------------------------------------------------------
: "${CONDA_ENVS:=/data1/tanseyw/fanj2/envs}"

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
# PY_PREP: DEVIATION from the runbook, which names the ALARMIST driver env `bptf`.
# No `bptf` recipe or lock exists anywhere in this tree (install_envs_iris.sh:144 explicitly
# refuses to build it), and the laptop it lived on is unavailable. Rather than guess at it,
# the six PY_PREP entry points were AST-scanned for their imports:
#   _common/prepare_luad_input.py  _common/dump_panel_genes.py  _common/merge_split_summaries.py
#   _common/fix_manifest.py        cytosignal/export_cs_input.py  cellchat/prepare_gbm_input.py
# The union is exactly four third-party packages -- anndata, numpy, pandas, scipy -- with no
# scanpy, no h5py import of their own (anndata pulls it), no subprocess/importlib/exec, and
# CRUCIALLY no import of the `alarmist` package itself. So PY_PREP does not need the ALARMIST
# runtime at all; it needs a scientific-stack env. `comp-prep` is that env, built as
#   conda create -n comp-prep -c conda-forge python=3.11 numpy anndata pandas scipy h5py
# Resolved: py 3.11.15, numpy 2.4.6, pandas 2.3.3, scipy 1.17.1, anndata 0.12.19, h5py 3.16.0.
#
# NUMPY IS DELIBERATELY *NOT* PINNED <= 1.24, despite CLAUDE.md. That pin's stated reason is
# numba ("Numba needs NumPy 1.24 or less") and applies to the ALARMIST `bptf` runtime.
# comp-prep has no numba and no llvmlite (verified), so the constraint does not apply here.
# It was applied on the first build and stage 00 died immediately at
# prepare_luad_input.py:159:
#     adata.obs_names = pd.Index(f"{section}_" + base.to_numpy().astype(str))
#     numpy.core._exceptions._UFuncNoLoopError: ufunc 'add' did not contain a loop with
#     signature matching types (dtype('<U8'), dtype('<U10'))
# `str + ndarray('<U')` needs the string ufuncs NumPy 2.0 added; under 1.24 it cannot work.
# That line is the ONLY such construct in the six PY_PREP entry points, and none of them uses
# a name NumPy 2 removed (np.float_, np.NaN, np.in1d, ...), so numpy 2.x is safe in both
# directions. It also matches the other four comparator envs, which are all numpy 2.x.
# NB this means the runbook's nominal `bptf` (numpy <= 1.24) could never have run this script.
: "${PY_PREP:=$CONDA_ENVS/comp-prep/bin/python}"
: "${PY_STLEARN:=$CONDA_ENVS/comp-stlearn/bin/python}"
: "${PY_SPATIALDM:=$CONDA_ENVS/comp-spatialdm/bin/python}"
: "${PY_LIANA:=$CONDA_ENVS/comp-liana/bin/python}"
: "${RS_CELLCHAT:=$CONDA_ENVS/comp-cellchat/bin/Rscript}"
: "${RS_CYTOSIGNAL:=$CONDA_ENVS/comp-cytosignal/bin/Rscript}"

# CytoSignal's differential stage needs a SECOND R that has the `nebula` CRAN package,
# because nebula does not build inside comp-cytosignal (see METHODS.md gotchas).
# run_nebula_stage.R reads these two env vars; set them to the comp-nebula env on iris.
# DEVIATION (iris, 2026-08-17): comp-nebula is built on R 4.4, not the r-base=4.3 that
# install_envs_iris.sh:82 asks for. Three reasons, in order of weight:
#   1. R 4.4 is what the SOURCE machine used for this second R. run_nebula_stage.R:56 defaults
#      SYS_RLIB to ~/Library/R/arm64/4.4/library and cytosignal/activate_env.sh describes the
#      fall-through R as 4.4.2. The installer's 4.3 is a copy-paste from the comp-cytosignal
#      recipe, not a considered pin.
#   2. CRAN's nebula (1.5.8) declares `Depends: R (>= 4.4.0)`, so on r-base=4.3 it is simply
#      "not available for this version of R" -- the installer's own fallback cannot work.
#   3. conda-forge ships no r43 build of r-seurat, and bioconda's SingleCellExperiment is
#      uninstallable here (its genomeinfodbdata post-link script is broken).
# Built: R 4.4.3, nebula 1.5.8, Seurat 5.5.1, all deps from conda-forge.
# SingleCellExperiment is deliberately ABSENT: nebula's NAMESPACE imports neither it nor
# Seurat (they back the scToNeb() converter only), and run_nebula_stage.R:155 loads just
# Matrix + nebula. Verified by running an actual nebula() NB mixed-model fit, not by import.
: "${SYS_RSCRIPT:=$CONDA_ENVS/comp-nebula/bin/Rscript}"
: "${SYS_RLIB:=$CONDA_ENVS/comp-nebula/lib/R/library}"

export REPO LUAD_SRC_DIR LUAD_PREPPED_DIR LUAD_RESULTS_DIR CONDA_ENVS LUAD_DB LUAD_TIER \
       LUAD_SEED SCRIPTS LUAD_SECTIONS PY_PREP PY_STLEARN PY_SPATIALDM PY_LIANA \
       RS_CELLCHAT RS_CYTOSIGNAL SYS_RSCRIPT SYS_RLIB

# ================================================================================
# INVARIANT: NO REQUESTED / NAMED LR PAIR ON LUAD-AIS.
#
# GRN_SORT1 and ANXA1_FPR1 are the two arms of ALARMIST motif 1 in the LGG/GBM TMA.
# They are a hypothesis about THAT dataset. They must never be requested, tracked,
# plotted into a requested_* directory, or "checked whether they are there" on a
# LUAD/AIS run. The correct value for every requested-LRI knob here is EMPTY.
#
# "The genes are on the LUAD panel, so it is testable" is NOT a reason to include
# them -- panel membership means the assay COULD measure it, not that the hypothesis
# belongs to this dataset. That argument was made in two comments in this tree and
# is wrong; both were removed 2026-08-24.
#
# Consequently, on a LUAD run these are CORRECT OUTPUT, not bugs to be fixed:
#     quant/requested_lr_status.csv    3 bytes (header only)
#     quant/requested_lr_in_DEA.csv    3 bytes
#     data/lr_of_interest_loadings.csv 1 byte
#     plots/requested*/                empty or absent
# Do not repopulate them.
#
# Each runner already neutralises its own knob explicitly -- 01_stlearn.sh and
# 02_spatialdm.sh pass `--requested-lrs ""`, 05_liana.sh passes a bare
# `--lr-of-interest` (argparse nargs="*" -> []), and run_cellchat_luad.R sets
# REQUESTED_LR <- character(0) via make_luad_variants.sh. The line below is the
# belt-and-braces layer for env-driven knobs, added after cytosignal/
# replot_from_quant.R was found defaulting CS_EXTRA to CCI-01109 (= GRN -> SORT1)
# with nothing anywhere overriding it, which put Rank_N_GRN-SORT1.png into all four
# LUAD sections. Every stage sources this file, so setting it here covers them all.
#
# If you ever add a knob of this kind, default it to EMPTY and add it here.
# ================================================================================
CS_EXTRA=""                # cytosignal/replot_from_quant.R -- extra CCI ids to force-plot
export CS_EXTRA

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

    # An executable Rscript is NOT enough for comp-cellchat. On 2026-08-18 the env existed with
    # 224 R packages and this check said OK, while CellChat itself was absent -- svglite had
    # failed to build inside a try() and every stage-04 run died at `library(CellChat)`. Probe
    # the package, not the interpreter. Skipped silently if the env is not built yet.
    if [ -x "$RS_CELLCHAT" ]; then
        echo "== comp-cellchat R packages =="
        if "$RS_CELLCHAT" -e 'q(status = !requireNamespace("CellChat", quietly = TRUE))' >/dev/null 2>&1; then
            printf "  OK    %-14s %s\n" CellChat \
                "$("$RS_CELLCHAT" -e 'cat(as.character(packageVersion("CellChat")))' 2>/dev/null)"
        else
            printf "  MISS  %-14s %s\n" CellChat "installed R library has no CellChat -- the env is HALF-BUILT."
            echo "        Re-run: bash scripts/comparators/cellchat/build_env_iris.sh  (it is idempotent)"
            bad=1
        fi
    fi

    if [ "$bad" -eq 0 ]; then echo "ALL OK"; else echo "SOMETHING IS MISSING -- fix luad_config.sh before running anything"; fi
    return $bad
}

# Executed directly (not sourced) -> run the self-check.
if [ "${BASH_SOURCE[0]}" = "${0}" ]; then luad_check; fi
