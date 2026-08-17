#!/usr/bin/env bash
# ======================================================================================
# Rebuild the six conda environments this LUAD/AIS run needs, on a fresh host (iris).
#
#   bash scripts/comparators/_common/install_envs_iris.sh [env ...]
#   bash scripts/comparators/_common/install_envs_iris.sh --check
#
# ONE ENV PER METHOD IS A HARD RULE, not a preference: their numpy / scanpy / Seurat pins
# conflict and co-installing two methods breaks both. See scripts/comparators/METHODS.md.
#
# Envs built here (COMMOT and NICHES are deliberately absent -- they are out of scope for
# this run because their dense O(N^2) edgelists need 0.5-9 TB on these sections):
#
#   comp-stlearn     python   from stlearn/env.lock.yml
#   comp-spatialdm   python   from spatialdm/env.lock.yml
#   comp-liana       python   from liana/env.lock.yml
#   comp-cellchat    R        from cellchat/env.lock.yml + install_env.R  (+ a CellChat clone)
#   comp-cytosignal  R        NO lock file exists -- built from the recipe below
#   comp-nebula      R        NEW: a minimal R that has `nebula`, which does not build inside
#                             comp-cytosignal. run_nebula_stage.R re-invokes itself here via
#                             SYS_RSCRIPT / SYS_RLIB.
#   bptf             python   the prep/driver env (anndata, scanpy, h5py, scipy, pandas)
#
# THREE THINGS THAT DO NOT COME FROM A LOCK FILE AND WILL BITE:
#
#   1. CellChat is installed as RemoteType: local from a working-copy clone. You must copy
#      the CellChat source tree to iris and point CELLCHAT_SRC at it. Version 2.2.0.9001.
#   2. CytoSignal is pinned to GitHub SHA b2c1a3cfda258fa57f6c8c3a92635dfe4b0505a0 on the
#      source machine. The old bundle installer called install_github() WITHOUT a ref, which
#      would silently give you upstream HEAD. This script pins the ref.
#   3. NICHES needs a vendored seurat-wrappers -- not relevant here, NICHES is out of scope.
#
# conda solver: libmamba is broken on the source machine, hence CONDA_SOLVER=classic. It is
# harmless elsewhere; drop it if iris's mamba is healthy and you want the speed.
# ======================================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="$(cd "$HERE/.." && pwd)"

: "${CONDA:=conda}"
: "${CONDA_SOLVER:=classic}"
: "${CHANNELS:=-c conda-forge --override-channels}"
: "${CELLCHAT_SRC:=}"          # path to the CellChat source clone, copied from the source host
: "${CYTOSIGNAL_SHA:=b2c1a3cfda258fa57f6c8c3a92635dfe4b0505a0}"
export CONDA_SOLVER

have_env() { $CONDA env list | awk '{print $1}' | grep -qx "$1"; }

from_lock() {  # from_lock <env> <lockfile>
    local env="$1" lock="$2"
    if have_env "$env"; then echo ">> $env already exists, skipping"; return 0; fi
    [ -f "$lock" ] || { echo "!! $lock not found -- cannot build $env"; return 1; }
    echo ">> creating $env from $lock"
    $CONDA env create -n "$env" -f "$lock"
}

build_cytosignal() {
    if have_env comp-cytosignal; then echo ">> comp-cytosignal already exists, skipping"; return 0; fi
    echo ">> creating comp-cytosignal (no lock file exists; recipe is in this script)"
    # shellcheck disable=SC2086
    $CONDA create -y -n comp-cytosignal $CHANNELS \
        r-base=4.3 r-remotes r-biocmanager \
        r-rcpp r-rcpparmadillo r-rcppprogress r-matrix r-rann r-irlba \
        r-dbscan r-glmnet r-ggrepel r-ggplot2 r-dplyr r-data.table r-reshape2 \
        r-gridextra r-scales r-plyr r-rlang r-scattermore r-cowplot \
        compilers make pkg-config
    $CONDA run -n comp-cytosignal Rscript -e "
      options(repos = c(CRAN='https://cloud.r-project.org'), Ncpus = 4)
      # PIN THE REF. Without it you get upstream HEAD, which is not what the GBM results used.
      remotes::install_github('welch-lab/cytosignal', ref='${CYTOSIGNAL_SHA}', upgrade='never')
      remotes::install_github('xzhoulab/SPARK', upgrade='never')   # rankIntrSpatialVar needs it
      cat('cytosignal:', as.character(packageVersion('cytosignal')),
          ' SPARK:', requireNamespace('SPARK', quietly=TRUE), '\n')"
}

build_nebula() {
    if have_env comp-nebula; then echo ">> comp-nebula already exists, skipping"; return 0; fi
    echo ">> creating comp-nebula (the second R, for run_nebula_stage.R stage 2)"
    # shellcheck disable=SC2086
    $CONDA create -y -n comp-nebula $CHANNELS \
        r-base=4.3 r-matrix r-nebula compilers make pkg-config || {
            echo "   r-nebula not on conda-forge here; falling back to install.packages"
            $CONDA create -y -n comp-nebula $CHANNELS r-base=4.3 r-matrix compilers make pkg-config
            $CONDA run -n comp-nebula Rscript -e "
              options(repos=c(CRAN='https://cloud.r-project.org'), Ncpus=4)
              install.packages('nebula')"
        }
    $CONDA run -n comp-nebula Rscript -e "
      library(nebula); cat('nebula:', as.character(packageVersion('nebula')), '\n')"
}

build_cellchat() {
    from_lock comp-cellchat "$SCRIPTS/cellchat/env.lock.yml" || return 1
    if [ -z "$CELLCHAT_SRC" ]; then
        cat <<'EOF'
!! CELLCHAT_SRC is not set.
   CellChat is installed as RemoteType: local -- the lock file does NOT contain it.
   Copy the CellChat source tree from the source host and re-run:
       CELLCHAT_SRC=/path/to/CellChat bash install_envs_iris.sh comp-cellchat
EOF
        return 1
    fi
    echo ">> installing CellChat from $CELLCHAT_SRC into comp-cellchat"
    $CONDA run -n comp-cellchat Rscript -e "
      options(repos=c(CRAN='https://cloud.r-project.org'), Ncpus=4)
      if (!requireNamespace('remotes', quietly=TRUE)) install.packages('remotes')
      remotes::install_local('${CELLCHAT_SRC}', upgrade='never')
      library(CellChat); cat('CellChat:', as.character(packageVersion('CellChat')), '\n')
      cat('presto:', requireNamespace('presto', quietly=TRUE), '\n')"
    # presto is not optional: computeCommunProb's DEA falls back to stats::wilcox.test
    # SILENTLY without it, and CellChat's thresh.fc was tuned for presto.
    $CONDA run -n comp-cellchat Rscript -e "
      if (!requireNamespace('presto', quietly=TRUE)) {
        remotes::install_github('immunogenomics/presto', upgrade='never')
      }
      cat('presto:', as.character(packageVersion('presto')), '\n')"
}

check() {
    echo "== environment check =="
    local bad=0
    for e in bptf comp-stlearn comp-spatialdm comp-liana comp-cellchat comp-cytosignal comp-nebula; do
        if have_env "$e"; then printf "  OK    %s\n" "$e"; else printf "  MISS  %s\n" "$e"; bad=1; fi
    done
    echo "== key packages =="
    $CONDA run -n bptf           python -c "import anndata,scanpy,h5py;print('  bptf       anndata',anndata.__version__)" 2>/dev/null || { echo "  bptf       FAIL"; bad=1; }
    $CONDA run -n comp-stlearn   python -c "import stlearn;print('  stlearn   ',stlearn.__version__)" 2>/dev/null || { echo "  stlearn    FAIL"; bad=1; }
    $CONDA run -n comp-spatialdm python -c "import spatialdm;print('  spatialdm ',spatialdm.__version__)" 2>/dev/null || { echo "  spatialdm  FAIL"; bad=1; }
    $CONDA run -n comp-liana     python -c "import liana,mofaflex;print('  liana     ',liana.__version__)" 2>/dev/null || { echo "  liana      FAIL"; bad=1; }
    $CONDA run -n comp-cellchat  Rscript -e "cat('  CellChat  ',as.character(packageVersion('CellChat')),'\n')" 2>/dev/null || { echo "  CellChat   FAIL"; bad=1; }
    $CONDA run -n comp-cytosignal Rscript -e "cat('  cytosignal',as.character(packageVersion('cytosignal')),'\n')" 2>/dev/null || { echo "  cytosignal FAIL"; bad=1; }
    $CONDA run -n comp-nebula    Rscript -e "cat('  nebula    ',as.character(packageVersion('nebula')),'\n')" 2>/dev/null || { echo "  nebula     FAIL"; bad=1; }
    [ "$bad" -eq 0 ] && echo "ALL OK" || echo "SOMETHING IS MISSING"
    return $bad
}

TARGETS=("$@")
if [ "${1:-}" = "--check" ]; then check; exit $?; fi
[ ${#TARGETS[@]} -eq 0 ] && TARGETS=(bptf comp-stlearn comp-spatialdm comp-liana comp-cellchat comp-cytosignal comp-nebula)

for t in "${TARGETS[@]}"; do
    case "$t" in
        bptf)            have_env bptf || echo "!! bptf must be created separately (anndata, scanpy, h5py, scipy, pandas; numpy<=1.24 for numba)" ;;
        comp-stlearn)    from_lock comp-stlearn   "$SCRIPTS/stlearn/env.lock.yml" ;;
        comp-spatialdm)  from_lock comp-spatialdm "$SCRIPTS/spatialdm/env.lock.yml" ;;
        comp-liana)      from_lock comp-liana     "$SCRIPTS/liana/env.lock.yml" ;;
        comp-cellchat)   build_cellchat ;;
        comp-cytosignal) build_cytosignal ;;
        comp-nebula)     build_nebula ;;
        *) echo "unknown target: $t" ;;
    esac
done

echo
check || true
