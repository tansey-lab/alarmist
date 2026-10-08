# Source this to use the comp-cellchat conda R env WITHOUT relying on `conda activate`.
#
# Setting PATH directly is the portable route: on the original MacBook `conda activate`
# and `conda run` were broken (libmamba/libarchive) and fell through to the *system* R,
# and on any host it is the one method that cannot silently pick up the wrong R. Also
# sources the conda-forge compiler activation scripts so CC/CXX/FC are set when building
# packages from source.
#
#   source scripts/comparators/cellchat/activate_env.sh
#   Rscript scripts/comparators/cellchat/run_cellchat.R ...
#
# The env location is NOT hardcoded: it is derived from CONDA_ENVS in
# _common/luad_config.sh, which is the single file to edit when moving machines.

_HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../_common/luad_config.sh
source "$_HERE/../_common/luad_config.sh"

ENV="$CONDA_ENVS/comp-cellchat"
[ -d "$ENV" ] || echo "WARNING: $ENV does not exist -- check CONDA_ENVS in _common/luad_config.sh"

export CONDA_PREFIX="$ENV"
export PATH="$ENV/bin:$PATH"

# Pin R's user library to the env's own. Without this, ~/R/x86_64-*-library/<ver> is shared
# across every conda R env on the host and packages built against a different R leak in --
# the classic "installed but won't load" failure. Verify with: Rscript -e '.libPaths()'
export R_LIBS_USER="$ENV/lib/R/library"

for f in "$ENV"/etc/conda/activate.d/*.sh; do source "$f" 2>/dev/null; done
unset _HERE
