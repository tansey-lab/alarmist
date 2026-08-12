# Source this to use the comp-cellchat conda R env WITHOUT relying on `conda activate`.
# This machine's conda is degraded (libmamba/libarchive broken), so `conda activate` and
# `conda run` fall through to the *system* R 4.4.2 instead of the env's R 4.3.3. Setting
# PATH directly is the only reliable route. Also sources the conda-forge compiler
# activation scripts so CC/CXX/FC/SDKROOT are set when building packages from source.
#
#   source scripts/comparators/cellchat/activate_env.sh
#   Rscript scripts/comparators/cellchat/run_cellchat.R ...
ENV=/Users/jiayifan/anaconda3/envs/comp-cellchat
export CONDA_PREFIX="$ENV"
export PATH="$ENV/bin:$PATH"
for f in "$ENV"/etc/conda/activate.d/*.sh; do source "$f" 2>/dev/null; done
