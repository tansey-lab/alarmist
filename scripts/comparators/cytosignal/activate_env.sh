# Source this to use the comp-cytosignal conda R env WITHOUT relying on `conda activate`
# (this machine's conda is degraded: libmamba/libarchive broken, so activate/run fall
#  through to system R 4.4.2). This sets PATH to the env's R 4.3.3 and sources the
# conda-forge compiler activation scripts so CC/CXX/FC/SDKROOT are set for building pkgs.
ENV=/Users/jiayifan/anaconda3/envs/comp-cytosignal
export CONDA_PREFIX="$ENV"
export PATH="$ENV/bin:$PATH"
for f in "$ENV"/etc/conda/activate.d/*.sh; do source "$f" 2>/dev/null; done
