#!/usr/bin/env bash
# ======================================================================================
# install_plot_deps.sh -- add CytoSignal's OPTIONAL plotting dependencies to comp-cytosignal.
#
#   bash scripts/comparators/cytosignal/install_plot_deps.sh
#
# WHY: `circlize` is in cytosignal's DESCRIPTION under Suggests, not Imports, and the env
# recipe (_common/install_envs_iris.sh:62-67 `build_cytosignal`) does not list it. So
# `plotCircosNIntr` -- the "which cell type talks to which" chord diagram, and the only
# cross-cell-type summary in CytoSignal's plot surface -- has never been able to run in
# this environment. It fails with the package's own
#     stop("Package 'circlize' is required for this function.")
# which plot_signif_rerun.R and replot_from_quant.R both swallow in a tryCatch, so it has
# been failing SILENTLY. Confirmed by absence: no circos output exists anywhere under
# results/, including from the completed GBM run.
#
# The other two Suggests that the plotting path actually reaches -- cowplot (plotSignif)
# and scattermore (plotIntrValue, raster = TRUE) -- ARE present, so nothing else is gated.
# plot3D / png / pheatmap are NOT needed: they appear only as string literals in
# `plot.fmt = c("png","pdf","svg")`-style arguments, never as namespace calls.
#
# Additive and idempotent: installs into the env's OWN library, upgrades nothing, and
# re-running is a no-op. Matches the convention already used for this env -- plain
# install.packages() against cloud.r-project.org (install_envs_iris.sh:69).
#
# Deliberately NOT folded into install_envs_iris.sh: that script records what actually
# built the existing envs, and editing its recipe would make a fresh build diverge from
# the env every current result came from. Declared in cytosignal/DEVIATIONS.md instead.
# ======================================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/activate_env.sh"
export TMPDIR="${TMPDIR:-/data1/tanseyw/projects/fanj2/tmp}"
mkdir -p "$TMPDIR"

echo ">> R           : $(command -v Rscript)"
echo ">> library     : $R_LIBS_USER"
echo ">> TMPDIR      : $TMPDIR"

Rscript - <<'RCODE'
options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 4)
lib <- Sys.getenv("R_LIBS_USER")
stopifnot(nzchar(lib), dir.exists(lib))
want <- c("circlize")
for (p in want) {
  if (requireNamespace(p, quietly = TRUE)) {
    cat(sprintf("   %-12s already present (%s)\n", p, as.character(packageVersion(p))))
  } else {
    cat(sprintf("   %-12s installing...\n", p))
    install.packages(p, lib = lib, dependencies = c("Depends", "Imports", "LinkingTo"))
  }
}
cat("\n== verify ==\n")
ok <- TRUE
for (p in want) {
  have <- requireNamespace(p, quietly = TRUE)
  cat(sprintf("   %-12s %s\n", p, if (have) as.character(packageVersion(p)) else "STILL MISSING"))
  ok <- ok && have
}
# Prove the actual entry point resolves, not just that the namespace loads.
if (ok) {
  suppressMessages(library(cytosignal))
  cat("   plotCircosNIntr reachable:",
      is.function(tryCatch(cytosignal::plotCircosNIntr, error = function(e) NULL)), "\n")
}
quit(status = if (ok) 0 else 1)
RCODE
echo ">> done"
