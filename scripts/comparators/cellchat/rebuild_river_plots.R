#!/usr/bin/env Rscript
# rebuild_river_plots.R -- regenerate only netAnalysis_river (VB:482,499).
#
# The first plotting pass lost these because `ggalluvial` was not ATTACHED: netAnalysis_river
# builds a ggalluvial plot with stat = "stratum", and ggplot2 resolves stats by name from the
# search path, so `ggalluvial::` is not enough -- it fails with "Can't find stat called
# 'stratum'". VB:469 loads it explicitly; plot_cellchat.R now does too. This script exists so
# the fix does not require re-rendering thousands of unaffected plots.
#
# k is re-derived from the persisted selectK measures with the same rule plot_cellchat.R uses,
# so the patterns are identical to the ones already written to quant/.
#
#   source scripts/comparators/cellchat/activate_env.sh
#   Rscript scripts/comparators/cellchat/rebuild_river_plots.R \
#     --out-dir results/comparators/cellchat/GBM/default --conditions low,high
suppressWarnings(suppressMessages({
  library(CellChat); library(NMF); library(ggalluvial); library(ComplexHeatmap)
}))
SCRIPT_DIR <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(SCRIPT_DIR, "cellchat_io.R"))

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(f, d = NULL) { i <- which(args == f); if (!length(i)) d else args[i[1] + 1] }
out_dir    <- getarg("--out-dir")
conditions <- strsplit(getarg("--conditions", "low,high"), ",")[[1]]
stopifnot(!is.null(out_dir))

pick_k <- function(meas) {
  co <- meas$cophenetic; ks <- meas$rank
  i <- which(diff(co) < -0.01)[1]
  if (is.na(i)) ks[which.max(co)] else ks[i]
}

for (cond in conditions) {
  obj <- readRDS(file.path(out_dir, "objects", paste0(cond, ".rds")))
  for (pat in c("outgoing", "incoming")) {
    mf <- file.path(out_dir, "quant", sprintf("%s_selectK_%s_measures.csv", cond, pat))
    if (!file.exists(mf)) { message("no measures for ", cond, " ", pat); next }
    k <- pick_k(utils::read.csv(mf))
    o2 <- identifyCommunicationPatterns(obj, pattern = pat, k = k, heatmap.show = FALSE)
    ok <- save_all_formats(function() netAnalysis_river(o2, pattern = pat),
                           file.path(out_dir, "plots", cond, "systems", paste0("river_", pat)),
                           9, 6)
    cat(sprintf("%s %s: k=%d river=%s\n", cond, pat, k, if (isTRUE(ok)) "OK" else "FAILED"))
  }
  rm(obj); gc(verbose = FALSE)
}
