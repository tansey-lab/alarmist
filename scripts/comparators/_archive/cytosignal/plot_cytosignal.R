#!/usr/bin/env Rscript
# Regenerate CytoSignal spatial plots from the pre-signif checkpoint (has score/imputation
# slots that purgeBeforeSave strips from cs_result.rds). Re-runs the cheap signif+rank tail.
# Usage: Rscript plot_cytosignal.R <out_dir>
suppressWarnings(suppressMessages({ library(Matrix); library(cytosignal) }))
args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[1] else "results/comparators/cytosignal/LUAD/default/run_crop_2p5mm"
log <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n")

cs <- readRDS(file.path(out_dir, "cs_checkpoint.rds"))
cs <- inferSignif(cs, p.thresh = 0.05, reads.thresh = 100, sig.thresh = 100)
cs <- rankIntrSpatialVar(cs, numCores = 4)
pdir <- file.path(out_dir, "plots"); dir.create(pdir, showWarnings = FALSE, recursive = TRUE)

# cluster map
tryCatch({ png(file.path(pdir, "cluster_map.png"), 1500, 1300, res = 150); print(plotCluster(cs)); dev.off() },
         error = function(e) log("plotCluster err:", conditionMessage(e)))

# top spatially-variable interactions per mode -> plotSignif writes into plot_dir
for (slot in c("diffusion-Raw_smooth", "contact-Raw_smooth")) {
  tag <- gsub("[^A-Za-z0-9]+", "_", slot)
  n <- tryCatch(length(showIntr(cs, slot.use = slot, signif.use = "result.spx", return.name = TRUE)),
                error = function(e) 0)
  ntop <- min(6, n); log(slot, "-> plotting top", ntop, "of", n)
  if (ntop > 0) tryCatch(
    plotSignif(cs, intr = seq_len(ntop), slot.use = slot, signif.use = "result.spx",
               plot_dir = file.path(pdir, paste0("signif_", tag, "/")),
               pt.size = 0.25, plot.fmt = "png", raster = TRUE),
    error = function(e) log("plotSignif err", slot, ":", conditionMessage(e)))
}
log("DONE. plots in", pdir)
