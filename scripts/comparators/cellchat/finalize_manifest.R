#!/usr/bin/env Rscript
# finalize_manifest.R -- repair two reporting fields in a run_manifest.json written by the
# first version of run_cellchat.R. NEITHER affects the inference; both are manifest metadata.
#
#  1. observed_min_cell_distance_um was NA for `high` because CellChat's computeCellDistance
#     builds a dense N x N matrix (~51 GB at 79,998 cells). Recomputed here with
#     BiocNeighbors::findKNN(k = 1), which returns the identical quantity in O(N log N).
#  2. peak_gc_mb read gc() column 6, but this R inserts a "limit (Mb)" column, so column 6 is
#     the raw object COUNT, not Mb. The true peak is not recoverable after the process exited,
#     so it is set to null; rss_mb_at_exit was measured correctly and is kept.
#
# run_cellchat.R has been fixed for both; this only backfills runs made before the fix.
#
#   Rscript scripts/comparators/cellchat/finalize_manifest.R --out-dir <tier dir>
suppressWarnings(suppressMessages({
  library(CellChat); library(jsonlite); library(BiocNeighbors)
}))

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(f, d = NULL) { i <- which(args == f); if (!length(i)) d else args[i[1] + 1] }
out_dir <- getarg("--out-dir")
stopifnot(!is.null(out_dir))

mf <- file.path(out_dir, "run_manifest.json")
m <- jsonlite::fromJSON(mf, simplifyVector = FALSE)

for (i in seq_along(m$conditions)) {
  cond <- m$conditions[[i]]$condition
  rds <- file.path(out_dir, "objects", paste0(cond, ".rds"))
  if (!file.exists(rds)) { cat("no object for", cond, "-- skipped\n"); next }
  obj <- readRDS(rds)
  ratio <- obj@images$spatial.factors$ratio[1]
  knn <- BiocNeighbors::findKNN(as.matrix(obj@images$coordinates), k = 1,
                                BNPARAM = BiocNeighbors::KmknnParam(), get.index = FALSE)
  d <- min(knn$distance[, 1]) * ratio
  m$conditions[[i]]$observed_min_cell_distance_um <- round(d, 4)
  m$conditions[[i]]$median_nn_distance_um <- round(stats::median(knn$distance[, 1]) * ratio, 4)
  cat(sprintf("%s: min NN distance %.3f um, median %.3f um (contact.range = %s)\n",
              cond, d, stats::median(knn$distance[, 1]) * ratio, m$parameters$contact.range))
  rm(obj); gc(verbose = FALSE)
}

m$peak_gc_mb <- NULL
m$manifest_patched <- list(
  by = "finalize_manifest.R",
  observed_min_cell_distance_um = "recomputed with BiocNeighbors::findKNN(k=1); CellChat's computeCellDistance is dense O(N^2) and OOMs at 79,998 cells",
  peak_gc_mb = "dropped -- original value indexed gc() column 6, which is an object count on this R (a 'limit (Mb)' column shifts the Mb column); true peak not recoverable post hoc. rss_mb_at_exit is unaffected.")

writeLines(jsonlite::toJSON(m, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null"), mf)
cat("patched", mf, "\n")
