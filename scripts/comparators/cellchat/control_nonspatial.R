#!/usr/bin/env Rscript
# control_nonspatial.R -- decisive control: does the spatial information change ANYTHING?
#
# Motivation: in the GBM runs the 250 um `interaction.range` filter excluded 0 of 81 cell-type
# pairs (d_spatial has no NaN in any run), and the governing multi-sample vignette sets
# `distance.use = FALSE`. Reading modeling.R:161-165, that combination makes P.spatial an
# all-ones matrix -- so for L-R pairs that take the "diffusible" branch the spatial constraint
# multiplies nothing. This script tests that claim empirically instead of by code-reading:
# it reruns the identical pipeline with `datatype = "RNA"` (no coordinates at all) and diffs
# the probability array against the spatial run.
#
# If they are bit-identical, the `default` tier used no spatial information whatsoever.
#
#   source scripts/comparators/cellchat/activate_env.sh
#   Rscript scripts/comparators/cellchat/control_nonspatial.R \
#     --input-dir results/comparators/cellchat/GBM/input \
#     --spatial-dir results/comparators/cellchat/GBM/default \
#     --out-dir results/comparators/cellchat/GBM/control_nonspatial \
#     --tier default --conditions low
suppressWarnings(suppressMessages({
  library(CellChat); library(Matrix); library(jsonlite); library(future)
}))
options(stringsAsFactors = FALSE)
SCRIPT_DIR <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(SCRIPT_DIR, "cellchat_io.R"))

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(f, d = NULL) { i <- which(args == f); if (!length(i)) d else args[i[1] + 1] }
input_dir   <- getarg("--input-dir")
spatial_dir <- getarg("--spatial-dir")
out_dir     <- getarg("--out-dir")
tier        <- getarg("--tier", "default")
conditions  <- strsplit(getarg("--conditions", "low"), ",")[[1]]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

CellChatDB <- CellChatDB.human
CellChatDB.use <- if (tier == "default")
  subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") else subsetDB(CellChatDB)

report <- list()
for (cond in conditions) {
  d <- file.path(input_dir, cond)
  m <- Matrix::readMM(file.path(d, "data.mtx"))
  genes <- readLines(file.path(d, "genes.tsv")); cells <- readLines(file.path(d, "barcodes.tsv"))
  meta <- utils::read.csv(file.path(d, "meta.csv"), colClasses = "character")
  data.input <- as(as(m, "CsparseMatrix"), "dgCMatrix")
  dimnames(data.input) <- list(genes, cells)
  meta.cc <- data.frame(labels = droplevels(factor(meta$labels)),
                        samples = droplevels(factor(meta$samples)), row.names = meta$cell_id)

  # THE ONLY DIFFERENCE: datatype = "RNA", no coordinates, no spatial.factors
  cc <- createCellChat(object = data.input, meta = meta.cc, group.by = "labels", datatype = "RNA")
  cc@DB <- CellChatDB.use
  cc <- subsetData(cc)
  future::plan("multisession", workers = 4)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, type = "truncatedMean", trim = 0.1, nboot = 100, seed.use = 1L)
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  save_cellchat_quant(cc, cond, out_dir)

  # diff against the spatial run
  A <- utils::read.csv(gzfile(file.path(spatial_dir, "quant", sprintf("%s_net_full.csv.gz", cond))))
  B <- utils::read.csv(gzfile(file.path(out_dir,     "quant", sprintf("%s_net_full.csv.gz", cond))))
  key <- c("source", "target", "interaction_name")
  mm <- merge(A[, c(key, "prob", "pval")], B[, c(key, "prob", "pval")], by = key,
              suffixes = c("_spatial", "_rna"))
  dp <- max(abs(mm$prob_spatial - mm$prob_rna)); dv <- max(abs(mm$pval_spatial - mm$pval_rna))
  report[[cond]] <- list(condition = cond,
    n_rows_spatial = nrow(A), n_rows_rna = nrow(B), n_shared = nrow(mm),
    only_spatial = nrow(A) - nrow(mm), only_rna = nrow(B) - nrow(mm),
    max_abs_delta_prob = dp, max_abs_delta_pval = dv,
    identical = (nrow(A) == nrow(B)) && (nrow(mm) == nrow(A)) && dp == 0 && dv == 0)
  cat(sprintf("%s: spatial %d rows vs RNA %d rows, shared %d | max|dprob|=%.3e max|dpval|=%.3e | IDENTICAL=%s\n",
              cond, nrow(A), nrow(B), nrow(mm), dp, dv, report[[cond]]$identical))
}
writeLines(jsonlite::toJSON(list(
  question = "Does the spatial information change the inferred communication at all?",
  method = "identical pipeline with datatype='RNA' (no coordinates) vs the spatial run",
  tier = tier, spatial_dir = spatial_dir, results = unname(report)),
  auto_unbox = TRUE, pretty = TRUE), file.path(out_dir, "spatial_vs_nonspatial.json"))
cat("wrote", file.path(out_dir, "spatial_vs_nonspatial.json"), "\n")
