# niches_io.R -- persist NICHES QUANTITATIVE outputs, not just plots and mechanism names.
#
# Per NICHES organization (NeighborhoodToCell / CellToCellSpatial), into out_dir:
#   objects/<org>.rds                  full Seurat object -- everything downstream re-reads this,
#                                      so a rerun never has to recompute the O(N^2) edgelist
#   quant/<org>_scores.mtx             sparse mechanisms x columns LR-score matrix
#   quant/<org>_features.tsv           mechanism names (rows), em-dash joined
#   quant/<org>_columns.tsv            column names (cells for N2C, cell-cell edges for C2CS)
#   quant/<org>_metadata.csv           the object's full meta.data
#   quant/<org>_mechanism_summary.csv  per mechanism: mean/max score, n and frac nonzero columns
suppressWarnings(suppressMessages({ library(Matrix); library(Seurat) }))

`%||%` <- function(a, b) if (is.null(a)) b else a

# NICHES writes its own assay data with both `counts` and `data` holding the raw LR products
# (RunNeighborhoodToCell.R:77-81), so either layer is the same numbers. We read `data`.
niches_matrix <- function(obj, org) {
  m <- if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    SeuratObject::LayerData(obj, assay = org, layer = "data")
  } else {
    Seurat::GetAssayData(obj, assay = org, slot = "data")
  }
  as(as(as(m, "CsparseMatrix"), "generalMatrix"), "dgCMatrix")
}

save_niches_quant <- function(obj, org, out_dir, log = message) {
  qdir <- file.path(out_dir, "quant")
  odir <- file.path(out_dir, "objects")
  dir.create(qdir, showWarnings = FALSE, recursive = TRUE)
  dir.create(odir, showWarnings = FALSE, recursive = TRUE)

  saveRDS(obj, file.path(odir, sprintf("%s.rds", org)), compress = TRUE)

  M <- niches_matrix(obj, org)
  Matrix::writeMM(M, file.path(qdir, sprintf("%s_scores.mtx", org)))
  writeLines(rownames(M), file.path(qdir, sprintf("%s_features.tsv", org)))
  writeLines(colnames(M), file.path(qdir, sprintf("%s_columns.tsv", org)))

  md <- obj@meta.data
  md$column_id <- rownames(md)
  utils::write.csv(md, file.path(qdir, sprintf("%s_metadata.csv", org)), row.names = FALSE)

  nz <- Matrix::rowSums(M > 0)
  summ <- data.frame(
    mechanism  = rownames(M),
    mean_score = as.numeric(Matrix::rowMeans(M)),
    max_score  = apply(M, 1, max),
    n_nonzero  = as.integer(nz),
    frac_nonzero = as.numeric(nz) / ncol(M),
    stringsAsFactors = FALSE
  )
  summ <- summ[order(-summ$frac_nonzero, -summ$mean_score), ]
  utils::write.csv(summ, file.path(qdir, sprintf("%s_mechanism_summary.csv", org)), row.names = FALSE)

  log(sprintf("[quant] %s: %d mechanisms x %d columns, %.1f%% nonzero entries",
              org, nrow(M), ncol(M), 100 * Matrix::nnzero(M) / prod(dim(M))))
  invisible(list(n_mechanisms = nrow(M), n_columns = ncol(M),
                 nnz = Matrix::nnzero(M), density = Matrix::nnzero(M) / prod(dim(M))))
}
