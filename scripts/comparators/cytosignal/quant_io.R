# quant_io.R — persist CytoSignal QUANTITATIVE outputs (not just interaction names).
# Per *_smooth score slot, into out_dir:
#   score_<slot>.mtx (+ .cells.tsv, .intr.tsv)   : cells x interactions LR-score (small runs, portable)
#   score_<slot>.rds                              : same dgCMatrix w/ dimnames (large runs; @score is ~dense)
#   reslist_<slot>.rds                            : full @res.list (result / result.hq / result.spx) =
#                                                   per-interaction character vector of significant cell barcodes
#   signif_summary_<slot>.csv                     : interaction_id, name, n_result, n_hq, n_spx
suppressWarnings(suppressMessages(library(Matrix)))
`%||%` <- function(a, b) if (is.null(a)) b else a
.ids <- function(x) {                     # interaction ids from a res.list element
  if (is.null(x)) return(character(0))
  if (is.list(x) && !is.null(names(x))) return(names(x))
  if (!is.null(dim(x)) && !is.null(rownames(x))) return(rownames(x))
  names(x) %||% as.character(x)
}

save_lrscore_quant <- function(cs, out_dir, intr_names = NULL, slots = NULL,
                               max_mtx_cells = 100000, log = message) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  if (is.null(slots)) slots <- grep("_smooth$", names(cs@lrscore), value = TRUE)
  for (slot in slots) {
    lr  <- cs@lrscore[[slot]]
    tag <- gsub("[^A-Za-z0-9]+", "_", slot)
    saveRDS(lr@res.list, file.path(out_dir, sprintf("reslist_%s.rds", tag)))

    sig <- unique(unlist(lapply(lr@res.list, .ids)))                # interaction ids
    S <- lr@score                                                  # orient to cells x interactions
    if (!is.null(rownames(S)) && length(intersect(rownames(S), sig)) &&
        !(!is.null(colnames(S)) && length(intersect(colnames(S), sig)))) S <- t(S)
    S <- as(as(S, "CsparseMatrix"), "dgCMatrix")
    if (nrow(S) <= max_mtx_cells) {
      Matrix::writeMM(S, file.path(out_dir, sprintf("score_%s.mtx", tag)))
      writeLines(rownames(S) %||% as.character(seq_len(nrow(S))), file.path(out_dir, sprintf("score_%s.cells.tsv", tag)))
      writeLines(colnames(S) %||% as.character(seq_len(ncol(S))), file.path(out_dir, sprintf("score_%s.intr.tsv", tag)))
      fmt <- "mtx"
    } else {
      saveRDS(S, file.path(out_dir, sprintf("score_%s.rds", tag)), compress = TRUE)   # dimnames embedded
      fmt <- "rds"
    }
    log(sprintf("[quant] %s: score %d cells x %d intr (%s)", slot, nrow(S), ncol(S), fmt))

    all_ids <- unique(unlist(lapply(lr@res.list, .ids)))
    cnt <- function(nm) { x <- lr@res.list[[nm]]; if (is.null(x) || !length(x)) integer(0) else vapply(x, length, integer(1)) }
    nr <- cnt("result"); nh <- cnt("result.hq"); ns <- cnt("result.spx")
    g  <- function(v, ids) { out <- v[ids]; out[is.na(out)] <- 0L; as.integer(out) }
    summ <- data.frame(interaction_id = all_ids,
                       name     = if (!is.null(intr_names)) unname(intr_names[all_ids]) else all_ids,
                       n_result = g(nr, all_ids), n_hq = g(nh, all_ids), n_spx = g(ns, all_ids),
                       stringsAsFactors = FALSE)
    summ <- summ[order(-summ$n_hq), ]
    write.csv(summ, file.path(out_dir, sprintf("signif_summary_%s.csv", tag)), row.names = FALSE)
  }
}
