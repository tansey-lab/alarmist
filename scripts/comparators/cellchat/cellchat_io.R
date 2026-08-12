# cellchat_io.R -- persist CellChat QUANTITATIVE outputs, not just plots and pathway names.
#
# CellChat's unit is a (sender cell type, receiver cell type, L-R pair) triple, so everything
# below is written in that shape. Per condition object, into <out_dir>:
#
#   objects/<cond>.rds              the full CellChat object -- plotting re-reads this, so a
#                                   re-plot never has to redo the 100-permutation test
#   quant/<cond>_net_full.csv.gz    EVERY source x target x LR cell: prob + pval (not only the
#                                   significant ones subsetCommunication returns)
#   quant/<cond>_net_significant.csv   subsetCommunication(object)          -- LR level
#   quant/<cond>_netP_significant.csv  subsetCommunication(slot.name="netP")-- pathway level
#   quant/<cond>_netP_prob.csv.gz   source x target x pathway probabilities
#   quant/<cond>_net_count.csv      cell type x cell type, # significant links (aggregateNet)
#   quant/<cond>_net_weight.csv     cell type x cell type, summed probability
#   quant/<cond>_centrality.csv     per pathway x cell type: outdeg/indeg/flowbet/info
#   quant/<cond>_d_spatial.csv      cell type x cell type mean centroid distance (um)
#   quant/<cond>_adj_contact.csv    cell type x cell type contact adjacency
#   quant/<cond>_group_sizes.csv    cells per cell type, and per cell type x sample
#   quant/<cond>_db_used.csv        the exact interaction table that was used
#   quant/<cond>_options.json       object@options$parameter -- every inference parameter
suppressWarnings(suppressMessages({ library(CellChat); library(jsonlite) }))

`%||%` <- function(a, b) if (is.null(a)) b else a

# 3-D array [source, target, third] -> long data.frame. Keeps zeros out to stay small.
array3_to_long <- function(arr, third.name, value.name, drop.zero = TRUE) {
  if (is.null(arr) || length(dim(arr)) != 3) return(NULL)
  dn <- dimnames(arr)
  idx <- which(if (drop.zero) arr != 0 else array(TRUE, dim(arr)), arr.ind = TRUE)
  if (nrow(idx) == 0) return(NULL)
  df <- data.frame(
    source = dn[[1]][idx[, 1]],
    target = dn[[2]][idx[, 2]],
    x      = dn[[3]][idx[, 3]],
    value  = arr[idx],
    stringsAsFactors = FALSE
  )
  colnames(df)[3:4] <- c(third.name, value.name)
  df
}

save_cellchat_quant <- function(object, cond, out_dir, log = message) {
  qdir <- file.path(out_dir, "quant"); odir <- file.path(out_dir, "objects")
  dir.create(qdir, showWarnings = FALSE, recursive = TRUE)
  dir.create(odir, showWarnings = FALSE, recursive = TRUE)
  p <- function(...) file.path(qdir, sprintf(...))

  saveRDS(object, file.path(odir, sprintf("%s.rds", cond)), compress = TRUE)

  # ---- the full probability / p-value arrays (NOT just the significant subset) ----
  prob <- object@net$prob; pval <- object@net$pval
  dfp <- array3_to_long(prob, "interaction_name", "prob", drop.zero = TRUE)
  if (!is.null(dfp)) {
    pv <- array3_to_long(pval, "interaction_name", "pval", drop.zero = FALSE)
    dfp <- merge(dfp, pv, by = c("source", "target", "interaction_name"), all.x = TRUE)
    # annotate with pathway / ligand / receptor from the DB actually used
    inter <- object@DB$interaction
    key <- match(dfp$interaction_name, inter$interaction_name)
    dfp$pathway_name <- inter$pathway_name[key]
    dfp$ligand       <- inter$ligand[key]
    dfp$receptor     <- inter$receptor[key]
    dfp$annotation   <- inter$annotation[key]
    dfp <- dfp[order(-dfp$prob), ]
    utils::write.csv(dfp, gzfile(p("%s_net_full.csv.gz", cond)), row.names = FALSE)
  }

  # ---- the significant communications, both levels ----
  df.net  <- tryCatch(subsetCommunication(object), error = function(e) NULL)
  if (!is.null(df.net)) utils::write.csv(df.net, p("%s_net_significant.csv", cond), row.names = FALSE)
  df.netP <- tryCatch(subsetCommunication(object, slot.name = "netP"), error = function(e) NULL)
  if (!is.null(df.netP)) utils::write.csv(df.netP, p("%s_netP_significant.csv", cond), row.names = FALSE)

  dfPp <- array3_to_long(object@netP$prob, "pathway_name", "prob", drop.zero = TRUE)
  if (!is.null(dfPp)) utils::write.csv(dfPp, gzfile(p("%s_netP_prob.csv.gz", cond)), row.names = FALSE)

  # ---- aggregated networks ----
  if (!is.null(object@net$count))  utils::write.csv(object@net$count,  p("%s_net_count.csv", cond))
  if (!is.null(object@net$weight)) utils::write.csv(object@net$weight, p("%s_net_weight.csv", cond))

  # ---- network centrality, per pathway ----
  centr <- object@netP$centr
  if (length(centr) > 0) {
    rows <- do.call(rbind, lapply(names(centr), function(pw) {
      m <- centr[[pw]]
      data.frame(pathway_name = pw, cell_type = names(m$outdeg),
                 outdeg = as.numeric(m$outdeg), indeg = as.numeric(m$indeg),
                 flowbet = as.numeric(m$flowbet %||% rep(NA, length(m$outdeg))),
                 info = as.numeric(m$info %||% rep(NA, length(m$outdeg))),
                 stringsAsFactors = FALSE)
    }))
    utils::write.csv(rows, p("%s_centrality.csv", cond), row.names = FALSE)
  }

  # ---- the spatial constraint CellChat actually used ----
  # computeCommunProb stores the cell-group distance matrix in object@images$distance
  # (modeling.R:322), NOT in options$parameter; adj.contact is not persisted by CellChat.
  par <- object@options$parameter
  if (!is.null(object@images$distance))
    utils::write.csv(as.matrix(object@images$distance), p("%s_d_spatial.csv", cond))

  # ---- group sizes, overall and per sample ----
  gs <- data.frame(cell_type = levels(object@idents), n_cells = as.numeric(table(object@idents)),
                   stringsAsFactors = FALSE)
  utils::write.csv(gs, p("%s_group_sizes.csv", cond), row.names = FALSE)
  if (!is.null(object@meta$samples)) {
    utils::write.csv(as.data.frame.matrix(table(object@idents, object@meta$samples)),
                     p("%s_group_sizes_by_sample.csv", cond))
  }

  # ---- the exact DB used, and every inference parameter ----
  utils::write.csv(object@DB$interaction, p("%s_db_used.csv", cond), row.names = FALSE)
  writeLines(jsonlite::toJSON(c(par, list(run.time.seconds = object@options$run.time)),
                              auto_unbox = TRUE, null = "null", na = "null", pretty = TRUE),
             p("%s_options.json", cond))

  log(sprintf("[quant] %s: %d cell types, %d LR pairs tested, %d significant LR links, %d pathways",
              cond, nlevels(object@idents), dim(prob)[3],
              if (is.null(df.net)) 0L else nrow(df.net), length(object@netP$pathways)))
  invisible(list(n_celltypes = nlevels(object@idents), n_lr_tested = dim(prob)[3],
                 n_significant = if (is.null(df.net)) 0L else nrow(df.net),
                 n_pathways = length(object@netP$pathways)))
}

# ---- comparison-level outputs (merged object) ----
save_comparison_quant <- function(cellchat, object.list, out_dir, log = message) {
  qdir <- file.path(out_dir, "quant"); odir <- file.path(out_dir, "objects")
  dir.create(qdir, showWarnings = FALSE, recursive = TRUE)
  dir.create(odir, showWarnings = FALSE, recursive = TRUE)
  p <- function(...) file.path(qdir, sprintf(...))

  saveRDS(cellchat, file.path(odir, "merged.rds"), compress = TRUE)
  saveRDS(object.list, file.path(odir, "object_list.rds"), compress = TRUE)

  nm <- names(object.list)
  # differential count / weight, second minus first
  d.count  <- object.list[[2]]@net$count  - object.list[[1]]@net$count
  d.weight <- object.list[[2]]@net$weight - object.list[[1]]@net$weight
  utils::write.csv(d.count,  p("diff_count_%s_minus_%s.csv",  nm[2], nm[1]))
  utils::write.csv(d.weight, p("diff_weight_%s_minus_%s.csv", nm[2], nm[1]))

  # total interactions / strength per condition
  utils::write.csv(data.frame(
    condition = nm,
    n_interactions = sapply(object.list, function(x) sum(x@net$count)),
    total_strength = sapply(object.list, function(x) sum(x@net$weight)),
    n_pathways = sapply(object.list, function(x) length(x@netP$pathways)),
    stringsAsFactors = FALSE), p("summary_by_condition.csv"), row.names = FALSE)

  # information flow per pathway per condition (what rankNet plots)
  flow <- do.call(rbind, lapply(nm, function(n) {
    o <- object.list[[n]]
    if (is.null(o@netP$prob)) return(NULL)
    data.frame(condition = n, pathway_name = dimnames(o@netP$prob)[[3]],
               information_flow = apply(o@netP$prob, 3, sum), stringsAsFactors = FALSE)
  }))
  if (!is.null(flow)) utils::write.csv(flow, p("information_flow_by_pathway.csv"), row.names = FALSE)

  log(sprintf("[quant] comparison: %s", paste(sprintf("%s=%d links", nm,
              sapply(object.list, function(x) sum(x@net$count))), collapse = ", ")))
  invisible(NULL)
}

# Save a ggplot/ComplexHeatmap/base plot in png + pdf + svg through ONE saver (CLAUDE.md).
# CellChat mixes three plot idioms -- base graphics (netVisual_circle, netVisual_aggregate,
# netVisual_chord_*), ggplot (netVisual_bubble, rankNet, netAnalysis_contribution,
# spatialFeaturePlot) and ComplexHeatmap (netVisual_heatmap) -- so the return value is
# dispatched on rather than declared by the caller.
save_all_formats <- function(plot_expr, path_no_ext, width = 7, height = 5) {
  dir.create(dirname(path_no_ext), showWarnings = FALSE, recursive = TRUE)
  ok <- TRUE
  for (fmt in c("png", "pdf", "svg")) {
    dev_ok <- tryCatch({
      switch(fmt,
        png = grDevices::png(paste0(path_no_ext, ".png"), width = width, height = height,
                             units = "in", res = 300),
        pdf = grDevices::pdf(paste0(path_no_ext, ".pdf"), width = width, height = height),
        svg = grDevices::svg(paste0(path_no_ext, ".svg"), width = width, height = height))
      TRUE
    }, error = function(e) FALSE)
    if (!dev_ok) { ok <- FALSE; next }
    res <- tryCatch({
      p <- plot_expr()
      # ggplot / ComplexHeatmap objects need an explicit print/draw; base plots self-draw
      if (inherits(p, c("Heatmap", "HeatmapList"))) ComplexHeatmap::draw(p)
      else if (inherits(p, c("ggplot", "gg", "patchwork"))) print(p)
      TRUE
    }, error = function(e) { message(sprintf("  [plot-fail] %s (%s): %s",
                                             basename(path_no_ext), fmt, conditionMessage(e))); FALSE })
    grDevices::dev.off()
    if (!res) { ok <- FALSE; unlink(paste0(path_no_ext, ".", fmt)) }
  }
  invisible(ok)
}
