#!/usr/bin/env Rscript
# plot_cellchat.R -- NOTES.md stages E (per-condition visualization) and F (high-vs-low
# comparison). Re-reads the objects written by run_cellchat.R, so replotting never redoes the
# 100-permutation test.
#
# Every plot goes through save_all_formats() -> png + pdf + svg (CLAUDE.md plotting rule).
# Per-LR plots land in TWO separate trees (skill invariant):
#   plots/top_lr/        -- the LRs CellChat itself ranks highest (by communication probability)
#   plots/requested_lr/  -- always GRN_SORT1 and ANXA1_FPR1, whatever their rank
#
#   source scripts/comparators/cellchat/activate_env.sh
#   Rscript scripts/comparators/cellchat/plot_cellchat.R \
#     --out-dir results/comparators/cellchat/GBM/default --conditions low,high
suppressWarnings(suppressMessages({
  library(CellChat); library(patchwork); library(ComplexHeatmap); library(jsonlite)
  # VB:468-469 loads these two explicitly before the communication-pattern section, and it
  # matters: netAnalysis_river builds a ggalluvial plot with stat = "stratum", which fails with
  # "Can't find stat called 'stratum'" unless ggalluvial is ATTACHED (ggalluvial:: is not enough).
  library(NMF); library(ggalluvial)
}))
options(stringsAsFactors = FALSE)

SCRIPT_DIR <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(SCRIPT_DIR, "cellchat_io.R"))

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(f, d = NULL) { i <- which(args == f); if (!length(i)) d else args[i[1] + 1] }
out_dir    <- getarg("--out-dir")
conditions <- strsplit(getarg("--conditions", "low,high"), ",")[[1]]
top_n_lr   <- as.numeric(getarg("--top-n-lr", "10"))
spatial_cores <- getarg("--spatial-cores", "all")   # all | first
max_pathways_spatial <- as.numeric(getarg("--max-pathways-spatial", "0"))  # 0 = no cap
stopifnot(!is.null(out_dir))

log_file <- file.path(out_dir, "plot.log")
log <- function(...) {
  msg <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...))
  cat(msg, "\n"); cat(msg, "\n", file = log_file, append = TRUE)
}
P <- function(...) file.path(out_dir, "plots", ...)
REQUESTED_LR <- c("GRN_SORT1", "ANXA1_FPR1")
dropped <- character(0)   # anything we deliberately do not draw gets recorded, never silent

# Wrap the saver so a plot that ERRORS is recorded too, not merely logged. Without this a
# failed render leaves no trace in plots_not_produced.txt and the run looks complete.
save_tracked <- function(expr, path_no_ext, width = 7, height = 5, what = NULL) {
  ok <- save_all_formats(expr, path_no_ext, width, height)
  if (!isTRUE(ok)) {
    dropped <<- c(dropped, sprintf("FAILED: %s", what %||% sub(paste0("^", out_dir, "/plots/"), "",
                                                               path_no_ext)))
  }
  ok
}
`%||%` <- function(a, b) if (is.null(a)) b else a

object.list <- lapply(conditions, function(c) readRDS(file.path(out_dir, "objects", paste0(c, ".rds"))))
names(object.list) <- conditions
log("loaded ", paste(conditions, collapse = ", "))

# ============================================================ Stage E: per condition
for (cond in conditions) {
  obj <- object.list[[cond]]
  ct  <- levels(obj@idents)
  nct <- length(ct)
  pathways <- obj@netP$pathways
  samples  <- levels(obj@meta$samples)
  log(cond, ": ", nct, " cell types, ", length(pathways), " pathways, ", length(samples), " samples")
  groupSize <- as.numeric(table(obj@idents))

  # ---- E1/E2: aggregated network ----
  A <- function(...) P(cond, "aggregate", ...)
  save_all_formats(function() netVisual_circle(obj@net$count, vertex.weight = rowSums(obj@net$count),
      weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions"),
      A("circle_count"), 7, 7)
  save_all_formats(function() netVisual_circle(obj@net$weight, vertex.weight = rowSums(obj@net$weight),
      weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction weights/strength"),
      A("circle_weight"), 7, 7)
  save_all_formats(function() netVisual_heatmap(obj, measure = "count", color.heatmap = "Blues"),
      A("heatmap_count"), 6, 5)
  save_all_formats(function() netVisual_heatmap(obj, measure = "weight", color.heatmap = "Blues"),
      A("heatmap_weight"), 6, 5)

  # ---- E3: one circle per sender ----
  mat <- obj@net$weight
  for (i in seq_len(nrow(mat))) {
    save_all_formats(function() {
      m2 <- matrix(0, nrow(mat), ncol(mat), dimnames = dimnames(mat)); m2[i, ] <- mat[i, ]
      netVisual_circle(m2, vertex.weight = groupSize, weight.scale = TRUE,
                       edge.weight.max = max(mat), title.name = rownames(mat)[i])
    }, A(sprintf("sender_%02d_%s", i, gsub("[^A-Za-z0-9]+", "_", rownames(mat)[i]))), 6, 6)
  }

  # ---- E12/E13: bubble + chord over all sources/targets ----
  save_all_formats(function() netVisual_bubble(obj, sources.use = seq_len(nct),
      targets.use = seq_len(nct), remove.isolate = FALSE),
      A("bubble_all"), max(8, nct * 1.6), max(6, 0.22 * nrow(obj@LR$LRsig)))
  # All-sources x all-targets is attempted, but on this data it is 161 ligand/receptor sectors
  # and circlize cannot lay that out at any small.gap. The tutorial's own scope is ONE sender
  # (VB:381 `sources.use = 4`), so per-sender chords are the covering set and are always drawn.
  save_tracked(function() netVisual_chord_gene(obj, sources.use = seq_len(nct),
      targets.use = seq_len(nct), lab.cex = 0.4, legend.pos.x = 10), A("chord_gene_all"), 10, 10,
      what = sprintf("%s chord_gene_all (all sources x all targets; circlize sector limit)", cond))
  save_tracked(function() netVisual_chord_gene(obj, sources.use = seq_len(nct),
      targets.use = seq_len(nct), slot.name = "netP", legend.pos.x = 10),
      A("chord_pathway_all"), 10, 10, what = sprintf("%s chord_pathway_all", cond))
  for (i in seq_len(nct)) {
    save_tracked(function() netVisual_chord_gene(obj, sources.use = i,
        targets.use = seq_len(nct), lab.cex = 0.4, legend.pos.x = 10,
        title.name = paste0("Signaling from ", ct[i], " - ", cond)),
        A(paste0("chord_gene_sender_", gsub("[^A-Za-z0-9]+", "_", ct[i]))), 10, 10,
        what = sprintf("%s chord_gene_sender_%s", cond, ct[i]))
  }

  # ---- E4-E7, E10, E14, E15, E22: per pathway ----
  vertex.receiver <- seq_len(min(4, nct))   # left panel of the hierarchy plot
  for (pw in pathways) {
    S <- function(...) P(cond, "pathways", gsub("[^A-Za-z0-9]+", "_", pw), ...)
    save_all_formats(function() netVisual_aggregate(obj, signaling = pw, layout = "circle"),
                     S("circle"), 7, 7)
    save_all_formats(function() netVisual_aggregate(obj, signaling = pw, layout = "chord"),
                     S("chord"), 8, 8)
    save_all_formats(function() netVisual_aggregate(obj, signaling = pw,
                     vertex.receiver = vertex.receiver), S("hierarchy"), 10, 6)
    save_all_formats(function() netVisual_heatmap(obj, signaling = pw, color.heatmap = "Reds"),
                     S("heatmap"), 6, 5)
    save_all_formats(function() netAnalysis_contribution(obj, signaling = pw),
                     S("LR_contribution"), 4, 3)
    save_all_formats(function() netAnalysis_signalingRole_network(obj, signaling = pw,
                     width = 8, height = 2.5, font.size = 10), S("signalingRole_network"), 8, 3)
    grp <- setNames(ct, ct)
    save_all_formats(function() netVisual_chord_cell(obj, signaling = pw, group = grp,
                     title.name = paste0(pw, " signaling network")), S("chord_cell"), 8, 8)
    save_all_formats(function() plotGeneExpression(obj, signaling = pw, enriched.only = TRUE,
                     type = "violin"), S("geneExpression_violin"), 7, 5)
    # E11: every enriched L-R pair of this pathway, circle + chord (CellChat's own ranking)
    lrs <- tryCatch(extractEnrichedLR(obj, signaling = pw, geneLR.return = FALSE),
                    error = function(e) NULL)
    if (!is.null(lrs) && nrow(lrs) > 0) {
      for (lr in as.character(lrs[[1]])) {
        L <- function(...) S("LR", gsub("[^A-Za-z0-9]+", "_", lr), ...)
        save_all_formats(function() netVisual_individual(obj, signaling = pw,
            pairLR.use = lr, layout = "circle"), L("circle"), 7, 7)
        save_all_formats(function() netVisual_individual(obj, signaling = pw,
            pairLR.use = lr, layout = "chord"), L("chord"), 8, 8)
      }
    }
  }

  # ---- E8/E9: spatial network plots, per pathway x per core ----
  cores_use <- if (spatial_cores == "first") samples[1] else samples
  pw_spatial <- if (max_pathways_spatial > 0 && length(pathways) > max_pathways_spatial)
    pathways[seq_len(max_pathways_spatial)] else pathways
  if (length(pw_spatial) < length(pathways))
    dropped <- c(dropped, sprintf("%s: spatial plots limited to %d/%d pathways",
                                  cond, length(pw_spatial), length(pathways)))
  if (length(cores_use) < length(samples))
    dropped <- c(dropped, sprintf("%s: spatial plots limited to core %s of %d",
                                  cond, cores_use[1], length(samples)))
  for (pw in pw_spatial) for (s in cores_use) {
    SP <- function(...) P(cond, "spatial", paste0("core", s), gsub("[^A-Za-z0-9]+", "_", pw), ...)
    save_all_formats(function() netVisual_aggregate(obj, signaling = pw, sample.use = s,
        layout = "spatial", edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.2,
        vertex.label.cex = 0), SP("spatial"), 8, 8)
    save_all_formats(function() netVisual_aggregate(obj, signaling = pw, sample.use = s,
        layout = "spatial", edge.width.max = 2, alpha.image = 0.2, vertex.weight = "incoming",
        vertex.size.max = 6, vertex.label.cex = 0), SP("spatial_incoming"), 8, 8)
  }

  # ---- E16-E19: systems analysis ----
  Y <- function(...) P(cond, "systems", ...)
  save_all_formats(function() netAnalysis_signalingRole_scatter(obj), Y("role_scatter_all"), 6, 5)
  save_all_formats(function() netAnalysis_signalingRole_heatmap(obj, pattern = "outgoing"),
                   Y("role_heatmap_outgoing"), 7, 8)
  save_all_formats(function() netAnalysis_signalingRole_heatmap(obj, pattern = "incoming",
                   color.heatmap = "GnBu"), Y("role_heatmap_incoming"), 7, 8)
  save_all_formats(function() netAnalysis_signalingRole_heatmap(obj, pattern = "all",
                   color.heatmap = "OrRd"), Y("role_heatmap_all"), 7, 8)

  # E18: communication patterns. selectK first, then k from the cophenetic/silhouette drop.
  # identifyCommunicationPatterns REQUIRES k (analysis.R:385-387 stop()s on NULL) and the
  # vignette picks it by eye off the selectK curve: "a suitable number of patterns is the one
  # at which Cophenetic and Silhouette values begin to drop suddenly" (VB:455). pick_k applies
  # exactly that rule to the same NMF::nmfEstimateRank measures selectK plots, and the measures
  # are written to disk so the choice is auditable rather than eyeballed.
  pick_k <- function(meas) {
    co <- meas$cophenetic; ks <- meas$rank
    i <- which(diff(co) < -0.01)[1]
    if (is.na(i)) ks[which.max(co)] else ks[i]
  }
  for (pat in c("outgoing", "incoming")) {
    # NMF rank must stay below both matrix dimensions; selectK's default seq(2,10) errors out
    # ("All the runs produced an error") whenever a condition has <=10 significant pathways.
    prob <- obj@netP$prob
    d0 <- if (pat == "outgoing") apply(prob, c(1, 3), sum) else apply(prob, c(2, 3), sum)
    d0 <- sweep(d0, 2L, apply(d0, 2, function(x) max(x, na.rm = TRUE)), "/", check.margin = FALSE)
    d0 <- as.matrix(d0); d0 <- d0[rowSums(d0) != 0, , drop = FALSE]
    kmax <- min(10, nrow(d0) - 1, ncol(d0) - 1)
    if (kmax < 2) {
      log("  [skip] patterns ", pat, ": matrix is ", nrow(d0), "x", ncol(d0), " -- no valid rank")
      dropped <- c(dropped, sprintf("%s: communication patterns (%s) -- only %d pathways, NMF needs >=3",
                                    cond, pat, ncol(d0)))
      next
    }
    krange <- seq(2, kmax)
    if (kmax < 10) log("  ", cond, " ", pat, ": k.range capped at ", kmax,
                       " (", nrow(d0), " cell types x ", ncol(d0), " pathways)")
    # CellChat's selectK() calls NMF::nmfEstimateRank with NMF's default (parallel) foreach
    # backend and exposes no way to override it. In this plotting session that backend is
    # already poisoned by earlier Seurat/ComplexHeatmap work and every NMF run dies with
    # "All the runs produced an error" -- the identical call succeeds in a clean session.
    # So the measures are computed here with .pbackend = "seq" (deterministic, same numbers)
    # and the selectK curve is redrawn from them; identifyCommunicationPatterns below is a
    # single-run NMF and is unaffected.
    kk <- tryCatch({
      est <- NMF::nmfEstimateRank(d0, range = krange, method = "lee", nrun = 30, seed = 10,
                                  .pbackend = "seq")
      mdf <- rbind(
        data.frame(k = est$measures$rank, score = est$measures$cophenetic, Measure = "Cophenetic"),
        data.frame(k = est$measures$rank, score = est$measures$silhouette.consensus,
                   Measure = "Silhouette"))
      save_all_formats(function() ggplot2::ggplot(mdf, ggplot2::aes(k, score)) +
          ggplot2::geom_line() + ggplot2::geom_point() +
          ggplot2::facet_wrap(~Measure, scales = "free_y") +
          ggplot2::scale_x_continuous(breaks = krange) +
          ggplot2::labs(title = paste0(pat, " signaling"), x = "Number of patterns (k)", y = "") +
          ggplot2::theme_bw(), Y(paste0("selectK_", pat)), 7, 3)
      utils::write.csv(est$measures, file.path(out_dir, "quant",
                       sprintf("%s_selectK_%s_measures.csv", cond, pat)), row.names = FALSE)
      k_sel <- pick_k(est$measures)
      log("  ", cond, " ", pat, " patterns: k = ", k_sel, " (first cophenetic drop)")
      obj2 <- identifyCommunicationPatterns(obj, pattern = pat, k = k_sel)
      save_all_formats(function() netAnalysis_river(obj2, pattern = pat),
                       Y(paste0("river_", pat)), 9, 6)
      save_all_formats(function() netAnalysis_dot(obj2, pattern = pat),
                       Y(paste0("dot_", pat)), 7, 4)
      utils::write.csv(obj2@netP$pattern[[pat]]$pattern$cell,
                       file.path(out_dir, "quant", sprintf("%s_pattern_%s_cell.csv", cond, pat)),
                       row.names = FALSE)
      utils::write.csv(obj2@netP$pattern[[pat]]$pattern$signaling,
                       file.path(out_dir, "quant", sprintf("%s_pattern_%s_signaling.csv", cond, pat)),
                       row.names = FALSE)
      k_sel
    }, error = function(e) { log("  [skip] patterns ", pat, ": ", conditionMessage(e)); NA })
    if (is.na(kk)) dropped <- c(dropped, sprintf("%s: communication patterns (%s) failed", cond, pat))
  }

  # E19: manifold learning of the pathway networks
  for (ty in c("functional", "structural")) {
    tryCatch({
      o <- computeNetSimilarity(obj, type = ty)
      o <- netEmbedding(o, type = ty, umap.method = "uwot")
      o <- netClustering(o, type = ty, do.parallel = FALSE)
      save_all_formats(function() netVisual_embedding(o, type = ty, label.size = 3.5),
                       Y(paste0("embedding_", ty)), 7, 6)
      utils::write.csv(data.frame(pathway = rownames(o@netP$similarity[[ty]]$dr[["single"]]),
                                  o@netP$similarity[[ty]]$dr[["single"]],
                                  group = o@netP$similarity[[ty]]$group[["single"]]),
                       file.path(out_dir, "quant", sprintf("%s_embedding_%s.csv", cond, ty)),
                       row.names = FALSE)
    }, error = function(e) {
      log("  [skip] embedding ", ty, ": ", conditionMessage(e))
      dropped <<- c(dropped, sprintf("%s: %s manifold embedding failed (%s)", cond, ty,
                                     conditionMessage(e)))
    })
  }

  # ---- top_lr: the LRs CellChat itself ranks highest ----
  df <- tryCatch(subsetCommunication(obj), error = function(e) NULL)
  if (!is.null(df) && nrow(df) > 0) {
    agg <- aggregate(prob ~ interaction_name, data = df, FUN = sum)
    agg <- agg[order(-agg$prob), ]
    top <- utils::head(as.character(agg$interaction_name), top_n_lr)
    utils::write.csv(agg, file.path(out_dir, "quant", sprintf("%s_lr_ranked.csv", cond)),
                     row.names = FALSE)
    for (lr in top) {
      pw <- df$pathway_name[match(lr, df$interaction_name)]
      T_ <- function(...) P("top_lr", cond, gsub("[^A-Za-z0-9]+", "_", lr), ...)
      save_all_formats(function() netVisual_individual(obj, signaling = pw, pairLR.use = lr,
                       layout = "circle"), T_("circle"), 7, 7)
      save_all_formats(function() netVisual_individual(obj, signaling = pw, pairLR.use = lr,
                       layout = "chord"), T_("chord"), 8, 8)
      save_all_formats(function() netVisual_bubble(obj, sources.use = seq_len(nct),
                       targets.use = seq_len(nct),
                       pairLR.use = data.frame(interaction_name = lr), remove.isolate = FALSE),
                       T_("bubble"), 8, 3)
      for (s in cores_use) {
        save_all_formats(function() spatialFeaturePlot(obj, pairLR.use = lr, sample.use = s,
            point.size = 0.5, do.binary = FALSE, cutoff = 0.05, enriched.only = FALSE,
            color.heatmap = "Reds", direction = 1), T_(paste0("spatial_core", s)), 7, 6)
      }
    }
  }

  # ---- requested_lr: always, whatever the rank ----
  # An absent / untested / non-significant requested LR IS the result (skill invariant), so the
  # status is written to requested_lr_status.csv either way and the network plots -- which
  # CellChat refuses to draw without a significant link -- are skipped with a reason.
  req_status <- list()
  for (lr in REQUESTED_LR) {
    R_ <- function(...) P("requested_lr", cond, gsub("[^A-Za-z0-9]+", "_", lr), ...)
    in_db  <- lr %in% obj@DB$interaction$interaction_name
    tested <- lr %in% obj@LR$LRsig$interaction_name
    pw <- if (in_db) obj@DB$interaction$pathway_name[match(lr, obj@DB$interaction$interaction_name)] else NA
    sig_pairs <- if (!is.null(df) && nrow(df) > 0) sum(df$interaction_name == lr) else 0L
    pathway_sig <- !is.na(pw) && pw %in% obj@netP$pathways
    req_status[[lr]] <- data.frame(condition = cond, interaction_name = lr, pathway = pw,
        in_db = in_db, tested = tested, n_significant_pairs = sig_pairs,
        pathway_significant = pathway_sig, stringsAsFactors = FALSE)
    log("  requested ", lr, " (", pw, "): in_db=", in_db, " tested=", tested,
        " significant_pairs=", sig_pairs, " pathway_significant=", pathway_sig)
    if (!tested) next
    genes <- unlist(strsplit(gsub("_", " ", lr), " "))
    # network views need a significant link; CellChat errors with "There is no significant
    # communication of <pathway>" otherwise. Gene/LR spatial maps below do NOT need one.
    if (sig_pairs > 0 && pathway_sig) {
      save_all_formats(function() netVisual_individual(obj, signaling = pw, pairLR.use = lr,
                       layout = "circle"), R_("circle"), 7, 7)
      save_all_formats(function() netVisual_individual(obj, signaling = pw, pairLR.use = lr,
                       layout = "chord"), R_("chord"), 8, 8)
      save_all_formats(function() netVisual_bubble(obj, sources.use = seq_len(nct),
                       targets.use = seq_len(nct),
                       pairLR.use = data.frame(interaction_name = lr), remove.isolate = FALSE),
                       R_("bubble"), 8, 3)
      save_all_formats(function() netAnalysis_contribution(obj, signaling = pw),
                       R_("pathway_LR_contribution"), 4, 3)
    } else {
      msg <- sprintf("%s: %s (%s) has no significant communication -- network plots skipped, spatial maps still written",
                     cond, lr, pw)
      log("  [note] ", msg); dropped <- c(dropped, msg)
    }
    for (s in cores_use) {
      save_all_formats(function() spatialFeaturePlot(obj, features = genes, sample.use = s,
          point.size = 0.8, color.heatmap = "Reds", direction = 1),
          R_(paste0("spatial_genes_core", s)), 9, 5)
      save_all_formats(function() spatialFeaturePlot(obj, pairLR.use = lr, sample.use = s,
          point.size = 0.5, do.binary = FALSE, cutoff = 0.05, enriched.only = FALSE,
          color.heatmap = "Reds", direction = 1), R_(paste0("spatial_LR_core", s)), 7, 6)
      save_all_formats(function() spatialFeaturePlot(obj, pairLR.use = lr, sample.use = s,
          point.size = 1, do.binary = TRUE, cutoff = 0.05, enriched.only = FALSE,
          color.heatmap = "Reds", direction = 1), R_(paste0("spatial_LR_binary_core", s)), 7, 6)
    }
    # F21 idiom, applied per condition
    tryCatch({
      d <- findEnrichedSignaling(obj, features = genes[1], idents = ct, pattern = "outgoing")
      utils::write.csv(d, file.path(out_dir, "quant",
                       sprintf("%s_findEnrichedSignaling_%s.csv", cond, genes[1])), row.names = FALSE)
    }, error = function(e) log("  [skip] findEnrichedSignaling ", genes[1], ": ", conditionMessage(e)))
  }
  if (length(req_status))
    utils::write.csv(do.call(rbind, req_status),
                     file.path(out_dir, "quant", sprintf("%s_requested_lr_plot_status.csv", cond)),
                     row.names = FALSE)
  log(cond, ": plots done")
}

# ============================================================ Stage F: comparison
if (length(conditions) == 2) {
  C <- function(...) P("comparison", ...)
  cellchat <- mergeCellChat(object.list, add.names = names(object.list))       # F1
  save_comparison_quant(cellchat, object.list, out_dir, log = log)
  nm <- names(object.list)
  log("merged: ", paste(nm, collapse = " -> "), " (positive change = higher in ", nm[2], ")")

  save_all_formats(function() compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2)) +
      compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2), measure = "weight"),
      C("compareInteractions"), 6, 3)                                          # F2
  save_all_formats(function() netVisual_diffInteraction(cellchat, weight.scale = TRUE),
      C("diffInteraction_count"), 7, 7)                                        # F3
  save_all_formats(function() netVisual_diffInteraction(cellchat, weight.scale = TRUE,
      measure = "weight"), C("diffInteraction_weight"), 7, 7)
  save_all_formats(function() netVisual_heatmap(cellchat), C("diff_heatmap_count"), 6, 5)  # F4
  save_all_formats(function() netVisual_heatmap(cellchat, measure = "weight"),
      C("diff_heatmap_weight"), 6, 5)

  weight.max <- getMaxWeight(object.list, attribute = c("idents", "count"))    # F5
  for (i in seq_along(object.list)) {
    save_all_formats(function() netVisual_circle(object.list[[i]]@net$count, weight.scale = TRUE,
        label.edge = FALSE, edge.weight.max = weight.max[2], edge.width.max = 12,
        title.name = paste0("Number of interactions - ", nm[i])),
        C(paste0("circle_count_", nm[i])), 7, 7)
  }

  num.link <- sapply(object.list, function(x) rowSums(x@net$count) + colSums(x@net$count) -
                       diag(x@net$count))                                       # F7
  wMinMax <- c(min(num.link), max(num.link))
  save_all_formats(function() patchwork::wrap_plots(plots = lapply(seq_along(object.list),
      function(i) netAnalysis_signalingRole_scatter(object.list[[i]], title = nm[i],
      weight.MinMax = wMinMax))), C("role_scatter_both"), 10, 5)

  for (ident in levels(object.list[[1]]@idents)) {                              # F8, all 9
    save_all_formats(function() netAnalysis_signalingChanges_scatter(cellchat, idents.use = ident),
        C("signalingChanges", gsub("[^A-Za-z0-9]+", "_", ident)), 6, 5)
  }

  for (ty in c("functional", "structural")) {                                   # F9/F10
    tryCatch({
      cc <- computeNetSimilarityPairwise(cellchat, type = ty)
      cc <- netEmbedding(cc, type = ty, umap.method = "uwot")
      cc <- netClustering(cc, type = ty, do.parallel = FALSE)
      save_all_formats(function() netVisual_embeddingPairwise(cc, type = ty, label.size = 3.5),
          C(paste0("embeddingPairwise_", ty)), 8, 7)
      save_all_formats(function() rankSimilarity(cc, type = ty),
          C(paste0("rankSimilarity_", ty)), 4, 6)
    }, error = function(e) {
      log("  [skip] pairwise ", ty, ": ", conditionMessage(e))
      dropped <<- c(dropped, sprintf("comparison: %s pairwise manifold failed (%s)", ty,
                                     conditionMessage(e)))
    })
  }

  save_all_formats(function() rankNet(cellchat, mode = "comparison", measure = "weight",
      stacked = TRUE, do.stat = TRUE), C("rankNet_stacked"), 6, 9)              # F11
  save_all_formats(function() rankNet(cellchat, mode = "comparison", measure = "weight",
      stacked = FALSE, do.stat = TRUE), C("rankNet_unstacked"), 6, 9)

  pathway.union <- union(object.list[[1]]@netP$pathways, object.list[[2]]@netP$pathways)
  for (pat in c("outgoing", "incoming", "all")) {                               # F12
    ch <- c(outgoing = "BuGn", incoming = "GnBu", all = "OrRd")[[pat]]
    save_all_formats(function() {
      h1 <- netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = pat,
              signaling = pathway.union, title = nm[1], width = 5, height = 9, color.heatmap = ch)
      h2 <- netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = pat,
              signaling = pathway.union, title = nm[2], width = 5, height = 9, color.heatmap = ch)
      h1 + h2
    }, C(paste0("role_heatmap_", pat, "_both")), 12, 10)
  }

  nct <- nlevels(object.list[[1]]@idents)
  save_all_formats(function() netVisual_bubble(cellchat, sources.use = seq_len(nct),
      targets.use = seq_len(nct), comparison = c(1, 2), angle.x = 45),
      C("bubble_comparison"), 14, 14)                                           # F13
  save_all_formats(function() netVisual_bubble(cellchat, sources.use = seq_len(nct),
      targets.use = seq_len(nct), comparison = c(1, 2), max.dataset = 2,
      title.name = paste0("Increased signaling in ", nm[2]), angle.x = 45, remove.isolate = TRUE),
      C("bubble_increased"), 14, 12)
  save_all_formats(function() netVisual_bubble(cellchat, sources.use = seq_len(nct),
      targets.use = seq_len(nct), comparison = c(1, 2), max.dataset = 1,
      title.name = paste0("Decreased signaling in ", nm[2]), angle.x = 45, remove.isolate = TRUE),
      C("bubble_decreased"), 14, 12)

  # ---- F14-F19: DEA-based up/down signaling ----
  tryCatch({
    pos.dataset <- nm[2]; features.name <- paste0(pos.dataset, ".merged")
    cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets",
        pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE,
        thresh.pc = 0.1, thresh.fc = 0.05, thresh.p = 0.05, group.DE.combined = FALSE)
    net <- netMappingDEG(cellchat, features.name = features.name, variable.all = TRUE)
    net.up   <- subsetCommunication(cellchat, net = net, datasets = nm[2],
                                    ligand.logFC = 0.05, receptor.logFC = NULL)
    net.down <- subsetCommunication(cellchat, net = net, datasets = nm[1],
                                    ligand.logFC = -0.05, receptor.logFC = NULL)
    utils::write.csv(net,      file.path(out_dir, "quant", "netMappingDEG.csv"), row.names = FALSE)
    utils::write.csv(net.up,   file.path(out_dir, "quant", sprintf("net_up_in_%s.csv", nm[2])),
                     row.names = FALSE)
    utils::write.csv(net.down, file.path(out_dir, "quant", sprintf("net_down_in_%s.csv", nm[2])),
                     row.names = FALSE)
    gene.up   <- extractGeneSubsetFromPair(net.up, cellchat)                    # F16
    gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
    writeLines(as.character(gene.up),   file.path(out_dir, "quant", "genes_up.txt"))
    writeLines(as.character(gene.down), file.path(out_dir, "quant", "genes_down.txt"))
    log("DEA: ", nrow(net.up), " up / ", nrow(net.down), " down LR pairs in ", nm[2])

    if (nrow(net.up) > 0) save_all_formats(function() netVisual_bubble(cellchat,
        pairLR.use = net.up[, "interaction_name", drop = FALSE], sources.use = seq_len(nct),
        targets.use = seq_len(nct), comparison = c(1, 2), angle.x = 90, remove.isolate = TRUE,
        title.name = paste0("Up-regulated signaling in ", nm[2])), C("bubble_up"), 14, 12)
    if (nrow(net.down) > 0) save_all_formats(function() netVisual_bubble(cellchat,
        pairLR.use = net.down[, "interaction_name", drop = FALSE], sources.use = seq_len(nct),
        targets.use = seq_len(nct), comparison = c(1, 2), angle.x = 90, remove.isolate = TRUE,
        title.name = paste0("Down-regulated signaling in ", nm[2])), C("bubble_down"), 14, 12)
    if (nrow(net.up) > 0) save_all_formats(function() netVisual_chord_gene(object.list[[2]],
        sources.use = seq_len(nct), targets.use = seq_len(nct), slot.name = "net", net = net.up,
        lab.cex = 0.6, small.gap = 3.5, title.name = paste0("Up-regulated in ", nm[2])),
        C("chord_up"), 10, 10)
    if (nrow(net.down) > 0) save_all_formats(function() netVisual_chord_gene(object.list[[1]],
        sources.use = seq_len(nct), targets.use = seq_len(nct), slot.name = "net", net = net.down,
        lab.cex = 0.6, small.gap = 3.5, title.name = paste0("Down-regulated in ", nm[2])),
        C("chord_down"), 10, 10)
    for (d in list(list(n = net.up, t = "up"), list(n = net.down, t = "down"))) {  # F19
      if (nrow(d$n) == 0) next
      save_all_formats(function() computeEnrichmentScore(d$n, species = "human",
          variable.both = TRUE), C(paste0("wordcloud_", d$t)), 5, 5)
    }
    # where do the requested LRIs sit in the up/down lists?
    req <- do.call(rbind, lapply(REQUESTED_LR, function(lr) data.frame(
      interaction_name = lr,
      in_net_up   = lr %in% net.up$interaction_name,
      in_net_down = lr %in% net.down$interaction_name, stringsAsFactors = FALSE)))
    utils::write.csv(req, file.path(out_dir, "quant", "requested_lr_in_DEA.csv"), row.names = FALSE)
  }, error = function(e) {
    log("  [skip] DEA comparison: ", conditionMessage(e))
    dropped <<- c(dropped, sprintf("comparison: DEA up/down failed (%s)", conditionMessage(e)))
  })

  # F20: signaling gene expression split by condition, for the requested LRIs' pathways
  cellchat@meta$datasets <- factor(cellchat@meta$datasets, levels = nm)
  req_pw <- unique(object.list[[1]]@DB$interaction$pathway_name[
    match(REQUESTED_LR, object.list[[1]]@DB$interaction$interaction_name)])
  for (pw in req_pw[!is.na(req_pw)]) {
    save_all_formats(function() plotGeneExpression(cellchat, signaling = pw,
        split.by = "datasets", colors.ggplot = TRUE, type = "violin"),
        C("geneExpression", gsub("[^A-Za-z0-9]+", "_", pw)), 8, 5)
  }
} else {
  log("only ", length(conditions), " condition(s) -- Stage F comparison skipped")
  dropped <- c(dropped, "Stage F comparison skipped (needs exactly 2 conditions)")
}

# Always rewrite (or clear) this file: leaving a previous run's copy in place makes a clean
# rerun look like it still had failures.
nf <- file.path(out_dir, "plots_not_produced.txt")
if (length(dropped)) {
  writeLines(dropped, nf)
  log("NOT produced (", length(dropped), "): ", paste(dropped, collapse = " | "))
} else {
  if (file.exists(nf)) unlink(nf)
  log("NOT produced: (none) -- every plot in the contract was written")
}
n_png <- length(list.files(file.path(out_dir, "plots"), pattern = "\\.png$", recursive = TRUE))
log("wrote ", n_png, " plots (x png/pdf/svg) under ", file.path(out_dir, "plots"))
log("=== plotting done ===")
