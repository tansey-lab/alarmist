#!/usr/bin/env Rscript
# Stages B and C of the NICHES workflow: merge the per-core objects, embed, and run the
# high-grade vs low-grade differential -- exactly the calls in vignette 07 (with vignette 04's
# idioms for the CellToCellSpatial object and vignette 01's marker/heatmap block).
# See scripts/comparators/niches/NOTES.md for the call-by-call contract.
#
# Usage:
#   Rscript analyze_niches.R --tier <alra|noimpute> --organization <NeighborhoodToCell|CellToCellSpatial> \
#                            [--root results/comparators/niches/GBM] [--force]
#
# Resumable: the merged+embedded object is cached as merged.rds and reloaded on a rerun, so
# the expensive merge/ScaleData/PCA/UMAP never repeats.
suppressWarnings(suppressMessages({
  library(Seurat); library(SeuratObject); library(Matrix); library(dplyr)
  library(ggplot2); library(cowplot); library(patchwork); library(jsonlite)
}))

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(f, d = NULL) { i <- match(f, args); if (is.na(i)) d else args[i + 1] }
`%||%` <- function(a, b) if (is.null(a)) b else a

tier  <- getarg("--tier")  %||% stop("need --tier alra|noimpute")
org   <- getarg("--organization", "NeighborhoodToCell")
root  <- getarg("--root", "results/comparators/niches/GBM")
force <- "--force" %in% args
stopifnot(tier %in% c("alra", "noimpute"),
          org  %in% c("NeighborhoodToCell", "CellToCellSpatial"))

tier_dir <- file.path(root, "cellchatdb2", tier)
out_dir  <- file.path(tier_dir, "_analysis", org)
plot_dir <- file.path(out_dir, "plots")
for (d in c(out_dir, plot_dir, file.path(plot_dir, "top_lr"), file.path(plot_dir, "requested_lr"),
            file.path(out_dir, "markers"))) dir.create(d, showWarnings = FALSE, recursive = TRUE)
log <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n")

# The two LRs the benchmark always asks for (ALARMIST motif-1 mGAM <-> MES-like loop).
# NICHES joins ligand and receptor with an em-dash (U+2014), not a hyphen.
REQUESTED <- c("GRN—SORT1", "ANXA1—FPR1")

## ---------------- plot saver: png + pdf (+ svg when the panel is not enormous) -----------
# Repo convention is png+pdf+svg through one saver; an svg of a 100k-point UMAP is hundreds
# of MB, so svg is skipped above a point budget and the skip is logged, never silent.
save_plot <- function(p, stem, w = 8, h = 6, npoints = 0, svg_max = 50000) {
  ggsave(paste0(stem, ".png"), p, width = w, height = h, dpi = 200, bg = "white")
  ggsave(paste0(stem, ".pdf"), p, width = w, height = h, device = grDevices::cairo_pdf)
  if (npoints <= svg_max) {
    ggsave(paste0(stem, ".svg"), p, width = w, height = h)
  } else {
    log(sprintf("   (svg skipped for %s: %d points > %d)", basename(stem), npoints, svg_max))
  }
  invisible(NULL)
}

## ---------------- embedding (vignette 07:86-93), with the disp overflow guard ------------
# Vignette 07 does ScaleData -> FindVariableFeatures(selection.method = "disp") -> RunPCA ->
# RunUMAP. On the ALRA tier that FindVariableFeatures call CRASHES, and the reason is a real
# incompatibility rather than anything we did:
#
#   NICHES deliberately stores the RAW ligand x receptor PRODUCTS in the Seurat `data` slot
#   (RunNeighborhoodToCell.R:77-81 copies `counts` into `data` with no normalisation), but
#   Seurat's "disp" method assumes `data` is log1p-space and un-logs it with ExpMean/LogVMR,
#   i.e. exp(x) - 1. Our ALRA-tier products reach 4,933.5; exp(709) already overflows a double,
#   so feature.mean becomes Inf and CalcDispersion dies in
#   seq.int(rx[1], rx[2], length.out = nb) : 'to' must be a finite number.
#
# Vignette 07 itself pairs ALRA with selection.method = "disp", so this is a latent break in
# the tutorial path that only bites at our dynamic range. We keep "disp" wherever it works
# (the whole noimpute tier) and fall back to Seurat's default "vst" -- which never exponentiates
# -- only where "disp" is numerically impossible. The choice is logged and recorded in the
# manifest, never silent. See DEVIATIONS.md D10.
HVF_METHOD <- "disp"
embed_merged <- function(obj) {
  log("ScaleData"); obj <- ScaleData(obj, verbose = FALSE)
  log("FindVariableFeatures(selection.method = 'disp')")
  ok <- tryCatch({ obj <- FindVariableFeatures(obj, selection.method = "disp", verbose = FALSE)
                   TRUE },
                 error = function(e) { log("  !! 'disp' failed:", conditionMessage(e)); FALSE })
  if (!ok) {
    mx <- max(LayerData(obj, assay = org, layer = "data"))
    log(sprintf("  !! max score = %.1f; exp() overflows above 709 -> falling back to 'vst'", mx))
    obj <- FindVariableFeatures(obj, selection.method = "vst", verbose = FALSE)
    HVF_METHOD <<- "vst"
  }
  npcs <- min(100, nrow(obj) - 1, ncol(obj) - 1)
  log("RunPCA(npcs =", npcs, ")"); obj <- RunPCA(obj, npcs = npcs, verbose = FALSE)
  dims <- min(50, npcs)
  log("RunUMAP(dims = 1:", dims, ")"); obj <- RunUMAP(obj, dims = 1:dims, verbose = FALSE)
  # persist the choice on the object so a resume from merged.rds reports it correctly
  obj@misc$hvf_selection_method <- HVF_METHOD
  obj
}

## ---------------- Stage B1-B3: load per-core objects, tag, merge -------------------------
merged_rds <- file.path(out_dir, "merged.rds")
raw_rds    <- file.path(out_dir, "merged_raw.rds")
if (file.exists(merged_rds) && !force) {
  log("loading cached merged object", merged_rds)
  obj <- readRDS(merged_rds)
} else if (file.exists(raw_rds) && !force) {
  # merge succeeded on an earlier attempt but the embedding did not -- resume from there
  log("loading cached pre-embedding merge", raw_rds)
  obj <- readRDS(raw_rds)
  obj <- embed_merged(obj)
  saveRDS(obj, merged_rds); log("cached ->", merged_rds)
} else {
  cores <- list.dirs(tier_dir, recursive = FALSE, full.names = TRUE)
  cores <- cores[file.exists(file.path(cores, "DONE"))]
  cores <- cores[grepl("core[0-9]+$", cores)]
  stopifnot(length(cores) > 0)
  # numeric core order, purely cosmetic
  cores <- cores[order(as.integer(sub(".*core", "", cores)))]
  log("merging", length(cores), "cores:", paste(basename(cores), collapse = ", "))

  objs <- list()
  for (cd in cores) {
    f <- file.path(cd, "objects", paste0(org, ".rds"))
    if (!file.exists(f)) { log("  missing", f, "-- skipping"); next }
    o <- readRDS(f)
    # NOTES B2: tag condition (vignette 07:80). grade is carried per-cell by meta.data.to.map;
    # for CellToCellSpatial it arrives as grade.Receiving (the receiving cell's grade).
    gcol <- if ("grade" %in% colnames(o@meta.data)) "grade" else "grade.Receiving"
    tcol <- if ("tma_id" %in% colnames(o@meta.data)) "tma_id" else "tma_id.Receiving"
    o$Condition <- as.character(o@meta.data[[gcol]])
    o$Core      <- paste0("core", as.character(o@meta.data[[tcol]]))
    objs[[basename(cd)]] <- o
    log(sprintf("  %s: %d columns, grade=%s", basename(cd), ncol(o), unique(o$Condition)[1]))
  }
  obj <- if (length(objs) == 1) objs[[1]] else
    merge(objs[[1]], y = objs[-1], add.cell.ids = names(objs))
  rm(objs); gc()

  # Seurat v5 keeps one layer per merged object (counts.1..counts.13); the vignettes predate
  # this. JoinLayers collapses them so ScaleData/FindVariableFeatures see one matrix.
  if (inherits(obj[[org]], "Assay5")) { log("JoinLayers (Seurat v5)"); obj <- JoinLayers(obj) }
  log(sprintf("merged: %d mechanisms x %d columns", nrow(obj), ncol(obj)))

  ## NOTES B4: low-information filter -- vignette 04:92, applied to the cell-cell object only
  if (org == "CellToCellSpatial") {
    nf <- paste0("nFeature_", org)
    before <- ncol(obj)
    obj <- subset(obj, cells = colnames(obj)[obj@meta.data[[nf]] > 5])
    log(sprintf("nFeature_%s > 5 filter: %d -> %d columns (%.1f%% kept)",
                org, before, ncol(obj), 100 * ncol(obj) / before))
  }

  # Checkpoint the merge BEFORE the embedding: the merge is the slow, memory-hungry part and
  # an embedding failure must not force it to be redone.
  saveRDS(obj, raw_rds); log("cached pre-embedding merge ->", raw_rds)

  ## NOTES B5-B10
  obj <- embed_merged(obj)
  saveRDS(obj, merged_rds)
  log("cached ->", merged_rds)
}

n_col <- ncol(obj)
npcs  <- ncol(Embeddings(obj, "pca"))
if (!is.null(obj@misc$hvf_selection_method)) HVF_METHOD <- obj@misc$hvf_selection_method

## ---------------- detection table: why FindAllMarkers does or doesn't fire ---------------
# FindAllMarkers' default min.pct = 0.25 needs a mechanism detected in >=25% of one group.
# Write the per-mechanism detection rate so any "no DE genes" result is backed by numbers
# rather than left as an unexplained blank.
{
  M <- LayerData(obj, assay = org, layer = "data")
  hi <- obj$Condition == "high"; lo <- obj$Condition == "low"
  det <- data.frame(
    mechanism      = rownames(M),
    frac_nonzero   = as.numeric(Matrix::rowMeans(M > 0)),
    frac_high      = as.numeric(Matrix::rowMeans(M[, hi, drop = FALSE] > 0)),
    frac_low       = as.numeric(Matrix::rowMeans(M[, lo, drop = FALSE] > 0)),
    mean_high      = as.numeric(Matrix::rowMeans(M[, hi, drop = FALSE])),
    mean_low       = as.numeric(Matrix::rowMeans(M[, lo, drop = FALSE])),
    stringsAsFactors = FALSE)
  det$max_group_pct <- pmax(det$frac_high, det$frac_low)
  det <- det[order(-det$max_group_pct), ]
  write.csv(det, file.path(out_dir, "mechanism_detection.csv"), row.names = FALSE)
  n_pass <- sum(det$max_group_pct >= 0.25)
  log(sprintf("detection: %d / %d mechanisms reach min.pct = 0.25 in at least one grade; "
              , n_pass, nrow(det)),
      sprintf("max detection rate = %.1f%%", 100 * max(det$max_group_pct)))
  rm(M); gc()
}

## ---------------- Stage B8-B11: diagnostics + embeddings --------------------------------
log("ElbowPlot / PCHeatmap / DimPlots")
save_plot(ElbowPlot(obj, ndims = npcs), file.path(plot_dir, "elbow"), 7, 5)
pc_hi <- min(48, npcs); pc_lo <- max(1, pc_hi - 8)
grDevices::png(file.path(plot_dir, "pc_heatmap.png"), 1400, 1600, res = 130)
PCHeatmap(obj, dims = pc_lo:pc_hi, cells = 100, balanced = TRUE); grDevices::dev.off()
grDevices::pdf(file.path(plot_dir, "pc_heatmap.pdf"), 11, 13)
PCHeatmap(obj, dims = pc_lo:pc_hi, cells = 100, balanced = TRUE); grDevices::dev.off()

group_vars <- if (org == "NeighborhoodToCell") {
  c("ReceivingType", "Condition", "Core")
} else {
  c("VectorType", "Condition", "Core", "SendingType", "ReceivingType")
}
for (g in group_vars) {
  if (!g %in% colnames(obj@meta.data)) next
  p <- DimPlot(obj, reduction = "umap", group.by = g, shuffle = TRUE, pt.size = 0.3) +
    ggtitle(sprintf("%s -- %s (%s)", org, g, tier))
  if (g == "VectorType") p <- p + NoLegend()
  save_plot(p, file.path(plot_dir, paste0("umap_", g)), 9, 7, npoints = n_col)
}

## composition: how much of each condition/celltype, written as a table not just a picture
comp_col <- if (org == "NeighborhoodToCell") "ReceivingType" else "VectorType"
comp <- as.data.frame(table(obj@meta.data[[comp_col]], obj$Condition))
colnames(comp) <- c(comp_col, "Condition", "n")
write.csv(comp, file.path(out_dir, sprintf("composition_%s_by_condition.csv", comp_col)),
          row.names = FALSE)

## ---------------- Stage C: differential signalling, high vs low grade -------------------
# C4 on the whole merged object first (vignette 01:108 finds markers per ident; here the
# ident of interest for the benchmark is grade).
run_markers <- function(o, tag) {
  f <- file.path(out_dir, "markers", sprintf("markers_%s.csv", tag))
  if (file.exists(f) && !force) { log("  cached markers", tag); return(read.csv(f)) }
  m <- tryCatch(
    FindAllMarkers(o, min.pct = 0.25, only.pos = TRUE, test.use = "roc", verbose = FALSE),
    error = function(e) { log("  FindAllMarkers failed for", tag, ":", conditionMessage(e)); NULL })
  if (is.null(m) || !nrow(m)) { write.csv(data.frame(), f, row.names = FALSE); return(NULL) }
  # NOTES C5: drop markers detected in only one group (infinite differential), vignette 07:130
  m$ratio <- m$pct.1 / m$pct.2
  m <- m[is.finite(m$ratio), ]
  write.csv(m, f, row.names = FALSE)
  m
}

# Vignette 01:108-111 -- markers of each cell type's NICHE, i.e. Idents = ReceivingType, not
# grade. This is a standard-workflow output in its own right AND the positive control for the
# grade contrast below: if this finds nothing either, the matrix is too sparse for
# FindAllMarkers' defaults rather than the grade contrast being null.
log("Stage C -- per-", loop_var0 <- if (org == "NeighborhoodToCell") "ReceivingType" else "VectorType",
    " niche markers (vignette 01:108)")
if (loop_var0 %in% colnames(obj@meta.data)) {
  Idents(obj) <- obj@meta.data[[loop_var0]]
  mk_ct <- run_markers(obj, paste0("ALL_by_", loop_var0))
  if (!is.null(mk_ct) && nrow(mk_ct)) {
    goi_ct <- mk_ct %>% group_by(cluster) %>% top_n(5, myAUC)
    p <- DoHeatmap(subset(obj, downsample = 300), features = unique(goi_ct$gene)) +
      scale_fill_gradientn(colors = c("grey", "white", "blue")) +
      ggtitle(sprintf("Niche mechanisms per %s -- %s (%s)", loop_var0, org, tier))
    save_plot(p, file.path(plot_dir, paste0("heatmap_ALL_by_", loop_var0)), 11, 14, npoints = 0)
  } else {
    log("  !! no markers by", loop_var0,
        "-- the matrix is too sparse for FindAllMarkers(min.pct = 0.25)")
  }
}

log("Stage C -- global grade contrast")
Idents(obj) <- obj$Condition
mk_all <- run_markers(obj, "ALL_by_grade")
if (!is.null(mk_all) && nrow(mk_all)) {
  goi <- mk_all %>% group_by(cluster) %>% top_n(20, myAUC)
  p <- DoHeatmap(subset(obj, downsample = 500), group.by = "ident",
                 features = unique(goi$gene)) +
    ggtitle(sprintf("Top DE mechanisms, high vs low grade -- %s (%s)", org, tier))
  save_plot(p, file.path(plot_dir, "heatmap_ALL_by_grade"), 11, 12, npoints = 0)
}

## C1: loop over every receiving population (vignette 07 hand-picks one; D6)
loop_var <- if (org == "NeighborhoodToCell") "ReceivingType" else "VectorType"
levels_all <- sort(unique(as.character(obj@meta.data[[loop_var]])))
if (org == "CellToCellSpatial") {
  # 81 possible VectorTypes is too many to embed one-by-one; take the most abundant that
  # actually have both conditions represented.
  tab <- table(obj@meta.data[[loop_var]], obj$Condition)
  keep <- rownames(tab)[apply(tab, 1, function(r) all(r >= 50))]
  levels_all <- head(keep[order(-rowSums(tab[keep, , drop = FALSE]))], 10)
  log("CellToCellSpatial: analysing top", length(levels_all), "VectorTypes with >=50 cells per grade")
}

summary_rows <- list()
for (lv in levels_all) {
  cells <- colnames(obj)[as.character(obj@meta.data[[loop_var]]) == lv]
  sub <- subset(obj, cells = cells)
  n_hi <- sum(sub$Condition == "high"); n_lo <- sum(sub$Condition == "low")
  safe <- gsub("[^A-Za-z0-9]+", "_", lv)
  log(sprintf("  %s = %s: %d columns (high %d / low %d)", loop_var, lv, ncol(sub), n_hi, n_lo))
  if (n_hi < 3 || n_lo < 3) { log("    too few in one grade -- skipping"); next }

  ## C2 embedding of the subset
  sub <- ScaleData(sub, verbose = FALSE)
  sub <- tryCatch(FindVariableFeatures(sub, selection.method = "disp", verbose = FALSE),
                  error = function(e) FindVariableFeatures(sub, selection.method = "vst",
                                                           verbose = FALSE))
  np <- min(50, nrow(sub) - 1, ncol(sub) - 1)
  sub <- RunPCA(sub, npcs = np, verbose = FALSE)
  sub <- RunUMAP(sub, dims = 1:min(40, np), verbose = FALSE)
  save_plot(DimPlot(sub, group.by = "Condition", shuffle = TRUE, pt.size = 0.4) +
              ggtitle(sprintf("%s -- %s", lv, tier)),
            file.path(plot_dir, paste0("umap_", safe, "_by_grade")), 7, 6, npoints = ncol(sub))

  ## C3-C7 differential
  Idents(sub) <- sub$Condition
  mk <- run_markers(sub, safe)
  if (!is.null(mk) && nrow(mk)) {
    goi <- mk %>% group_by(cluster) %>% top_n(20, myAUC)
    p <- DoHeatmap(subset(sub, downsample = 500), group.by = "ident",
                   features = unique(goi$gene)) +
      ggtitle(sprintf("Top DE mechanisms, high vs low grade: %s", lv))
    save_plot(p, file.path(plot_dir, paste0("heatmap_", safe)), 10, 10, npoints = 0)

    ## C8 (a) the method's OWN top-ranked LRs for this population
    for (gene in head(unique(goi$gene), 4)) {
      gs <- gsub("[^A-Za-z0-9]+", "_", gene)
      save_plot(FeaturePlot(sub, features = gene) + ggtitle(paste(lv, "--", gene)),
                file.path(plot_dir, "top_lr", paste0(safe, "__", gs)), 6.5, 5.5,
                npoints = ncol(sub))
    }
    summary_rows[[lv]] <- data.frame(
      level = lv, n_columns = ncol(sub), n_high = n_hi, n_low = n_lo,
      n_markers = nrow(mk),
      top_high = paste(head(mk$gene[mk$cluster == "high"], 5), collapse = "; "),
      top_low  = paste(head(mk$gene[mk$cluster == "low"],  5), collapse = "; "),
      stringsAsFactors = FALSE)
  }
  rm(sub); gc()
}
if (length(summary_rows))
  write.csv(bind_rows(summary_rows), file.path(out_dir, "differential_summary.csv"),
            row.names = FALSE)

## ---------------- C8 (b): the REQUESTED LRs, whatever their rank ------------------------
feat_present <- rownames(obj)
req_status <- list()
for (lrname in REQUESTED) {
  gs <- gsub("[^A-Za-z0-9]+", "_", lrname)
  if (!lrname %in% feat_present) {
    log("REQUESTED", lrname, "-- ABSENT from the NICHES output for this tier")
    req_status[[lrname]] <- list(present = FALSE, reason = "not in filtered ground truth")
    next
  }
  v <- FetchData(obj, vars = lrname)[, 1]
  nz <- sum(v > 0)
  log(sprintf("REQUESTED %s: nonzero in %d / %d columns (%.2f%%), mean %.4f",
              lrname, nz, length(v), 100 * nz / length(v), mean(v)))
  save_plot(FeaturePlot(obj, features = lrname, order = TRUE) +
              ggtitle(sprintf("%s -- %s (%s)", lrname, org, tier)),
            file.path(plot_dir, "requested_lr", paste0(gs, "_umap")), 7, 6, npoints = n_col)
  save_plot(VlnPlot(obj, features = lrname, group.by = "Condition", pt.size = 0) +
              ggtitle(sprintf("%s by grade", lrname)),
            file.path(plot_dir, "requested_lr", paste0(gs, "_vln_grade")), 6, 5, npoints = 0)
  gb <- if (org == "NeighborhoodToCell") "ReceivingType" else "ReceivingType"
  if (gb %in% colnames(obj@meta.data)) {
    save_plot(VlnPlot(obj, features = lrname, group.by = gb, split.by = "Condition",
                      pt.size = 0) + ggtitle(sprintf("%s by %s and grade", lrname, gb)),
              file.path(plot_dir, "requested_lr", paste0(gs, "_vln_", gb)), 10, 5, npoints = 0)
  }
  # per-grade quantification, so the picture is backed by numbers
  df <- data.frame(score = v, Condition = obj$Condition,
                   Group = as.character(obj@meta.data[[gb]]))
  agg <- df %>% group_by(Group, Condition) %>%
    summarise(n = n(), frac_nonzero = mean(score > 0), mean_score = mean(score),
              .groups = "drop")
  write.csv(agg, file.path(out_dir, sprintf("requested_%s_by_%s_grade.csv", gs, gb)),
            row.names = FALSE)
  req_status[[lrname]] <- list(present = TRUE, n_nonzero = nz, n_columns = length(v),
                               frac_nonzero = nz / length(v), mean_score = mean(v))
}

## ---------------- C12: spatial map of the requested LRs (D7: ggplot, no @images) ---------
if (org == "NeighborhoodToCell" && all(c("x", "y") %in% colnames(obj@meta.data))) {
  for (lrname in REQUESTED) {
    if (!lrname %in% feat_present) next
    gs <- gsub("[^A-Za-z0-9]+", "_", lrname)
    d <- data.frame(x = as.numeric(obj$x), y = as.numeric(obj$y),
                    score = FetchData(obj, vars = lrname)[, 1],
                    Core = obj$Core, Condition = obj$Condition)
    # Cores sit at different absolute positions on the slide. Centre each one so a single
    # fixed scale works -- facet_wrap(scales="free") and coord_fixed() are mutually exclusive
    # in ggplot2, and we want the aspect ratio preserved (these are real microns).
    d <- d %>% group_by(Core) %>% mutate(x = x - mean(x), y = y - mean(y)) %>% ungroup()
    d$panel <- paste0(d$Core, " (", d$Condition, ")")
    d <- d[order(d$score), ]
    p <- ggplot(d, aes(x = x, y = y, colour = score)) +
      geom_point(size = 0.25) +
      scale_colour_viridis_c(option = "magma") +
      facet_wrap(~ panel) +
      coord_fixed() + theme_cowplot(9) +
      labs(title = sprintf("%s -- NeighborhoodToCell score in situ (%s)", lrname, tier),
           x = "x (µm, centred)", y = "y (µm, centred)")
    save_plot(p, file.path(plot_dir, "requested_lr", paste0(gs, "_spatial")), 14, 11,
              npoints = n_col)
  }
}

## ---------------- manifest ---------------------------------------------------------------
write(jsonlite::toJSON(list(
  tier = tier, organization = org, n_mechanisms = nrow(obj), n_columns = n_col,
  n_pcs = npcs, loop_variable = loop_var, levels_analysed = levels_all,
  hvf_selection_method = HVF_METHOD,
  hvf_note = if (HVF_METHOD == "vst") {
    paste("vignette 07's selection.method='disp' overflows on this tier: NICHES stores raw",
          "LR products in the `data` slot and Seurat's disp method exponentiates them;",
          "max product exceeds exp()'s 709 domain. Fell back to Seurat's default 'vst'.")
  } else "vignette 07 value (selection.method = 'disp')",
  requested_lr = req_status,
  seurat_version = as.character(packageVersion("Seurat")),
  niches_version = as.character(packageVersion("NICHES"))
), auto_unbox = TRUE, pretty = TRUE, null = "null"),
file.path(out_dir, "analysis_manifest.json"))

log("DONE ->", out_dir)
