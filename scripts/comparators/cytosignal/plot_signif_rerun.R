#!/usr/bin/env Rscript
# Produce CytoSignal's stock `plotSignif` combination figures for a run that was executed in
# `nosave` mode (no cs_checkpoint.rds / cs_result.rds to replot from — the GBM full run).
#
# plotSignif needs @score + @imputation + @counts + @cells.loc, all of which `nosave` discards,
# so this script REBUILDS the object from the input dir using EXACTLY the parameters baked into
# run_cytosignal.R, then verifies the rebuild reproduced the original before plotting.
#
# The object is never written to disk (a 89k-cell cs object is ~10-15 GB; the box has ~16 GB free).
#
# Usage: source activate_env.sh
#        Rscript plot_signif_rerun.R <input_dir> <out_dir> [db_rds] [n_top]
#   input_dir : counts.mtx, genes.tsv, barcodes.tsv, meta.csv   (as run_cytosignal.R)
#   out_dir   : the EXISTING run dir — must contain quant/signif_summary_*.csv to verify against.
#               Plots are written to <out_dir>/plots/signif_<slot>/ ; nothing existing is overwritten.
#   db_rds    : custom DB from build_cellchat_db.R (omit -> bundled CellPhoneDB v2)
#   n_top     : interactions to plot per mode (default 6, matching the LUAD crop run)
#   extra     : comma-separated interaction IDs to ALSO plot, regardless of rank, into
#               plots/requested_<slot>/. Default = the ALARMIST motif-1 (mGAM) loop:
#               CCI-01109 (GRN -> SORT1) and CCI-01088 (ANXA1 -> FPR1).
suppressWarnings(suppressMessages({ library(Matrix); library(cytosignal) }))
args    <- commandArgs(trailingOnly = TRUE)
in_dir  <- if (length(args) >= 1) args[1] else stop("need input_dir")
out_dir <- if (length(args) >= 2) args[2] else stop("need out_dir")
db_rds  <- if (length(args) >= 3) args[3] else ""
n_top   <- if (length(args) >= 4) as.integer(args[4]) else 6L
extra_intr <- if (length(args) >= 5) trimws(strsplit(args[5], ",")[[1]]) else
                                     c("CCI-01109", "CCI-01088")   # GRN->SORT1, ANXA1->FPR1
# Slots to rank (SPARK) and plot. Omitted -> plot diffusion+contact and rank ALL slots, so the
# reproduction check can cover all three. Restricting slots makes a targeted second pass cheaper
# (SPARK is the dominant cost) at the price of a partial check.
plot_slots <- if (length(args) >= 6) trimws(strsplit(args[6], ",")[[1]]) else
                                     c("diffusion-Raw_smooth", "contact-Raw_smooth")
rank_slots <- if (length(args) >= 6) plot_slots else NULL
t0  <- Sys.time()
log <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n")
cdb <- if (nzchar(db_rds)) readRDS(db_rds) else NULL

## -------- load inputs (identical to run_cytosignal.R) --------
log("reading inputs from", in_dir)
dge <- as(as(Matrix::readMM(file.path(in_dir, "counts.mtx")), "CsparseMatrix"), "dgCMatrix")
genes <- readLines(file.path(in_dir, "genes.tsv")); cells <- readLines(file.path(in_dir, "barcodes.tsv"))
rownames(dge) <- genes; colnames(dge) <- cells
meta <- read.csv(file.path(in_dir, "meta.csv"), stringsAsFactors = FALSE)
stopifnot(all(meta$cell_id == cells))
loc <- as.matrix(meta[, c("x", "y")]); rownames(loc) <- meta$cell_id; colnames(loc) <- c("x", "y")
clust <- factor(meta$celltype); names(clust) <- meta$cell_id
log("dge:", nrow(dge), "genes x", ncol(dge), "cells; clusters:", nlevels(clust),
    "; DB:", if (is.null(cdb)) "bundled CellPhoneDBv2" else basename(db_rds))

## -------- rebuild: SAME parameters as run_cytosignal.R --------
cs <- createCytoSignal(raw.data = dge, cells.loc = loc, clusters = clust)
cs <- if (!is.null(cdb)) addIntrDB(cs, cdb$g_to_u, cdb$db.diff, cdb$db.cont, cdb$inter.index) else
                         addIntrDB(cs, g_to_u, db.diff, db.cont, inter.index)
cs <- removeLowQuality(cs, counts.thresh = 100, gene.thresh = 20)
cs <- changeUniprot(cs)
log("cells after QC:", ncol(cs@counts))
cs <- inferEpsParams(cs, scale.factor = 1, r.eps.real = 200)
log("findNN");   cs <- findNN(cs)
log("imputeLR"); cs <- imputeLR(cs)
log("inferIntrScore (perm.size=1e5, numCores=4)"); set.seed(42)
cs <- inferIntrScore(cs, perm.size = 1e5, numCores = 4)
log("inferSignif"); cs <- inferSignif(cs, p.thresh = 0.05, reads.thresh = 100, sig.thresh = 100)
tryCatch({ log("rankIntrSpatialVar (SPARK) on:", if (is.null(rank_slots)) "all slots" else paste(rank_slots, collapse = ", "))
           cs <- rankIntrSpatialVar(cs, slot.use = rank_slots, numCores = 4) },
         error = function(e) log("rankIntrSpatialVar FAILED:", conditionMessage(e)))

## -------- VERIFY the rebuild reproduces the original run --------
# The plots are only trustworthy if this object matches the one that produced quant/.
# inferIntrScore permutes under numCores=4, so an exact match is not guaranteed a priori —
# report the comparison rather than assume it.
pdir <- file.path(out_dir, "plots"); dir.create(pdir, showWarnings = FALSE, recursive = TRUE)
verify <- data.frame()
verify_slots <- if (is.null(rank_slots)) grep("_smooth$", names(cs@lrscore), value = TRUE) else rank_slots
for (slot in verify_slots) {
  tag <- gsub("[^A-Za-z0-9]+", "_", slot)
  ref_f <- file.path(out_dir, "quant", sprintf("signif_summary_%s.csv", tag))
  if (!file.exists(ref_f)) { log("verify: no reference for", slot, "- skipping"); next }
  ref <- read.csv(ref_f, stringsAsFactors = FALSE)
  cnt <- function(nm) { x <- cs@lrscore[[slot]]@res.list[[nm]]
                        if (is.null(x) || !length(x)) integer(0) else vapply(x, length, integer(1)) }
  nr <- cnt("result"); nh <- cnt("result.hq"); ns <- cnt("result.spx")
  g  <- function(v, ids) { o <- v[ids]; o[is.na(o)] <- 0L; as.integer(o) }
  new_hq  <- g(nh, ref$interaction_id); new_spx <- g(ns, ref$interaction_id)
  verify <- rbind(verify, data.frame(
    slot = slot, n_ref_intr = nrow(ref),
    n_new_intr = length(unique(unlist(lapply(cs@lrscore[[slot]]@res.list, names)))),
    hq_exact   = sum(new_hq == ref$n_hq),  hq_cor  = suppressWarnings(cor(new_hq,  ref$n_hq)),
    spx_exact  = sum(new_spx == ref$n_spx), spx_cor = suppressWarnings(cor(new_spx, ref$n_spx)),
    stringsAsFactors = FALSE))
}
cat("\n== REPRODUCTION CHECK vs existing quant/ ==\n"); print(verify, row.names = FALSE)
# Never clobber a passing check from an earlier, more complete pass.
vf <- file.path(pdir, "reproduction_check.csv")
if (file.exists(vf)) vf <- file.path(pdir, "reproduction_check_pass2.csv")
write.csv(verify, vf, row.names = FALSE)

## -------- plotSignif: the stock per-interaction combination figure --------
# Raw-Raw_smooth is deliberately skipped: it is the multi-cell-spot readout (Visium-like) and is
# not the intended slot for single-cell-resolution Xenium. Matches the LUAD crop run.
for (slot in plot_slots) {
  tag <- gsub("[^A-Za-z0-9]+", "_", slot)
  sig_ids <- tryCatch(showIntr(cs, slot.use = slot, signif.use = "result.spx", return.name = FALSE),
                      error = function(e) { log("showIntr err", slot, ":", conditionMessage(e)); character(0) })
  n <- length(sig_ids)

  # (a) the method's OWN most significant interactions: top-N by SPARK-X spatial variability
  ntop <- min(n_top, n); log(slot, "-> plotting top", ntop, "of", n, "spatially-variable interactions")
  if (ntop > 0) tryCatch(
    plotSignif(cs, intr = seq_len(ntop), slot.use = slot, signif.use = "result.spx",
               plot_dir = file.path(pdir, paste0("signif_", tag, "/")),
               pt.size = 0.25, plot.fmt = "png", raster = TRUE),
    error = function(e) log("plotSignif FAILED", slot, ":", conditionMessage(e)))

  # (b) explicitly requested interactions (ALARMIST motif-1 mGAM loop: GRN->SORT1, ANXA1->FPR1).
  #     Plotted into a separate dir so they are never confused with the method's own ranking.
  want <- intersect(extra_intr, sig_ids)
  miss <- setdiff(extra_intr, sig_ids)
  if (length(miss)) log(slot, "-> requested but NOT significant here, skipped:", paste(miss, collapse = ", "))
  if (length(want)) {
    log(slot, "-> plotting requested:", paste(want, collapse = ", "))
    tryCatch(
      plotSignif(cs, intr = want, slot.use = slot, signif.use = "result.spx",
                 plot_dir = file.path(pdir, paste0("requested_", tag, "/")),
                 pt.size = 0.25, plot.fmt = "png", raster = TRUE),
      error = function(e) log("plotSignif FAILED (requested)", slot, ":", conditionMessage(e)))
  }
}

## -------- COMPLETENESS PASS: every other applicable CytoSignal plot --------
## First pass produced only plotCluster + plotSignif. These are the rest.
## Produced in THIS session so the (very large) object need not be persisted first.
for (slot in plot_slots) {
  tag <- gsub("[^A-Za-z0-9]+", "_", slot)
  sig_ids <- tryCatch(showIntr(cs, slot.use = slot, signif.use = "result.spx", return.name = FALSE),
                      error = function(e) character(0))
  if (!length(sig_ids)) next
  pick <- c(head(sig_ids, 3), intersect(extra_intr, sig_ids))
  ## plotEdge: the 3D sender -> receiver edge view the tutorial demonstrates
  for (ii in pick) {
    for (ty in c("sender", "receiver")) {
      tryCatch({
        plotEdge(cs, intr = ii, type = ty, slot.use = slot, signif.use = "result.spx",
                 plot_dir = file.path(pdir, paste0("edge_", tag, "/")),
                 filename = paste0(ii, "_", ty, ".png"), plot.fmt = "png",
                 return.plot = FALSE, pt.size = 0.05, edge.size = 300)
        log("  plotEdge ok:", ii, ty)
      }, error = function(e) log("  plotEdge FAILED", ii, ty, ":", conditionMessage(e)))
    }
  }
  ## plotIntrValue: ligand / receptor / score / null panels for one interaction
  for (ii in head(pick, 2)) {
    tryCatch({
      png(file.path(pdir, sprintf("intrValue_%s_%s.png", tag, ii)), 1600, 1200, res = 150)
      print(plotIntrValue(cs, intr = ii, slot.use = slot, signif.use = "result.spx",
                          pt.size = 0.15, raster = TRUE)); dev.off()
      log("  plotIntrValue ok:", ii)
    }, error = function(e) { try(dev.off(), silent = TRUE)
      log("  plotIntrValue FAILED", ii, ":", conditionMessage(e)) })
  }
  ## plotCircosNIntr: circos of interaction COUNTS between cell types -- the
  ## "which cell type talks to which" summary that was entirely missing
  tryCatch({
    png(file.path(pdir, sprintf("circos_nIntr_%s.png", tag)), 1600, 1600, res = 170)
    print(plotCircosNIntr(cs, lrscore.use = slot, intr.use = head(sig_ids, 30))); dev.off()
    log("  plotCircosNIntr ok:", slot)
  }, error = function(e) { try(dev.off(), silent = TRUE)
    log("  plotCircosNIntr FAILED", slot, ":", conditionMessage(e)) })
  ## plotSigCluster: top interactions per cluster, sender and receiver
  for (ty in c("sender", "receiver")) {
    tryCatch({
      plotSigCluster(cs, plot_dir = file.path(pdir, paste0("sigCluster_", tag, "_", ty, "/")),
                     intr.num = 5, type = ty, slot.use = slot, signif.use = "result.spx",
                     plot.fmt = "png", use.cex = 0.05, edge.size = 500)
      log("  plotSigCluster ok:", slot, ty)
    }, error = function(e) log("  plotSigCluster FAILED", slot, ty, ":", conditionMessage(e)))
  }
}

## -------- PERSIST the object if disk allows (skill: results must be re-readable) --------
free_gb <- tryCatch(as.numeric(system(sprintf("df -g %s | tail -1 | awk '{print $4}'", pdir),
                                      intern = TRUE)), error = function(e) NA)
log("free disk:", free_gb, "GB")
if (!is.na(free_gb) && free_gb >= 20) {
  tryCatch({ saveRDS(purgeBeforeSave(cs), file.path(out_dir, "cs_result.rds"))
             log("saved cs_result.rds") },
           error = function(e) log("saveRDS FAILED:", conditionMessage(e)))
} else {
  log("NOT saving cs_result.rds: needs ~10-15 GB for 89k cells and only", free_gb,
      "GB free. Recorded in DEVIATIONS.md as a disk limitation, not an oversight.")
}

man <- c(script = "plot_signif_rerun.R", input_dir = in_dir, out_dir = out_dir,
         db = if (nzchar(db_rds)) basename(db_rds) else "bundled CellPhoneDBv2",
         n_cells_in = ncol(dge), n_cells_qc = ncol(cs@counts), n_top = n_top, seed = 42,
         params = "counts.thresh=100 gene.thresh=20 scale.factor=1 r.eps.real=200 perm.size=1e5",
         cytosignal = as.character(packageVersion("cytosignal")),
         R = R.version.string,
         wall_min = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1))
writeLines(c("{", paste0(sprintf('  "%s": "%s"', names(man), man), collapse = ",\n"), "}"),
           file.path(out_dir, "plots", "plot_signif_manifest.json"))
log("DONE in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min. plots in", pdir)
