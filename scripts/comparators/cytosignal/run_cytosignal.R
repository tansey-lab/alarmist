#!/usr/bin/env Rscript
# CytoSignal run on a CytoSignal-format input dir (counts.mtx + meta.csv).
# Usage: Rscript run_cytosignal.R <input_dir> <output_dir> [save|nosave] [custom_db.rds]
#   input_dir     : counts.mtx (genes x cells, integer), genes.tsv, barcodes.tsv, meta.csv(cell_id,x,y,celltype)
#   save|nosave   : whether to save the multi-GB cs objects (nosave = disk-safe for large runs)
#   custom_db.rds : optional list(g_to_u, db.diff, db.cont, inter.index) from build_cellchat_db.R
# Env: comp-cytosignal. Xenium coords in microns -> scale.factor = 1.
# Persists QUANTITATIVE outputs (quant/): per-slot score matrix + full @res.list + significance summary.
suppressWarnings(suppressMessages({ library(Matrix); library(cytosignal) }))
args <- commandArgs(trailingOnly = TRUE)
script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(script_dir, "quant_io.R"))

in_dir   <- if (length(args) >= 1) args[1] else stop("need input_dir")
out_dir  <- if (length(args) >= 2) args[2] else file.path(in_dir, "..", "run")
save_rds <- if (length(args) >= 3) tolower(args[3]) else "save"
db_rds   <- if (length(args) >= 4) args[4] else ""
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
ckpt <- file.path(out_dir, "cs_checkpoint.rds")
log  <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n")
cdb  <- if (nzchar(db_rds)) readRDS(db_rds) else NULL

## -------- load inputs --------
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
print(summary(Matrix::colSums(dge)))

## -------- pipeline --------
if (file.exists(ckpt) && save_rds == "save") {
  log("resuming from checkpoint", ckpt); cs <- readRDS(ckpt)
} else {
  cs <- createCytoSignal(raw.data = dge, cells.loc = loc, clusters = clust)
  if (!is.null(cdb)) cs <- addIntrDB(cs, cdb$g_to_u, cdb$db.diff, cdb$db.cont, cdb$inter.index)
  else               cs <- addIntrDB(cs, g_to_u, db.diff, db.cont, inter.index)
  log("removeLowQuality (counts>=100, genes>=20) — tuned for Xenium panel")
  cs <- removeLowQuality(cs, counts.thresh = 100, gene.thresh = 20)
  cs <- changeUniprot(cs)
  log("cells after QC:", ncol(cs@counts))
  cs <- inferEpsParams(cs, scale.factor = 1, r.eps.real = 200)     # Xenium coords already in microns
  log("findNN"); cs <- findNN(cs)
  log("imputeLR"); cs <- imputeLR(cs)
  log("inferIntrScore (perm.size=1e5, numCores=4)"); set.seed(42)
  cs <- inferIntrScore(cs, perm.size = 1e5, numCores = 4)
  if (save_rds == "save") { saveRDS(cs, ckpt); log("checkpoint saved") } else log("disk-safe: skip checkpoint")
}

log("lrscore slots:", paste(names(cs@lrscore), collapse = ", "))
log("inferSignif"); cs <- inferSignif(cs, p.thresh = 0.05, reads.thresh = 100, sig.thresh = 100)
## -------- interaction names (HOISTED above the SPARK step: save_lrscore_quant needs them,
##          and the whole point below is to call it before SPARK can die) --------
used_ii <- if (!is.null(cdb)) cdb$inter.index else inter.index
intr_names <- setNames(paste(sub("-ligand$", "", used_ii$partner_a),
                             sub("-receptor$", "", used_ii$partner_b), sep = " - "),
                       used_ii$id_cp_interaction)

## -------- BANK the quantitative outputs BEFORE rankIntrSpatialVar (added 2026-08-18) --------
## On 2026-08-18 P17_LUAD was OOM-killed inside rankIntrSpatialVar after inferSignif had already
## run to completion -- 6 h 42 min of finished work, four statements from disk, all discarded.
## MEASURED, not inferred: sstat -j 8869477.0 reports MaxRSS 100,662,144 K = 95.9989 GiB against
## a --mem=96G cgroup cap, i.e. 1.12 MiB of headroom left, so the kernel cgroup OOM killer
## SIGKILLed the R process. The tryCatch below CANNOT help with that: it catches R conditions,
## and a SIGKILL is not one. Neither can the checkpoint at :52 -- that is written BEFORE
## inferSignif, so even `save` mode would have protected only the preceding ~53 min.
## Ordering was the real defect. This call writes everything except the n_spx column;
## save_lrscore_quant tolerates a missing result.spx (quant_io.R:45-47 -> cnt() returns
## integer(0) -> g() fills 0L), and the summary CSV is sorted by -n_hq, not n_spx, so the
## pre-SPARK file is already correctly ordered. The post-SPARK call below overwrites it with
## n_spx filled in. Cost: one extra write of quant/ per section.
log("banking quant/ BEFORE the SPARK step (n_spx will be 0 until the rewrite below)")
save_lrscore_quant(cs, file.path(out_dir, "quant"), intr_names = intr_names, log = log)

## -------- SPARK-X ranking -- SKIPPABLE via CS_SKIP_SPX=1 --------
## What the step actually does, VERIFIED on the completed P17_AIS outputs rather than assumed:
## result.spx is a RE-ORDERING of result.hq, not a filter. Across all three slots it holds the
## same interactions (146/146, 642/642, 68/68) with byte-identical per-interaction cell sets,
## only in a different order, and the package discards the SPARK-X adjusted p-values. So
## skipping it costs the ordering, the n_spx column and the result.spx element -- no
## significance call, no score matrix, no cell set, and not the NEBULA differential, which
## rebuilds from input/ and never reads these outputs. It IS still a deviation from the authors'
## documented step order; declare it if you use it. Default is OFF: the step runs.
if (nzchar(Sys.getenv("CS_SKIP_SPX"))) {
  log("rankIntrSpatialVar SKIPPED (CS_SKIP_SPX set); n_spx stays 0, result.spx absent")
} else {
  tryCatch({ log("rankIntrSpatialVar (SPARK)"); cs <- rankIntrSpatialVar(cs, numCores = 4) },
           error = function(e) log("rankIntrSpatialVar FAILED:", conditionMessage(e)))
}
if (save_rds == "save") saveRDS(purgeBeforeSave(cs), file.path(out_dir, "cs_result.rds")) else log("disk-safe: skip cs_result.rds")

## -------- rewrite quant/ with n_spx filled in (both modes) --------
save_lrscore_quant(cs, file.path(out_dir, "quant"), intr_names = intr_names, log = log)

## -------- plots (only when small enough; heavy at >200k cells) --------
if (ncol(cs@counts) <= 200000) {
  tryCatch({ png(file.path(out_dir, "cluster_map.png"), 1400, 1200, res = 150); print(plotCluster(cs)); dev.off() },
           error = function(e) log("plotCluster err:", conditionMessage(e)))
  ds <- grep("diffusion.*_smooth", names(cs@lrscore), value = TRUE)[1]
  nsig <- tryCatch(length(showIntr(cs, ds, "result.spx", return.name = FALSE)), error = function(e) 0)
  if (length(nsig) && !is.na(nsig) && nsig > 0) tryCatch(
    plotSignif(cs, intr = seq_len(min(6, nsig)), slot.use = ds, signif.use = "result.spx",
               plot_dir = file.path(out_dir, "signif_plots/"), pt.size = 0.3, plot.fmt = "png", raster = TRUE),
    error = function(e) log("plotSignif err:", conditionMessage(e)))
}
log("DONE. outputs in", out_dir)
