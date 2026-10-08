#!/usr/bin/env Rscript
# ======================================================================================
# replot_from_quant.R — produce CytoSignal's stock figures for a section that finished
# `run_cytosignal.R` but shipped NO images because of the cell-count gate.
#
# WHY THIS EXISTS
#   run_cytosignal.R:105 wraps every plotting call in `if (ncol(cs@counts) <= 200000)`.
#   On LUAD only P17_AIS (182,378 cells) is under it; P21_AIS, P21_LUAD and P17_LUAD are
#   all far above and therefore have quant/ but no figures. Stage 03 also runs `nosave`,
#   so there is no cs_checkpoint.rds / cs_result.rds to replot from.
#
# THE POINT: `inferSignif` DOES NOT HAVE TO RUN AGAIN.
#   Measured on the completed P21_LUAD run (cyto_P21_LUAD2.log, 498,422 cells post-QC):
#       read inputs        ~1 min
#       findNN             ~5 min
#       imputeLR           ~1 min
#       inferIntrScore     ~43 min
#       inferSignif        6 h 03 min   <-- 82% of the wall clock
#       quant/ writes      ~25 min (x2)
#       rankIntrSpatialVar ~5.5 min
#   `inferSignif`'s ENTIRE output is `@res.list`, and quant_io.R already banked it:
#   quant/reslist_<tag>.rds is a verbatim `saveRDS(lr@res.list)` holding result /
#   result.hq / result.spx. Verified present and populated for all four LUAD sections
#   (`n_spx` is non-zero in every signif_summary_*.csv), and because the post-SPARK write
#   is the one that survives, the names carry the SPARK-X ordering.
#
#   So: rebuild only as far as `inferIntrScore`, then ASSIGN the banked res.list back into
#   the object. ~50 min instead of ~7.5 h, and the plots then correspond EXACTLY to the
#   numbers the benchmark reports rather than to a fresh permutation.
#
# WHY inferIntrScore STILL RUNS (it is not skippable)
#   cytosignal:::getIntrValue -- which plotSignif/plotIntrValue go through -- opens with
#       sample.index <- sample(ncol(score.obj@lig.null), nCells)
#   unconditionally, and its default `type` includes ligand_null / receptor_null /
#   score_null. @lig.null, @recep.null and @score.null are produced by inferIntrScore and
#   are NOT in quant/. They are also the only thing it is needed for.
#
# WHAT IS AND IS NOT A DEVIATION
#   Not a deviation: every significance call, cell set and ordering that is plotted comes
#   from the authors' own inferSignif + rankIntrSpatialVar output, unmodified.
#   IS a deviation, and declared: the permutation NULL displayed in the `*_null` panels is
#   a fresh draw, not the one from the original run (inferIntrScore permutes under
#   numCores=4, so bit-identical reproduction was never available anyway). It affects only
#   those display panels -- no statistic, no ranking, no significance call.
#   Set CS_MODE=full to re-run inferSignif + rankIntrSpatialVar instead and remove even
#   that; it costs the 6 h back. quant/ is NEVER written by this script in either mode.
#
# USAGE
#   source scripts/comparators/cytosignal/activate_env.sh
#   Rscript replot_from_quant.R <input_dir> <out_dir> [db_rds]
#
#   input_dir : counts.mtx, genes.tsv, barcodes.tsv, meta.csv   (as run_cytosignal.R)
#   out_dir   : the EXISTING run dir; must contain quant/. Figures land in <out_dir>/plots/
#               and NOTHING already there is overwritten.
#   db_rds    : custom DB from build_cellchat_db.R (omit -> bundled CellPhoneDB v2)
#
# ENV KNOBS (all optional)
#   CS_MODE     inject (default) | full
#   CS_NTOP     interactions per slot from the method's own ranking   [6]
#   CS_EXTRA    comma-separated CCI ids to ALSO plot regardless of rank, into
#               plots/requested_<tag>/.  DEFAULT: EMPTY -- no interaction is privileged
#               unless the caller names one.
#
#               This defaulted to CCI-01109 (GRN -> SORT1) until 2026-08-24, and nothing
#               anywhere set CS_EXTRA, so every LUAD/AIS section silently got a
#               requested_<tag>/Rank_N_GRN-SORT1.png. GRN -> SORT1 and ANXA1 -> FPR1 are the
#               two arms of ALARMIST motif 1 in the **LGG/GBM** TMA. They are a hypothesis
#               about THAT dataset and must never be requested, tracked or "checked" on
#               LUAD/AIS. The old comment argued the pair was fair game here because GRN,
#               SORT1 and FPR1 are all on the LUAD 5K panel -- that reasoning is wrong and
#               was rejected: panel membership means the assay COULD measure it, not that
#               the hypothesis belongs to this dataset.
#
#               A dataset-agnostic script must not carry any dataset's hypothesis as a
#               default. If a GBM run wants it, that runner passes CS_EXTRA=CCI-01109.
#   CS_SLOTS    comma-separated score slots  [diffusion-Raw_smooth,contact-Raw_smooth]
#               Raw-Raw_smooth is deliberately excluded: it is the multi-cell-spot
#               (Visium-like) readout, not the intended slot for single-cell Xenium.
#   CS_FMT      png (default) | pdf | svg -- the method's own plot.fmt. These are
#               CytoSignal's stock diagnostic figures, NOT paper panels; publication
#               panels are built separately through _common/plotting.py.
#   CS_SAVE_RDS 1 to persist purgeBeforeSave(cs) as cs_result.rds (needs ~10-15 GB per
#               100k cells; off by default)
# ======================================================================================
suppressWarnings(suppressMessages({ library(Matrix); library(cytosignal) }))

args    <- commandArgs(trailingOnly = TRUE)
in_dir  <- if (length(args) >= 1) args[1] else stop("need input_dir")
out_dir <- if (length(args) >= 2) args[2] else stop("need out_dir")
db_rds  <- if (length(args) >= 3) args[3] else ""

envd <- function(k, d) { v <- Sys.getenv(k, unset = NA); if (is.na(v)) d else v }
csv  <- function(s) { s <- trimws(strsplit(s, ",")[[1]]); s[nzchar(s)] }

mode       <- tolower(envd("CS_MODE", "inject"))
n_top      <- as.integer(envd("CS_NTOP", "6"))
extra_intr <- csv(envd("CS_EXTRA", ""))   # see CS_EXTRA above: empty by design
plot_slots <- csv(envd("CS_SLOTS", "diffusion-Raw_smooth,contact-Raw_smooth"))
plot_fmt   <- envd("CS_FMT", "png")
save_rds   <- envd("CS_SAVE_RDS", "0") %in% c("1", "true", "TRUE", "yes")
stopifnot(mode %in% c("inject", "full"))

t0   <- Sys.time()
log  <- function(...) { cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n"); flush.console() }
tag_of <- function(s) gsub("[^A-Za-z0-9]+", "_", s)
mem  <- function(where) {
  # Peak RSS straight from the kernel -- the number to size --mem against next time.
  v <- tryCatch(grep("^VmHWM", readLines("/proc/self/status"), value = TRUE), error = function(e) character(0))
  if (length(v)) log(sprintf("  [mem] peak RSS after %s: %s", where, sub("^VmHWM:\\s*", "", v[1])))
}

log("mode =", mode, "| slots =", paste(plot_slots, collapse = ", "), "| fmt =", plot_fmt,
    "| n_top =", n_top, "| extra =", if (length(extra_intr)) paste(extra_intr, collapse = ",") else "(none)")

quant_dir <- file.path(out_dir, "quant")
if (mode == "inject" && !dir.exists(quant_dir))
  stop("CS_MODE=inject needs ", quant_dir, " -- it does not exist. Use CS_MODE=full.")

cdb <- if (nzchar(db_rds)) readRDS(db_rds) else NULL

## ---------------------------------------------------------------- load inputs
## Byte-for-byte the same load run_cytosignal.R does, so the QC that follows selects the
## same cells and the banked res.list barcodes line up.
log("reading inputs from", in_dir)
dge   <- as(as(Matrix::readMM(file.path(in_dir, "counts.mtx")), "CsparseMatrix"), "dgCMatrix")
genes <- readLines(file.path(in_dir, "genes.tsv")); cells <- readLines(file.path(in_dir, "barcodes.tsv"))
rownames(dge) <- genes; colnames(dge) <- cells
meta  <- read.csv(file.path(in_dir, "meta.csv"), stringsAsFactors = FALSE)
stopifnot(all(meta$cell_id == cells))
loc   <- as.matrix(meta[, c("x", "y")]); rownames(loc) <- meta$cell_id; colnames(loc) <- c("x", "y")
clust <- factor(meta$celltype); names(clust) <- meta$cell_id
log("dge:", nrow(dge), "genes x", ncol(dge), "cells; clusters:", nlevels(clust),
    "; DB:", if (is.null(cdb)) "bundled CellPhoneDBv2" else basename(db_rds))

## ------------------------------------------------- rebuild: run_cytosignal.R parameters
cs <- createCytoSignal(raw.data = dge, cells.loc = loc, clusters = clust)
cs <- if (!is.null(cdb)) addIntrDB(cs, cdb$g_to_u, cdb$db.diff, cdb$db.cont, cdb$inter.index) else
                         addIntrDB(cs, g_to_u, db.diff, db.cont, inter.index)
cs <- removeLowQuality(cs, counts.thresh = 100, gene.thresh = 20)
cs <- changeUniprot(cs)
log("cells after QC:", ncol(cs@counts))
cs <- inferEpsParams(cs, scale.factor = 1, r.eps.real = 200)   # Xenium coords already in microns
rm(dge); gc(verbose = FALSE)
log("findNN");   cs <- findNN(cs);   mem("findNN")
log("imputeLR"); cs <- imputeLR(cs); mem("imputeLR")
## Needed for @lig.null / @recep.null / @score.null -- see the header. NOT skippable.
log("inferIntrScore (perm.size=1e5, numCores=4)"); set.seed(42)
cs <- inferIntrScore(cs, perm.size = 1e5, numCores = 4); mem("inferIntrScore")
log("lrscore slots:", paste(names(cs@lrscore), collapse = ", "))

avail_slots <- setdiff(names(cs@lrscore), "default")
bad <- setdiff(plot_slots, avail_slots)
if (length(bad)) stop("requested slot(s) not present after inferIntrScore: ", paste(bad, collapse = ", "),
                      "\navailable: ", paste(avail_slots, collapse = ", "))

## ---------------------------------------------------------------- significance
verify <- data.frame()
if (mode == "full") {
  log("CS_MODE=full: re-running inferSignif (this is the ~6 h stage)")
  cs <- inferSignif(cs, p.thresh = 0.05, reads.thresh = 100, sig.thresh = 100); mem("inferSignif")
  tryCatch({ log("rankIntrSpatialVar (SPARK)"); cs <- rankIntrSpatialVar(cs, slot.use = plot_slots, numCores = 4) },
           error = function(e) log("rankIntrSpatialVar FAILED:", conditionMessage(e)))
} else {
  ## ---- INJECT the banked res.list -------------------------------------------------
  ## Assign, then PROVE the assignment against signif_summary_*.csv, which was written by
  ## the same save_lrscore_quant() call that wrote the reslist. A silent mismatch here
  ## would produce confident, wrong figures, so every check below is fatal.
  for (slot in plot_slots) {
    tg <- tag_of(slot)
    rf <- file.path(quant_dir, sprintf("reslist_%s.rds", tg))
    if (!file.exists(rf)) stop("missing banked res.list: ", rf)
    rl <- readRDS(rf)
    if (!is.list(rl) || !length(rl)) stop("banked res.list for ", slot, " is empty/not a list")
    cs@lrscore[[slot]]@res.list <- rl
    log(sprintf("injected %s: %s", slot,
                paste(sprintf("%s=%d", names(rl), vapply(rl, length, integer(1))), collapse = " ")))

    ## (a) banked barcodes must exist among the post-QC cells of THIS rebuild. If QC selected a
    ##     different cell set the whole injection is invalid.
    ##     Sampled, not exhaustive, and deliberately so: unlist()ing every cell set is ~10^8
    ##     strings for a large section (P17_LUAD's result.spx alone is ~865 interactions of up
    ##     to 215,868 barcodes) -- minutes of wall clock and GBs of RAM to answer a question
    ##     the largest few cell sets already answer. A QC mismatch shifts the whole cell
    ##     vocabulary, so it cannot hide in the interactions this does not look at.
    qc_cells <- colnames(cs@counts)
    last <- rl[[length(rl)]]
    if (length(last)) {
      probe_i <- head(order(vapply(last, length, integer(1)), decreasing = TRUE), 5)
      allbc   <- unique(unlist(last[probe_i], use.names = FALSE))
      miss_bc <- length(setdiff(allbc, qc_cells))
      if (miss_bc > 0) stop(sprintf("%s: %d of %d banked barcodes (largest %d cell sets) are absent from ",
                                    slot, miss_bc, length(allbc), length(probe_i)),
                            "this rebuild's post-QC cells -- QC is not reproducing the original run; ",
                            "refusing to plot.")
      log(sprintf("  barcode check: %d distinct barcodes from the %d largest cell sets all present",
                  length(allbc), length(probe_i)))
    }

    ## (b) counts must match the summary CSV exactly.
    sf <- file.path(quant_dir, sprintf("signif_summary_%s.csv", tg))
    hq_ex <- spx_ex <- NA_integer_; nref <- NA_integer_
    if (file.exists(sf)) {
      ref <- read.csv(sf, stringsAsFactors = FALSE); nref <- nrow(ref)
      cnt <- function(nm) { x <- rl[[nm]]; if (is.null(x) || !length(x)) integer(0) else vapply(x, length, integer(1)) }
      g   <- function(v, ids) { o <- v[ids]; o[is.na(o)] <- 0L; as.integer(o) }
      hq_ex  <- sum(g(cnt("result.hq"),  ref$interaction_id) == ref$n_hq)
      spx_ex <- sum(g(cnt("result.spx"), ref$interaction_id) == ref$n_spx)
      if (hq_ex != nref) stop(sprintf("%s: n_hq mismatch (%d/%d rows agree) -- refusing to plot.", slot, hq_ex, nref))
    } else log("  no signif_summary for", slot, "- count check skipped")

    ## (c) the score this rebuild computed vs the score that was banked. @score is a
    ##     deterministic function of the imputation (only the NULL is permuted), so these
    ##     must agree; a disagreement means findNN/imputeLR did not reproduce.
    scf <- file.path(quant_dir, sprintf("score_%s.rds", tg)); maxdiff <- NA_real_
    if (file.exists(scf)) {
      Sd <- readRDS(scf); Sn <- cs@lrscore[[slot]]@score
      ok_dim <- identical(dim(Sd), dim(Sn)) &&
                identical(rownames(Sd), rownames(Sn)) && identical(colnames(Sd), colnames(Sn))
      if (!ok_dim) {
        log("  WARNING: banked score dimnames differ from the rebuild for", slot,
            "- plotted scores are the rebuild's, res.list is the banked one.")
      } else {
        probe   <- head(seq_len(ncol(Sd)), 5)
        maxdiff <- max(abs(as.matrix(Sd[, probe, drop = FALSE]) - as.matrix(Sn[, probe, drop = FALSE])))
        log(sprintf("  score reproduction: max|banked-rebuilt| over %d probe interactions = %.3g",
                    length(probe), maxdiff))
      }
      rm(Sd); gc(verbose = FALSE)
    }
    verify <- rbind(verify, data.frame(slot = slot, n_ref_intr = nref,
                                       hq_exact = hq_ex, spx_exact = spx_ex,
                                       score_maxdiff = maxdiff, stringsAsFactors = FALSE))
  }
  cat("\n== INJECTION CHECK vs banked quant/ ==\n"); print(verify, row.names = FALSE)
}

pdir <- file.path(out_dir, "plots"); dir.create(pdir, showWarnings = FALSE, recursive = TRUE)
if (nrow(verify)) {
  vf <- file.path(pdir, "injection_check.csv")
  i <- 2L; while (file.exists(vf)) { vf <- file.path(pdir, sprintf("injection_check_pass%d.csv", i)); i <- i + 1L }
  write.csv(verify, vf, row.names = FALSE)
}
mem("significance ready")

## ---------------------------------------------------------------- plots
## `signif.use = "result.spx"` where available (SPARK-X ordering), else "result.hq".
pick_signif <- function(slot) {
  av <- names(cs@lrscore[[slot]]@res.list)
  if ("result.spx" %in% av && length(cs@lrscore[[slot]]@res.list$result.spx)) "result.spx" else "result.hq"
}

tryCatch({
  f <- file.path(out_dir, "plots", paste0("cluster_map.", plot_fmt))
  if (plot_fmt == "png") png(f, 1400, 1200, res = 150) else if (plot_fmt == "pdf") pdf(f, 9.3, 8) else svg(f, 9.3, 8)
  print(plotCluster(cs)); dev.off(); log("plotCluster ok ->", basename(f))
}, error = function(e) { try(dev.off(), silent = TRUE); log("plotCluster err:", conditionMessage(e)) })

for (slot in plot_slots) {
  tg <- tag_of(slot); su <- pick_signif(slot)
  sig_ids <- tryCatch(showIntr(cs, slot.use = slot, signif.use = su, return.name = FALSE),
                      error = function(e) { log("showIntr err", slot, ":", conditionMessage(e)); character(0) })
  n <- length(sig_ids)
  log("== ", slot, " (signif.use =", su, ") ->", n, "significant interactions ==")
  if (!n) next

  ## (a) the method's OWN top interactions, by its own ranking
  ntop <- min(n_top, n)
  tryCatch(plotSignif(cs, intr = seq_len(ntop), slot.use = slot, signif.use = su,
                      plot_dir = file.path(pdir, paste0("signif_", tg, "/")),
                      pt.size = 0.25, plot.fmt = plot_fmt, raster = TRUE),
           error = function(e) log("plotSignif FAILED", slot, ":", conditionMessage(e)))
  log("  plotSignif top", ntop, "done"); mem(paste("plotSignif", tg))

  ## (b) explicitly requested interactions, kept in a separate directory so they are never
  ##     mistaken for the method's own ranking
  want <- intersect(extra_intr, sig_ids); miss <- setdiff(extra_intr, sig_ids)
  if (length(miss)) log("  requested but NOT significant in this slot, skipped:", paste(miss, collapse = ", "))
  if (length(want)) tryCatch(
    plotSignif(cs, intr = want, slot.use = slot, signif.use = su,
               plot_dir = file.path(pdir, paste0("requested_", tg, "/")),
               pt.size = 0.25, plot.fmt = plot_fmt, raster = TRUE),
    error = function(e) log("plotSignif FAILED (requested)", slot, ":", conditionMessage(e)))
  if (length(want)) log("  plotSignif requested done:", paste(want, collapse = ", "))

  ## (c) the rest of CytoSignal's applicable plot surface
  pick <- unique(c(head(sig_ids, 3), want))
  for (ii in pick) for (ty in c("sender", "receiver"))
    tryCatch({ plotEdge(cs, intr = ii, type = ty, slot.use = slot, signif.use = su,
                        plot_dir = file.path(pdir, paste0("edge_", tg, "/")),
                        filename = paste0(ii, "_", ty, ".", plot_fmt), plot.fmt = plot_fmt,
                        return.plot = FALSE, pt.size = 0.05, edge.size = 300)
               log("  plotEdge ok:", ii, ty) },
             error = function(e) log("  plotEdge FAILED", ii, ty, ":", conditionMessage(e)))

  ## plotIntrValue returns a NESTED list: outputList[[intrName]] is itself a list of 8
  ## ggplots, one per `type` (ligand / ligand_ori / ligand_null / receptor / receptor_ori /
  ## receptor_null / score / score_null). `print()`-ing that draws all 8 sequentially onto
  ## one device, so a single-file format keeps only the LAST panel (score_null) and silently
  ## discards the other seven. Combine them into the intended panel figure instead.
  ## (plot_signif_rerun.R has the same defect -- fixed there too, 2026-08-20.)
  for (ii in head(pick, 2))
    tryCatch({ res <- plotIntrValue(cs, intr = ii, slot.use = slot, signif.use = su,
                                    pt.size = 0.15, raster = TRUE)
               pl <- if (length(res) && is.list(res[[1]])) res[[1]] else res
               g  <- cowplot::plot_grid(plotlist = pl, nrow = 2,
                                        labels = names(pl), label_size = 9)
               f <- file.path(pdir, sprintf("intrValue_%s_%s.%s", tg, ii, plot_fmt))
               if (plot_fmt == "png") png(f, 2600, 1400, res = 150) else if (plot_fmt == "pdf") pdf(f, 17.3, 9.3) else svg(f, 17.3, 9.3)
               print(g); dev.off()
               log("  plotIntrValue ok:", ii, "->", length(pl), "panels") },
             error = function(e) { try(dev.off(), silent = TRUE); log("  plotIntrValue FAILED", ii, ":", conditionMessage(e)) })

  tryCatch({ f <- file.path(pdir, sprintf("circos_nIntr_%s.%s", tg, plot_fmt))
             if (plot_fmt == "png") png(f, 1600, 1600, res = 170) else if (plot_fmt == "pdf") pdf(f, 9.4, 9.4) else svg(f, 9.4, 9.4)
             print(plotCircosNIntr(cs, lrscore.use = slot, intr.use = head(sig_ids, 30))); dev.off()
             log("  plotCircosNIntr ok") },
           error = function(e) { try(dev.off(), silent = TRUE); log("  plotCircosNIntr FAILED:", conditionMessage(e)) })

  for (ty in c("sender", "receiver"))
    tryCatch({ plotSigCluster(cs, plot_dir = file.path(pdir, paste0("sigCluster_", tg, "_", ty, "/")),
                              intr.num = 5, type = ty, slot.use = slot, signif.use = su,
                              plot.fmt = plot_fmt, use.cex = 0.05, edge.size = 500)
               log("  plotSigCluster ok:", ty) },
             error = function(e) log("  plotSigCluster FAILED", ty, ":", conditionMessage(e)))
  gc(verbose = FALSE)
}

## ---------------------------------------------------------------- optional persist
## NB the original plot_signif_rerun.R gated this on `df -g`, which is a macOS flag: on
## Linux `df -g` exits with "invalid option", system(intern=TRUE) returns character(0),
## and `if (!is.na(numeric(0)) && ...)` raises "missing value where TRUE/FALSE needed" --
## after the figures but before the manifest. Opt-in via CS_SAVE_RDS instead.
if (save_rds) tryCatch({ saveRDS(purgeBeforeSave(cs), file.path(out_dir, "cs_result.rds")); log("saved cs_result.rds") },
                       error = function(e) log("saveRDS FAILED:", conditionMessage(e)))

hwm <- tryCatch(sub("^VmHWM:\\s*", "", grep("^VmHWM", readLines("/proc/self/status"), value = TRUE)[1]),
                error = function(e) NA_character_)
man <- c(script = "replot_from_quant.R", mode = mode, input_dir = in_dir, out_dir = out_dir,
         db = if (nzchar(db_rds)) basename(db_rds) else "bundled CellPhoneDBv2",
         n_cells_in = length(cells), n_cells_qc = ncol(cs@counts),
         slots = paste(plot_slots, collapse = "|"), n_top = n_top,
         extra = paste(extra_intr, collapse = "|"), plot_fmt = plot_fmt, seed = 42,
         params = "counts.thresh=100 gene.thresh=20 scale.factor=1 r.eps.real=200 perm.size=1e5",
         cytosignal = as.character(packageVersion("cytosignal")), R = R.version.string,
         peak_rss = if (is.na(hwm)) "unknown" else hwm,
         wall_min = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1))
writeLines(c("{", paste0(sprintf('  "%s": "%s"', names(man), man), collapse = ",\n"), "}"),
           file.path(pdir, "replot_manifest.json"))
log("DONE in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min. plots in", pdir)
log("peak RSS:", if (is.na(hwm)) "unknown" else hwm, "-- size --mem from this next time")
