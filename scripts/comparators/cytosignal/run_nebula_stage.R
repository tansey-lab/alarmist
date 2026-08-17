#!/usr/bin/env Rscript
# ======================================================================================
# CytoSignal's NATIVE cross-condition differential test on the LUAD/AIS four-section set:
# AIS vs LUAD via mergeCytoSignal + runNEBULA.
#
# This is the LUAD sibling of run_nebula_grade.R (the signed-off GBM record, which is NOT
# touched). Differences from it, all forced by the data:
#   - the sample unit is a SECTION (4 of them), not a TMA core (13);
#   - the covariate is `stage` (AIS reference => logFC reads LUAD-vs-AIS), not `grade`;
#   - inputs are four per-section directories written by cytosignal/export_cs_input.py,
#     not one directory plus a meta_grade.csv;
#   - gene.thresh = 20, matching run_cytosignal.R. run_nebula_grade.R lowered it to 10 as a
#     concession to small TMA cores; LUAD sections are 182k-641k cells and need no such
#     concession, so the consistent value is used.
#
# TWO R INSTALLS, exactly as in run_nebula_grade.R, because `nebula` does not build inside
# comp-cytosignal:
#   Stage 1 (comp-cytosignal): build one CytoSignal object per section, mergeCytoSignal,
#           extract runNEBULA's exact inputs via cytosignal:::.setup.model.
#   Stage 2 (an R that HAS nebula): nebula::nebula on those inputs -- the NB mixed model
#           with section as the random effect and total counts as offset.
# Stage 1 re-invokes this same file in the second R, gated by env var NEBULA_STAGE.
# Point SYS_RSCRIPT / SYS_RLIB at that second R (on iris: the comp-nebula env).
#
# WHY THIS PATH FITS IN MEMORY WHEN run_cytosignal.R DOES NOT:
#   mergeCytoSignal never calls inferIntrScore. inferIntrScore is what floors perm.size at
#   ncol(dge.raw) (see cytosignal:::permuteLR), which is what peaked at ~57 GB on the full
#   P21_LUAD section. The merge path does findNN + imputeLR + inferScoreLR only.
#
# KNOWN RISK, READ BEFORE RUNNING: the merged LRscore matrix is n_cells x n_interactions =
# 1,676,162 x ~2,683 = ~4.5e9 elements. A dgCMatrix stores its @i slot as int32, capped at
# 2,147,483,647 non-zeros. If the LRscore matrix is near-dense (METHODS.md records that it
# is), the merge will fail with a length/allocation error. If that happens, the fallback is
# to merge PAIRWISE within patient (P17_AIS+P17_LUAD, then P21_AIS+P21_LUAD) and report two
# within-patient contrasts instead of one pooled one -- which is arguably the better design
# anyway, since the four sections are a matched 2x2. --pairs does exactly that.
#
# STATISTICAL HEALTH WARNING: the random effect has FOUR levels and `stage` is a
# section-level covariate with 2 vs 2. The ~1.5M cell-level observations will make p-values
# look far more confident than the design supports. Treat the output as exploratory
# material, not as a test.
#
# Usage:
#   source scripts/comparators/cytosignal/activate_env.sh
#   Rscript run_nebula_stage.R --input-root <dir with one subdir per section> \
#                              --db <cellchat_db_human.rds> --out-dir <dir> [--pairs]
# ======================================================================================

args  <- commandArgs(trailingOnly = TRUE)
stage <- Sys.getenv("NEBULA_STAGE", "1")
self  <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
log   <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n")
getarg <- function(f, d = NULL) { i <- which(args == f); if (!length(i)) d else args[i[1] + 1] }

SYS_R   <- Sys.getenv("SYS_RSCRIPT", "/usr/local/bin/Rscript")
SYS_LIB <- Sys.getenv("SYS_RLIB",    path.expand("~/Library/R/arm64/4.4/library"))

# Section order defines the factor levels; AIS first so it is the reference.
SECTIONS <- c("P17_AIS", "P17_LUAD", "P21_AIS", "P21_LUAD")

## ================= STAGE 1 — comp-cytosignal: build + merge + extract =================
if (stage != "2") {
  in_root <- getarg("--input-root"); db_rds <- getarg("--db"); out_dir <- getarg("--out-dir")
  do_pairs <- "--pairs" %in% args
  gene_thresh   <- as.numeric(getarg("--gene-thresh", "20"))
  counts_thresh <- as.numeric(getarg("--counts-thresh", "100"))
  if (is.null(in_root) || is.null(db_rds) || is.null(out_dir))
    stop("need --input-root, --db and --out-dir")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  suppressWarnings(suppressMessages({ library(Matrix); library(cytosignal) }))

  cdb <- readRDS(db_rds)

  build_one <- function(sec) {
    d <- file.path(in_root, sec)
    log("building ", sec, " from ", d)
    dge <- as(as(Matrix::readMM(file.path(d, "counts.mtx")), "CsparseMatrix"), "dgCMatrix")
    genes <- readLines(file.path(d, "genes.tsv"))
    cells <- readLines(file.path(d, "barcodes.tsv"))
    rownames(dge) <- genes; colnames(dge) <- cells
    meta <- read.csv(file.path(d, "meta.csv"), stringsAsFactors = FALSE)
    meta$cell_id <- as.character(meta$cell_id)
    stopifnot(identical(meta$cell_id, cells))

    loc <- as.matrix(meta[, c("x", "y")])
    rownames(loc) <- meta$cell_id; colnames(loc) <- c("x", "y")
    clust <- factor(meta$celltype); names(clust) <- meta$cell_id

    cs <- createCytoSignal(raw.data = dge, cells.loc = loc, clusters = clust)
    cs <- addIntrDB(cs, cdb$g_to_u, cdb$db.diff, cdb$db.cont, cdb$inter.index)
    cs <- removeLowQuality(cs, counts.thresh = counts_thresh, gene.thresh = gene_thresh)
    cs <- changeUniprot(cs)
    cs <- inferEpsParams(cs, scale.factor = 1, r.eps.real = 200)  # microns already
    cs <- findNN(cs)
    cs <- imputeLR(cs)      # the DT imputation mergeCytoSignal requires; not optional
    log("  ", sec, ": ", ncol(cs@counts), " cells after QC")
    list(cs = cs,
         row = data.frame(sample  = sec,
                          stage   = meta$stage[1],
                          patient = meta$patient[1],
                          stringsAsFactors = FALSE))
  }

  merge_and_extract <- function(secs, tag) {
    objList <- list(); mrows <- list()
    for (s in secs) {
      b <- build_one(s)
      objList[[s]] <- b$cs
      mrows[[length(mrows) + 1]] <- b$row
    }
    metadf <- do.call(rbind, mrows)
    metadf$stage   <- factor(metadf$stage, levels = c("AIS", "LUAD"))  # AIS = reference
    metadf$patient <- factor(metadf$patient)
    cat("== metadata for ", tag, " ==\n", sep = ""); print(metadf)

    log("mergeCytoSignal: ", length(objList), " sections (", tag, ")")
    merged <- mergeCytoSignal(objList, metadata = metadf, name.by = "sample")
    saveRDS(merged, file.path(out_dir, sprintf("merged_%s.rds", tag)))

    ms <- cytosignal:::.setup.model(merged, c("stage"))
    saveRDS(list(diff.lrscore = merged@diff.lrscore, cont.lrscore = merged@cont.lrscore,
                 dataset = merged@metadata$dataset,
                 total_counts = merged@metadata[["total_counts"]],
                 model = ms[[1]], cov.use = ms[[2]],
                 diff.intr = rownames(merged@diff.lrscore),
                 cont.intr = rownames(merged@cont.lrscore)),
            file.path(out_dir, sprintf("nebula_inputs_%s.rds", tag)))
    log("  wrote merged_", tag, ".rds + nebula_inputs_", tag, ".rds")
    rm(merged, objList); gc()
  }

  if (do_pairs) {
    log("--pairs: two WITHIN-PATIENT merges instead of one pooled merge")
    merge_and_extract(c("P17_AIS", "P17_LUAD"), "P17")
    merge_and_extract(c("P21_AIS", "P21_LUAD"), "P21")
    tags <- c("P17", "P21")
  } else {
    merge_and_extract(SECTIONS, "pooled")
    tags <- "pooled"
  }
  writeLines(tags, file.path(out_dir, "tags.txt"))

  log("Stage 1 done. Re-invoking ", SYS_R, " for Stage 2 (nebula) ...")
  if (!file.exists(SYS_R))
    stop(sprintf(paste0("second R not found at %s. Set SYS_RSCRIPT/SYS_RLIB to an R that has ",
                        "the `nebula` package (on iris: the comp-nebula env). Stage 1 outputs ",
                        "are already on disk, so you can re-run Stage 2 alone with ",
                        "NEBULA_STAGE=2 <Rscript> %s %s"), SYS_R, self, out_dir))
  code <- system2(SYS_R, args = c(shQuote(self), shQuote(out_dir)),
                  env = c("R_LIBS=", paste0("R_LIBS_USER=", SYS_LIB), "NEBULA_STAGE=2"))
  quit(status = code)
}

## ================= STAGE 2 — the R that has nebula: the NB mixed model =================
suppressWarnings(suppressMessages({ library(Matrix); library(nebula) }))
out_dir <- args[1]
log("Stage 2: nebula ", as.character(packageVersion("nebula")))
tags <- readLines(file.path(out_dir, "tags.txt"))

run1 <- function(count, id, model, off)
  nebula(count, id, pred = model, offset = off, cpc = 0.001, ncore = 4)  # cpc per the vignette

parse <- function(res, mode) {
  s <- res$summary
  cn <- grep("^logFC_stage", colnames(s), value = TRUE)
  if (!length(cn)) cn <- grep("^logFC_", colnames(s), value = TRUE)[2]
  cov <- sub("^logFC_", "", cn[1])
  df <- data.frame(interaction = s$gene, logFC = s[[paste0("logFC_", cov)]],
                   se = s[[paste0("se_", cov)]], p = s[[paste0("p_", cov)]],
                   mode = mode, covariate = cov, stringsAsFactors = FALSE)
  df$padj <- p.adjust(df$p, "BH")   # BH per mode, matching runNEBULA
  df
}

for (tag in tags) {
  inp <- readRDS(file.path(out_dir, sprintf("nebula_inputs_%s.rds", tag)))
  log("nebula on ", tag, ": ", length(inp$diff.intr), " diff / ",
      length(inp$cont.intr), " cont interactions")
  dres <- run1(inp$diff.lrscore, inp$dataset, inp$model, inp$total_counts)
  cres <- run1(inp$cont.lrscore, inp$dataset, inp$model, inp$total_counts)
  saveRDS(list(diff = dres, cont = cres),
          file.path(out_dir, sprintf("nebula_results_raw_%s.rds", tag)))

  out <- rbind(parse(dres, "diff"), parse(cres, "cont"))
  out <- out[order(out$p), ]
  f <- file.path(out_dir, sprintf("nebula_stage_results_%s.csv", tag))
  write.csv(out, f, row.names = FALSE)
  log(sprintf("  %s: %d interactions | raw p<0.05: %d | BH<0.05: %d -> %s",
              tag, nrow(out), sum(out$p < 0.05, na.rm = TRUE),
              sum(out$padj < 0.05, na.rm = TRUE), basename(f)))
  print(head(out[, c("interaction", "logFC", "p", "padj")], 10), row.names = FALSE, digits = 3)
}
log("DONE. logFC is LUAD relative to AIS (AIS is the reference level).")
