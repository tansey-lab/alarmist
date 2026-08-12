#!/usr/bin/env Rscript
# CytoSignal-native two-condition (high vs low grade) differential LR test via runNEBULA — ONE script,
# but it spans TWO R installs because no single R here has both cytosignal AND nebula:
#   Stage 1 (this conda R, has cytosignal): build one CytoSignal object per TMA core, mergeCytoSignal,
#           and extract runNEBULA's exact inputs (via cytosignal:::.setup.model).
#   Stage 2 (system R 4.4.2, has the nebula CRAN binary): call nebula::nebula on those inputs — this is
#           runNEBULA's actual computation (NB mixed model, core = subject, ~grade, offset = counts).
# Stage 1 auto-re-invokes this same file in system R for Stage 2 (gated by env var NEBULA_STAGE).
# Usage: source activate_env.sh; Rscript run_nebula_grade.R <input_dir> <db_rds> <out_dir>
args  <- commandArgs(trailingOnly = TRUE)
stage <- Sys.getenv("NEBULA_STAGE", "1")
self  <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
log   <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n")

# system R (nebula binary) — override via env if your paths differ
SYS_R   <- Sys.getenv("SYS_RSCRIPT", "/usr/local/bin/Rscript")
SYS_LIB <- Sys.getenv("SYS_RLIB",    path.expand("~/Library/R/arm64/4.4/library"))

## ============================ STAGE 1 — conda R (cytosignal): build + merge + extract ============================
if (stage != "2") {
  in_dir <- args[1]; db_rds <- args[2]; out_dir <- args[3]
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  suppressWarnings(suppressMessages({ library(Matrix); library(cytosignal) }))

  dge <- as(as(Matrix::readMM(file.path(in_dir, "counts.mtx")), "CsparseMatrix"), "dgCMatrix")
  genes <- readLines(file.path(in_dir, "genes.tsv")); cells <- readLines(file.path(in_dir, "barcodes.tsv"))
  rownames(dge) <- genes; colnames(dge) <- cells
  meta <- read.csv(file.path(in_dir, "meta_grade.csv"), stringsAsFactors = FALSE)
  meta$cell_id <- as.character(meta$cell_id); stopifnot(all(meta$cell_id == cells))
  cdb <- readRDS(db_rds)

  # one CytoSignal object per TMA core (grade is constant within a core); merge does the LR scoring
  cores <- sort(unique(meta$tma_id)); objList <- list(); mrows <- list()
  for (t in cores) {
    idx <- which(meta$tma_id == t)
    loc <- as.matrix(meta[idx, c("x", "y")]); rownames(loc) <- meta$cell_id[idx]; colnames(loc) <- c("x", "y")
    clust <- factor(meta$celltype[idx]); names(clust) <- meta$cell_id[idx]
    log("core", t, ":", length(idx), "cells, grade", meta$grade[idx][1])
    cs <- createCytoSignal(raw.data = dge[, idx, drop = FALSE], cells.loc = loc, clusters = clust)
    cs <- addIntrDB(cs, cdb$g_to_u, cdb$db.diff, cdb$db.cont, cdb$inter.index)
    cs <- removeLowQuality(cs, counts.thresh = 100, gene.thresh = 10)   # lower for small cores -> keep shared interactions
    cs <- changeUniprot(cs)
    cs <- inferEpsParams(cs, scale.factor = 1, r.eps.real = 200)
    cs <- findNN(cs); cs <- imputeLR(cs)          # DT (contact) imputation is required by mergeCytoSignal
    objList[[paste0("core", t)]] <- cs
    mrows[[length(mrows) + 1]] <- data.frame(sample = paste0("core", t), grade = meta$grade[idx][1],
                                              patient = as.character(meta$patient[idx][1]), stringsAsFactors = FALSE)
  }
  metadf <- do.call(rbind, mrows)
  metadf$grade <- factor(metadf$grade, levels = c("low", "high"))   # low = reference => logFC is high-vs-low
  metadf$patient <- factor(metadf$patient)
  cat("== per-core metadata ==\n"); print(metadf)

  log("mergeCytoSignal:", length(objList), "cores")
  merged <- mergeCytoSignal(objList, metadata = metadf, name.by = "sample")
  saveRDS(merged, file.path(out_dir, "merged.rds"))

  # extract runNEBULA's exact inputs (cytosignal:::.setup.model builds the ~grade design used by runNEBULA)
  ms <- cytosignal:::.setup.model(merged, c("grade")); model <- ms[[1]]
  saveRDS(list(diff.lrscore = merged@diff.lrscore, cont.lrscore = merged@cont.lrscore,
               dataset = merged@metadata$dataset, total_counts = merged@metadata[["total_counts"]],
               model = model, cov.use = ms[[2]],
               diff.intr = rownames(merged@diff.lrscore), cont.intr = rownames(merged@cont.lrscore)),
          file.path(out_dir, "nebula_inputs.rds"))
  log("Stage 1 done (merged.rds + nebula_inputs.rds). Re-invoking system R for Stage 2 (nebula)...")

  if (!file.exists(SYS_R)) stop(sprintf("system Rscript not found at %s (set SYS_RSCRIPT); nebula stage cannot run", SYS_R))
  code <- system2(SYS_R, args = c(shQuote(self), shQuote(out_dir)),
                  env = c("R_LIBS=", paste0("R_LIBS_USER=", SYS_LIB), "NEBULA_STAGE=2"))
  quit(status = code)
}

## ============================ STAGE 2 — system R 4.4.2 (nebula): the NB mixed model ============================
suppressWarnings(suppressMessages({ library(Matrix); library(nebula) }))
out_dir <- args[1]
log("Stage 2: nebula", as.character(packageVersion("nebula")))
inp <- readRDS(file.path(out_dir, "nebula_inputs.rds"))

run1 <- function(count, id, model, off) nebula(count, id, pred = model, offset = off, cpc = 0.001, ncore = 4)  # cpc=0.001 per CytoSignal tutorial
dres <- run1(inp$diff.lrscore, inp$dataset, inp$model, inp$total_counts)
cres <- run1(inp$cont.lrscore, inp$dataset, inp$model, inp$total_counts)
saveRDS(list(diff = dres, cont = cres), file.path(out_dir, "nebula_results_raw.rds"))

parse <- function(res, mode) {
  s <- res$summary
  cn <- grep("^logFC_grade", colnames(s), value = TRUE); cov <- sub("^logFC_", "", cn[1])
  df <- data.frame(interaction = s$gene, logFC = s[[paste0("logFC_", cov)]], se = s[[paste0("se_", cov)]],
                   p = s[[paste0("p_", cov)]], mode = mode, stringsAsFactors = FALSE)
  df$padj <- p.adjust(df$p, "BH")   # BH per mode (matches runNEBULA, which p.adjusts diff & cont separately)
  df
}
out <- rbind(parse(dres, "diff"), parse(cres, "cont")); out <- out[order(out$p), ]
write.csv(out, file.path(out_dir, "nebula_grade_results.csv"), row.names = FALSE)
log(sprintf("DONE. runNEBULA: %d interactions | raw p<0.05: %d | BH<0.05: %d -> nebula_grade_results.csv",
            nrow(out), sum(out$p < 0.05, na.rm = TRUE), sum(out$padj < 0.05, na.rm = TRUE)))
print(head(out[, c("interaction", "logFC", "p", "padj")], 10), row.names = FALSE, digits = 3)
