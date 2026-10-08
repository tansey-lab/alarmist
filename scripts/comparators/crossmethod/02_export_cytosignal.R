#!/usr/bin/env Rscript
# ======================================================================================
# Export CytoSignal per-cell LR scores to a form the Python AUC engine can read.
#
# NOTHING IS RECOMPUTED.  The banked quant/score_<slot>.rds objects written by
# run_cytosignal.R are read as-is and re-serialised; no CytoSignal function is called.
#
# WHY DENSE.  The score objects are dgCMatrix but 70-80% dense (P17_AIS: contact
# 0.797, diffusion 0.698), so MatrixMarket text would be ~100M lines for one
# section.  A flat float32 binary is smaller than the .rds and maps straight into
# numpy.  Written COLUMN BY COLUMN, so the file is column-major (Fortran order) --
# the Python side must read it with order='F'.
#
# Slots (CytoSignal's own names, after quant_io.R's gsub):
#   diffusion_Raw_smooth  secreted signalling, diffusion-imputed ligand
#   contact_Raw_smooth    contact-dependent signalling
#   Raw_Raw_smooth        same interactions as `contact`, no imputation (control)
#
# Usage:  Rscript 02_export_cytosignal.R <out_dir> [section ...]
# ======================================================================================
suppressMessages(library(Matrix))

args    <- commandArgs(trailingOnly = TRUE)
out_root <- if (length(args) >= 1) args[[1]] else
    "/data1/tanseyw/projects/fanj2/results/_crossmethod/label_auc/cytosignal_export"
sections <- if (length(args) >= 2) args[-1] else
    c("P17_AIS", "P17_LUAD", "P21_AIS", "P21_LUAD")

QUANT <- "/data1/tanseyw/projects/fanj2/results/cytosignal/LUAD/cellchatdb2"
SLOTS <- c("diffusion_Raw_smooth", "contact_Raw_smooth", "Raw_Raw_smooth")

log <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

for (sec in sections) {
    out_dir <- file.path(out_root, sec)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    cells_written <- NULL

    for (slot in SLOTS) {
        src <- file.path(QUANT, sec, "quant", sprintf("score_%s.rds", slot))
        if (!file.exists(src)) { log(sec, " ", slot, ": MISSING ", src); next }

        t0 <- Sys.time()
        S  <- readRDS(src)
        n  <- nrow(S); p <- ncol(S)

        # Cell order must be identical across slots within a section, so it is
        # written once and checked thereafter.
        if (is.null(cells_written)) {
            writeLines(rownames(S), file.path(out_dir, "cells.tsv"))
            cells_written <- rownames(S)
        } else if (!identical(cells_written, rownames(S))) {
            stop(sprintf("%s/%s: cell order differs from the first slot", sec, slot))
        }
        writeLines(colnames(S), file.path(out_dir, sprintf("intr_%s.tsv", slot)))

        con <- file(file.path(out_dir, sprintf("score_%s.f32", slot)), "wb")
        for (j in seq_len(p)) writeBin(as.numeric(S[, j]), con, size = 4)
        close(con)

        writeLines(sprintf('{"n_cells": %d, "n_intr": %d, "order": "F", "dtype": "float32"}',
                           n, p),
                   file.path(out_dir, sprintf("shape_%s.json", slot)))
        log(sprintf("%-9s %-22s %d cells x %d intr  (%.1f min)", sec, slot, n, p,
                    as.numeric(difftime(Sys.time(), t0, units = "mins"))))
        rm(S); gc(verbose = FALSE)
    }
}
log("done")
