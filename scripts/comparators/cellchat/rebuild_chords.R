#!/usr/bin/env Rscript
# rebuild_chords.R -- per-SENDER chord-gene diagrams (VB:381, VC:397).
#
# Why this exists: the benchmark asks for every plot at full scope, so plot_cellchat.R first
# tries netVisual_chord_gene over ALL sources x ALL targets. On this data that is 161 distinct
# ligand/receptor sectors, and circlize cannot lay that out at ANY gap -- verified at
# small.gap = 1, 0.5, 0.2 and 0.1, all "Maybe your `gap.degree` is too large so that there is
# no space to allocate sectors". The tutorial never draws it at that scope either: VB:381 and
# VC:397 both restrict to a single sender (`sources.use = 4`). So the covering set is one chord
# per sender, which is the tutorial's own call shape and renders for all 9 cell types.
#
# Also redraws the up/down-regulated chords per sender for the same reason.
#
#   source scripts/comparators/cellchat/activate_env.sh
#   Rscript scripts/comparators/cellchat/rebuild_chords.R \
#     --out-dir results/comparators/cellchat/GBM/default --conditions low,high --pos high
suppressWarnings(suppressMessages({ library(CellChat); library(ComplexHeatmap) }))
SCRIPT_DIR <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(SCRIPT_DIR, "cellchat_io.R"))

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(f, d = NULL) { i <- which(args == f); if (!length(i)) d else args[i[1] + 1] }
out_dir    <- getarg("--out-dir")
conditions <- strsplit(getarg("--conditions", "low,high"), ",")[[1]]
pos        <- getarg("--pos", "high")
stopifnot(!is.null(out_dir))
slug <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
notes <- character(0)

for (cond in conditions) {
  obj <- readRDS(file.path(out_dir, "objects", paste0(cond, ".rds")))
  ct <- levels(obj@idents); nct <- length(ct)
  for (i in seq_len(nct)) {
    ok <- save_all_formats(function() netVisual_chord_gene(obj, sources.use = i,
             targets.use = seq_len(nct), lab.cex = 0.4, legend.pos.x = 10,
             title.name = paste0("Signaling from ", ct[i], " - ", cond)),
             file.path(out_dir, "plots", cond, "aggregate",
                       paste0("chord_gene_sender_", slug(ct[i]))), 10, 10)
    if (!isTRUE(ok)) notes <- c(notes, sprintf("FAILED: %s chord_gene_sender_%s", cond, slug(ct[i])))
  }
  cat(cond, ": per-sender chords done\n")
  rm(obj); gc(verbose = FALSE)
}

# up / down regulated, per sender, on the positive-condition object
for (dir_ in c("up", "down")) {
  f <- file.path(out_dir, "quant", sprintf("net_%s_in_%s.csv", dir_, pos))
  if (!file.exists(f)) { cat("no", f, "\n"); next }
  net <- utils::read.csv(f)
  if (nrow(net) == 0) next
  src_cond <- if (dir_ == "up") pos else setdiff(conditions, pos)[1]
  obj <- readRDS(file.path(out_dir, "objects", paste0(src_cond, ".rds")))
  ct <- levels(obj@idents); nct <- length(ct)
  for (i in seq_len(nct)) {
    if (!ct[i] %in% net$source) next
    ok <- save_all_formats(function() netVisual_chord_gene(obj, sources.use = i,
             targets.use = seq_len(nct), slot.name = "net", net = net, lab.cex = 0.5,
             small.gap = 3.5, title.name = paste0(dir_, "-regulated in ", pos, " from ", ct[i])),
             file.path(out_dir, "plots", "comparison",
                       paste0("chord_", dir_, "_sender_", slug(ct[i]))), 10, 10)
    if (!isTRUE(ok)) notes <- c(notes, sprintf("FAILED: chord_%s_sender_%s", dir_, slug(ct[i])))
  }
  cat(dir_, ": per-sender chords done (", nrow(net), "pairs )\n")
  rm(obj); gc(verbose = FALSE)
}

nf <- file.path(out_dir, "plots_not_produced.txt")
existing <- if (file.exists(nf)) readLines(nf) else character(0)
allnotes <- unique(c(existing,
  "netVisual_chord_gene at ALL sources x ALL targets: not renderable -- 161 ligand/receptor sectors exceed circlize's layout capacity at every small.gap tried (1, 0.5, 0.2, 0.1). Covered instead by one chord per sender, which is the tutorial's own scope (VB:381, VC:397).",
  notes))
writeLines(allnotes, nf)
cat("recorded", length(allnotes), "notes in", nf, "\n")
