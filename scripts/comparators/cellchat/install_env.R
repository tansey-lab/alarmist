#!/usr/bin/env Rscript
# install_env.R -- finish building the comp-cellchat R library and install CellChat itself.
#
# Two things make a plain install.packages() run fail on this box, both recorded as env
# deviations in METHODS.md:
#
#  1. The env is R 4.3.3 (cloned from comp-cytosignal), but CRAN HEAD has moved on: several
#     current sources now require R >= 4.4 (e.g. Deriv 4.2.0 uses Rf_allocLang, which is not
#     API in 4.3.x). So package sources are taken from a DATED CRAN SNAPSHOT contemporaneous
#     with R 4.3.3 rather than from CRAN HEAD.
#  2. Base anaconda's `xml2-config` and headers sit ahead of the env on PATH, so igraph linked
#     against base's libxml2.2.dylib, which the env does not ship -> "loading failed". igraph
#     is therefore built with --disable-graphml (CellChat never reads GraphML), which removes
#     the libxml2 dependency entirely.
#
#   source scripts/comparators/cellchat/activate_env.sh
#   Rscript scripts/comparators/cellchat/install_env.R [--cellchat-src /path/to/CellChat]
SNAPSHOT <- "https://packagemanager.posit.co/cran/2024-06-01"
options(repos = c(CRAN = SNAPSHOT), Ncpus = max(1, parallel::detectCores() - 2))

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(f, d) { i <- which(args == f); if (!length(i)) d else args[i[1] + 1] }
src <- getarg("--cellchat-src", "/Users/jiayifan/tansey_lab/CellChat")

cat("R:", R.version.string, "\nrepo:", SNAPSHOT, "\n")

have <- function(p) requireNamespace(p, quietly = TRUE)
need <- function(p) if (!have(p)) { cat("=== installing", p, "===\n"); TRUE } else FALSE

# --- 1. igraph without GraphML (drops the libxml2 link that fails on this box) ---
# `configure.args = "--disable-graphml"` is NOT enough: igraph is also pulled in as a
# dependency of ggnetwork/Seurat, where the flag is not forwarded. The root cause is that
# igraph's configure probes `xml2-config`, finds BASE anaconda's (the env ships no such
# script), and links -lxml2 against base's libxml2.2.dylib, which is not on the env's rpath.
# Shadowing xml2-config with a failing stub makes configure conclude libxml2 is unavailable
# and build without GraphML I/O -- which CellChat never uses -- for every install path.
shim <- file.path(tempdir(), "no-xml2-bin")
dir.create(shim, showWarnings = FALSE)
writeLines(c("#!/bin/sh", "exit 1"), file.path(shim, "xml2-config"))
Sys.chmod(file.path(shim, "xml2-config"), "0755")
Sys.setenv(PATH = paste(shim, Sys.getenv("PATH"), sep = ":"))
cat("xml2-config now resolves to:", system("command -v xml2-config", intern = TRUE), "\n")

if (need("igraph")) try(install.packages("igraph", configure.args = "--disable-graphml"))

# --- 2. svglite needs systemfonts only; textshaping is a ragg dependency CellChat never uses
for (p in c("systemfonts", "svglite")) if (need(p)) try(install.packages(p))

# --- 3. everything else CellChat declares ---
rest <- c("dplyr", "ggplot2", "future", "future.apply", "pbapply", "irlba", "NMF",
          "ggalluvial", "stringr", "ggrepel", "circlize", "RColorBrewer", "cowplot",
          "RSpectra", "Rcpp", "RcppEigen", "reticulate", "scales", "sna", "reshape2",
          "FNN", "shape", "magrittr", "patchwork", "colorspace", "plyr", "ggpubr",
          "ggnetwork", "plotly", "shiny", "bslib", "collapse", "uwot", "wordcloud",
          "jsonlite", "tidyr", "purrr", "Seurat")
for (p in rest) if (need(p)) try(install.packages(p))

# --- 4. Bioconductor pieces. Biobase is listed FIRST because NMF depends on it and CRAN
#        cannot resolve it -- that is why the NMF install failed on the previous pass.
if (!have("BiocManager")) install.packages("BiocManager")
for (p in c("BiocGenerics", "Biobase", "BiocNeighbors", "ComplexHeatmap")) {
  if (need(p)) try(BiocManager::install(p, ask = FALSE, update = FALSE))
}
# NMF/ggnetwork/Seurat all failed earlier only because igraph or Biobase was absent; retry now
for (p in c("NMF", "ggnetwork", "Seurat")) if (need(p)) try(install.packages(p))

# --- 5. presto: NOT on CRAN. do.fast=TRUE silently falls back to stats::wilcox.test without
#         it and returns systematically larger logFC, while the comparison vignette's
#         thresh.fc = 0.05 was tuned FOR presto (utilities.R:434-445, VC:291).
if (need("presto")) try(remotes::install_github("immunogenomics/presto", upgrade = "never"))

# --- 6. CellChat itself, from the local clone (no network, exact version) ---
if (!have("CellChat")) {
  cat("=== installing CellChat from", src, "===\n")
  try(remotes::install_local(src, dependencies = FALSE, upgrade = "never", force = TRUE))
}

all_needed <- unique(c("igraph", "systemfonts", "textshaping", "svglite", rest,
                       "BiocGenerics", "BiocNeighbors", "ComplexHeatmap", "presto", "CellChat"))
missing <- all_needed[!vapply(all_needed, have, logical(1))]
cat("\n=== STILL MISSING:", if (length(missing)) paste(missing, collapse = " ") else "(none)", "===\n")
if (have("CellChat")) {
  suppressMessages(library(CellChat))
  cat("CellChat", as.character(packageVersion("CellChat")), "loaded OK\n")
  cat("CellChatDB.human$interaction rows:", nrow(CellChatDB.human$interaction), "\n")
}
