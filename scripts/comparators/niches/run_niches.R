#!/usr/bin/env Rscript
# NICHES on one TMA core of the GBM/LGG Xenium slide, following the NICHES vignettes
# (07 Spatiotemporal for the per-sample spatial call, 01 Spatial for the spatial idioms).
# Every argument is either the vignette's literal value or the package default -- see
# scripts/comparators/niches/NOTES.md for the call-by-call contract.
#
# Usage:
#   Rscript run_niches.R --input <core_dir> --output <out_dir> --lr-db <csv> \
#                        --impute <alra|none> [--force]
#
#   <core_dir> : counts.mtx (genes x cells, raw integer), genes.tsv, barcodes.tsv,
#                meta.csv (cell_id,x,y,celltype,grade,tma_id)   -- from prepare_gbm_input.py
#   <csv>      : 2-column custom_LR_database (ligand,receptor), subunits '_'-separated
#
# Resumable: writes a DONE sentinel and saves every NICHES object as .rds, so a rerun
# skips finished cores outright and downstream work never recomputes the O(N^2) edgelist.
# Env: comp-niches (source scripts/comparators/niches/activate_env.sh first).
suppressWarnings(suppressMessages({
  library(Matrix); library(Seurat); library(NICHES); library(jsonlite)
}))

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(script_dir, "niches_io.R"))

## -------- args --------
args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  args[i + 1]
}
in_dir  <- getarg("--input")  %||% stop("need --input")
out_dir <- getarg("--output") %||% stop("need --output")
lr_csv  <- getarg("--lr-db")  %||% stop("need --lr-db")
impute  <- getarg("--impute", "none")
force   <- "--force" %in% args
stopifnot(impute %in% c("alra", "none"))

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
sentinel <- file.path(out_dir, "DONE")
log <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n")

if (file.exists(sentinel) && !force) {
  log("SKIP -- already complete:", sentinel, "(pass --force to redo)")
  quit(save = "no", status = 0)
}

t0 <- Sys.time()
gc(reset = TRUE)

## -------- load one core --------
log("reading", in_dir)
dge <- as(as(Matrix::readMM(file.path(in_dir, "counts.mtx")), "CsparseMatrix"), "dgCMatrix")
rownames(dge) <- readLines(file.path(in_dir, "genes.tsv"))
colnames(dge) <- readLines(file.path(in_dir, "barcodes.tsv"))
meta <- read.csv(file.path(in_dir, "meta.csv"), stringsAsFactors = FALSE)
stopifnot(all(meta$cell_id == colnames(dge)))
rownames(meta) <- meta$cell_id

# NOTES A1/A2: Seurat object + explicit x/y metadata columns (vignette 01:61-62)
meta.df <- data.frame(
  cell_type = meta$celltype,
  grade     = meta$grade,
  tma_id    = meta$tma_id,
  x         = as.numeric(meta$x),
  y         = as.numeric(meta$y),
  row.names = meta$cell_id,
  stringsAsFactors = FALSE
)
obj <- Seurat::CreateSeuratObject(counts = dge, meta.data = meta.df, assay = "RNA")
log(sprintf("%d genes x %d cells; %d cell types; grade=%s",
            nrow(obj), ncol(obj), length(unique(meta.df$cell_type)), unique(meta.df$grade)[1]))

## -------- NOTES A3: NormalizeData, all defaults (vignette 01:65, 04:46) --------
log("NormalizeData (LogNormalize, scale.factor = 1e4)")
obj <- Seurat::NormalizeData(obj, verbose = FALSE)

## -------- NOTES A4: optional ALRA imputation (vignette 01:72, 07 uses ALRA data) --------
niches_assay <- "RNA"
if (impute == "alra") {
  # SeuratWrappers 0.4.0 will not install in this env: its Imports list pulls in Banksy,
  # whose Bioconductor dependency chain fails to build here (16 packages, incl. rhdf5 /
  # edgeR / HDF5Array / SpatialExperiment). RunALRA itself needs none of that -- only
  # rsvd + Matrix + Seurat, all present. So we source the authors' UNMODIFIED R/alra.R
  # and R/internal.R from seurat-wrappers @ 8df83435cf10ca76f28ddfc658ad3f6bec61824c.
  # Same function, same defaults, same code; see DEVIATIONS.md.
  log("RunALRA (vendored SeuratWrappers alra.R) -- vignette 01:72")
  suppressWarnings(suppressMessages({ library(rsvd) }))
  source(file.path(script_dir, "vendor", "seurat-wrappers", "internal.R"))
  source(file.path(script_dir, "vendor", "seurat-wrappers", "alra.R"))
  obj <- RunALRA(obj)                    # all defaults: k=NULL (auto), q=10, use.mkl=FALSE
  niches_assay <- "alra"
  stopifnot(!is.null(obj@assays[["alra"]]))
  log("alra assay:", nrow(obj@assays[["alra"]]), "genes;",
      "chosen k =", Seurat::Tool(obj, slot = "RunALRA")$k %||% NA)
}

## -------- NOTES A5/A6: RunNICHES --------
lr <- read.csv(lr_csv, stringsAsFactors = FALSE)
stopifnot(ncol(lr) >= 2)
lr <- lr[, 1:2]
log("custom LR database:", basename(lr_csv), nrow(lr), "pairs")

log("RunNICHES: LR.database='custom', assay=", niches_assay,
    ", k=4 (default), blend='mean' (default), CellToCellSpatial=T, NeighborhoodToCell=T, CellToCell=F")
set.seed(42)
niches_out <- NICHES::RunNICHES(
  object             = obj,
  assay              = niches_assay,
  LR.database        = "custom",
  custom_LR_database = lr,
  species            = "human",          # ignored when LR.database='custom' (LoadCustom.R:12)
  cell_types         = "cell_type",
  position.x         = "x",
  position.y         = "y",
  k                  = 4,                # package default
  rad.set            = NULL,             # package default; ignored while k is set
  blend              = "mean",           # package default
  min.cells.per.ident = NULL,            # package default
  min.cells.per.gene  = NULL,            # package default
  meta.data.to.map   = c("orig.ident", "cell_type", "grade", "tma_id", "x", "y"),
  output_format      = "seurat",         # package default
  CellToCellSpatial  = TRUE,             # vignette 07
  NeighborhoodToCell = TRUE,             # vignettes 01 + 07
  CellToCell         = FALSE             # vignette 07
)
log("organizations returned:", paste(names(niches_out), collapse = ", "))

## -------- persist everything --------
shapes <- list()
for (org in names(niches_out)) {
  shapes[[org]] <- save_niches_quant(niches_out[[org]], org, out_dir, log = log)
}

## -------- manifest --------
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
g <- gc()
peak_gb <- sum(g[, "max used"] * c(56, 8)) / 1024^3   # Ncells*56B + Vcells*8B, R heap only
manifest <- list(
  method = "NICHES", version = as.character(packageVersion("NICHES")),
  seurat_version = as.character(packageVersion("Seurat")),
  r_version = paste(R.version$major, R.version$minor, sep = "."),
  dataset = "GBM", tier = "cellchatdb2", imputation = impute,
  input_dir = in_dir, output_dir = out_dir,
  lr_database = normalizePath(lr_csv), n_lr_pairs_input = nrow(lr),
  n_cells = ncol(obj), n_genes = nrow(obj),
  cell_type_column = "cell_type", n_cell_types = length(unique(meta.df$cell_type)),
  grade = unique(meta.df$grade)[1], tma_id = unique(meta.df$tma_id)[1],
  coordinate_units = "microns",
  params = list(assay = niches_assay, LR.database = "custom", k = 4, rad.set = NULL,
                blend = "mean", min.cells.per.ident = NULL, min.cells.per.gene = NULL,
                CellToCellSpatial = TRUE, NeighborhoodToCell = TRUE, CellToCell = FALSE,
                CellToNeighborhood = FALSE, CellToSystem = FALSE, SystemToCell = FALSE),
  seed = 42,
  organizations = shapes,
  wall_minutes = round(elapsed, 2),
  r_heap_peak_gb = round(peak_gb, 2),
  git_sha_niches = "d698e37b8c38ebd103c34acd0b35b03b48d3c5a3"
)
write(jsonlite::toJSON(manifest, auto_unbox = TRUE, pretty = TRUE, null = "null"),
      file.path(out_dir, "run_manifest.json"))

writeLines(format(Sys.time()), sentinel)
log(sprintf("DONE in %.1f min (R heap peak %.1f GB) -> %s", elapsed, peak_gb, out_dir))
