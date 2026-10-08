#!/usr/bin/env Rscript
# run_cellchat.R -- CellChat inference (NOTES.md stages A-D + centrality), one object per
# condition, every core of that condition carried as a `samples` level.
#
# Implements the call table in scripts/comparators/cellchat/NOTES.md one-to-one. Nothing is
# hardcoded: input dir, output dir, tier, and every inference parameter are arguments.
#
#   source scripts/comparators/cellchat/activate_env.sh
#   Rscript scripts/comparators/cellchat/run_cellchat.R \
#     --input-dir results/comparators/cellchat/GBM/input \
#     --out-dir   results/comparators/cellchat/GBM/default \
#     --tier default
suppressWarnings(suppressMessages({
  library(CellChat); library(Matrix); library(jsonlite); library(future); library(BiocNeighbors)
}))
options(stringsAsFactors = FALSE)

SCRIPT_DIR <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(SCRIPT_DIR, "cellchat_io.R"))

# ------------------------------------------------------------------ args
args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  args[i[1] + 1]
}
as_logical_arg <- function(x) toupper(as.character(x)) %in% c("TRUE", "T", "1", "YES")

input_dir  <- getarg("--input-dir")
out_dir    <- getarg("--out-dir")
tier       <- getarg("--tier", "default")
conditions <- strsplit(getarg("--conditions", "low,high"), ",")[[1]]
samples_keep <- getarg("--samples", "")          # smoke test: restrict to these tma_ids
max_cells  <- as.numeric(getarg("--max-cells", "0"))  # smoke test: downsample per condition

interaction_range <- as.numeric(getarg("--interaction-range", "250"))
contact_range     <- as.numeric(getarg("--contact-range", "10"))
distance_use      <- as_logical_arg(getarg("--distance-use", "FALSE"))
scale_distance_in <- getarg("--scale-distance", "auto")
mean_type         <- getarg("--type", "truncatedMean")
trim              <- as.numeric(getarg("--trim", "0.1"))
nboot             <- as.numeric(getarg("--nboot", "100"))
seed              <- as.numeric(getarg("--seed", "1"))
workers           <- as.numeric(getarg("--workers", "4"))
plan_arg          <- getarg("--plan", "multisession")   # multisession | multicore | sequential
min_cells         <- as.numeric(getarg("--min-cells", "10"))
# RESUME. Added 2026-08-20. Default FALSE => behaviour is byte-identical to before.
# The loop at the bottom of this file has no skip logic: it recomputes every condition named
# by --conditions and overwrites its outputs. That is correct for a clean run and destructive
# for a resume -- the 2026-08-20 LUAD run completed AIS in 2166.7 s and was then OOM-killed
# inside the LUAD condition, so a naive re-run would have thrown the good AIS away. With
# --reuse-existing TRUE a condition whose objects/<cond>.rds already exists is LOADED instead
# of recomputed, which also keeps run_manifest.json complete (it is written once, after the
# loop, from every element of `results`). Everything about a reused condition is recomputed
# FROM the saved object -- nothing is copied out of a previous manifest.
reuse_existing    <- as_logical_arg(getarg("--reuse-existing", "FALSE"))

stopifnot(!is.null(input_dir), !is.null(out_dir))
stopifnot(tier %in% c("default", "cellchatdb2"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(out_dir, "run.log")
log <- function(...) {
  msg <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...))
  cat(msg, "\n"); cat(msg, "\n", file = log_file, append = TRUE)
}
log("=== CellChat ", tier, " tier | input=", input_dir, " out=", out_dir, " ===")
log("CellChat ", as.character(packageVersion("CellChat")), " | R ", R.version.string)
presto_available <- rlang::is_installed("presto")
log("presto installed: ", presto_available,
    if (!presto_available) "  <-- do.fast=TRUE will fall back to stats::wilcox.test" else "")

# ------------------------------------------------------------------ Stage B: the DB
# B1-B3 (NOTES.md). The bundled CellChatDB.human IS CellChatDB v2; the tier chooses which
# annotation categories are used, which is the only knob the vignettes expose.
CellChatDB <- CellChatDB.human
CellChatDB.use <- if (tier == "default") {
  subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation")   # V2:131
} else {
  subsetDB(CellChatDB)                                                      # V1:133
}
log("DB tier '", tier, "': ", nrow(CellChatDB.use$interaction), " interactions, ",
    length(unique(CellChatDB.use$interaction$pathway_name)), " pathways; categories: ",
    paste(sort(unique(CellChatDB.use$interaction$annotation)), collapse = ", "))

REQUESTED_LR <- c("GRN_SORT1", "ANXA1_FPR1")
log("requested LRIs in this tier: ",
    paste(sprintf("%s=%s", REQUESTED_LR, REQUESTED_LR %in% CellChatDB.use$interaction$interaction_name),
          collapse = ", "))

# ------------------------------------------------------------------ Stage A: build objects
read_condition <- function(cond) {
  d <- file.path(input_dir, cond)
  log("reading ", d)
  m     <- Matrix::readMM(file.path(d, "data.mtx"))          # genes x cells, log-normalized
  genes <- readLines(file.path(d, "genes.tsv"))
  cells <- readLines(file.path(d, "barcodes.tsv"))
  meta  <- utils::read.csv(file.path(d, "meta.csv"), colClasses = "character")
  sf    <- utils::read.csv(file.path(d, "spatial_factors.csv"), colClasses = "character")
  meta$x <- as.numeric(meta$x); meta$y <- as.numeric(meta$y)
  sf$ratio <- as.numeric(sf$ratio); sf$tol <- as.numeric(sf$tol)
  stopifnot(nrow(meta) == length(cells), identical(meta$cell_id, cells))

  data.input <- as(as(m, "CsparseMatrix"), "dgCMatrix")
  dimnames(data.input) <- list(genes, cells)

  # smoke-test subsetting
  keep <- rep(TRUE, ncol(data.input))
  if (nzchar(samples_keep)) keep <- keep & meta$samples %in% strsplit(samples_keep, ",")[[1]]
  if (max_cells > 0 && sum(keep) > max_cells) {
    set.seed(seed)
    idx <- sample(which(keep), max_cells); keep <- rep(FALSE, length(keep)); keep[idx] <- TRUE
  }
  if (!all(keep)) {
    data.input <- data.input[, keep, drop = FALSE]; meta <- meta[keep, , drop = FALSE]
    sf <- sf[sf$sample %in% unique(meta$samples), , drop = FALSE]
    log("  SUBSET to ", ncol(data.input), " cells / ", nrow(sf), " samples")
  }

  # A2: labels + samples as factors. Unused levels abort computeCommunProb (modeling.R:105).
  meta.cc <- data.frame(labels  = factor(meta$labels),
                        samples = factor(meta$samples),
                        row.names = meta$cell_id)
  meta.cc$labels  <- droplevels(meta.cc$labels)
  meta.cc$samples <- droplevels(meta.cc$samples)

  # A3: coordinates, row-aligned to the columns of data.input
  coords <- data.frame(x = meta$x, y = meta$y, row.names = meta$cell_id)

  # A4: spatial.factors -- ONE ROW PER SAMPLE, in levels(samples) order (modeling.R:1165,1212)
  sf <- sf[match(levels(meta.cc$samples), sf$sample), , drop = FALSE]
  spatial.factors <- data.frame(ratio = sf$ratio, tol = sf$tol, row.names = sf$sample)
  stopifnot(identical(rownames(spatial.factors), levels(meta.cc$samples)))

  list(data.input = data.input, meta = meta.cc, coords = coords, spatial.factors = spatial.factors)
}

# Where do the requested LRIs land? Lifted verbatim out of run_condition on 2026-08-20 so
# that the --reuse-existing path can produce the identical table from a loaded object. It
# reads only slots that survive saveRDS, so a reused condition and a freshly computed one
# give the same answer.
requested_lr_status <- function(cellchat, cond) {
  df.net <- tryCatch(subsetCommunication(cellchat), error = function(e) NULL)
  do.call(rbind, lapply(REQUESTED_LR, function(lr) {
    in_db   <- lr %in% CellChatDB.use$interaction$interaction_name
    in_test <- lr %in% cellchat@LR$LRsig$interaction_name
    hits    <- if (!is.null(df.net)) df.net[df.net$interaction_name == lr, , drop = FALSE] else NULL
    data.frame(condition = cond, interaction_name = lr, in_db = in_db,
               tested = in_test, n_significant_pairs = if (is.null(hits)) 0L else nrow(hits),
               max_prob = if (is.null(hits) || nrow(hits) == 0) NA_real_ else max(hits$prob),
               stringsAsFactors = FALSE)
  }))
}

run_condition <- function(cond) {
  t0 <- Sys.time()

  # ---- resume: load a completed condition instead of recomputing it -------------------
  rds_path <- file.path(out_dir, "objects", sprintf("%s.rds", cond))
  if (reuse_existing && file.exists(rds_path)) {
    log(cond, ": --reuse-existing -> loading ", rds_path, " (NOT recomputing)")
    cellchat <- readRDS(rds_path)
    stopifnot(inherits(cellchat, "CellChat"))
    # d.obs below omits inp$spatial.factors$ratio, which saveRDS does not keep on the object
    # (@images$scale.factors is NULL). For this dataset 04_cellchat.sh passes --ratio 1.0
    # (coordinates already in microns) so the two agree; on a dataset with ratio != 1 the
    # reused d_obs would be in native units. It is a reported sanity check, not an input to
    # any model -- nothing downstream consumes it.
    d.obs <- tryCatch({
      knn <- BiocNeighbors::findKNN(as.matrix(cellchat@images$coordinates), k = 1,
                                    BNPARAM = BiocNeighbors::KmknnParam(), get.index = FALSE)
      min(knn$distance[, 1])
    }, error = function(e) NA_real_)
    # save_rds = FALSE: we just read this file. Rewriting it buys nothing and a kill during
    # the re-serialisation would destroy the very condition we are resuming around.
    stats <- save_cellchat_quant(cellchat, cond, out_dir, log = log, save_rds = FALSE)
    req   <- requested_lr_status(cellchat, cond)
    log(cond, ": reused ", nrow(cellchat@meta), " cells, ",
        nrow(cellchat@LR$LRsig), " LR pairs tested, ",
        nlevels(cellchat@meta$samples), " samples")
    return(list(object = cellchat, stats = stats, wall = NA_real_, req = req,
                scale_distance = cellchat@options$parameter$scale.distance,
                d_obs = d.obs, n_cells = nrow(cellchat@meta),
                n_samples = nlevels(cellchat@meta$samples), reused = TRUE))
  }

  inp <- read_condition(cond)
  log(cond, ": ", ncol(inp$data.input), " cells x ", nrow(inp$data.input), " genes, ",
      nlevels(inp$meta$samples), " samples, ", nlevels(inp$meta$labels), " cell types")

  # observed nearest-neighbour distance -- sanity check on contact.range (VF:35,39).
  # CellChat's own computeCellDistance materialises a dense N x N matrix, which is fine for a
  # Visium section but blows up at 79,998 Xenium cells (~51 GB). BiocNeighbors::findKNN gives
  # the identical quantity -- the minimum cell centre-to-centre distance -- in O(N log N).
  d.obs <- tryCatch({
    knn <- BiocNeighbors::findKNN(as.matrix(inp$coords), k = 1,
                                  BNPARAM = BiocNeighbors::KmknnParam(), get.index = FALSE)
    min(knn$distance[, 1]) * inp$spatial.factors$ratio[1]
  }, error = function(e) NA_real_)
  log(cond, ": observed min cell centre-to-centre distance = ", round(d.obs, 3), " um ",
      "(contact.range = ", contact_range, ")")

  # A5
  cellchat <- createCellChat(object = inp$data.input, meta = inp$meta, group.by = "labels",
                             datatype = "spatial", coordinates = inp$coords,
                             spatial.factors = inp$spatial.factors)
  cellchat@DB <- CellChatDB.use                                             # B4

  # Stage C
  cellchat <- subsetData(cellchat)                                          # C1
  # C2 -- the tutorial's call is future::plan("multisession", workers = 4) (NOTES.md:92).
  # DEFAULT IS UNCHANGED. `--plan multicore` is an opt-in added 2026-08-19 after the LUAD run
  # was OOM-killed at MaxRSS 192.0 GiB against --mem=192G (SLURM job 9051862, oom_kill event).
  #
  # WHY multicore, and why it cannot change a single number:
  #   computeCommunProb calls future_sapply INSIDE the per-LR-pair loop (deparse line 214,
  #   loop closes at 271), i.e. ONCE PER LR PAIR -- 709 times for LUAD/AIS. With multisession
  #   every one of those 709 calls serialises the globals to all 4 PSOCK workers; future's own
  #   error named the payload: "the 11 globals exported is 2.81 GiB ... data.use (2.63 GiB)".
  #   data.use is nSignalingGenes x nCells dense: 2.6 GiB for AIS (475,240 cells) and 6.6 GiB
  #   for LUAD (1,200,922). fork/multicore is copy-on-write, so nothing is serialised at all.
  #   This is a mechanism deviation, not a modelling one. Declare it as such.
  #
  # WHY THE BACKEND CANNOT CHANGE A NUMBER. Full write-up in cellchat/DEVIATIONS.md D-9/D-10;
  # the short version, all MEASURED 2026-08-20:
  #   1. THE ONLY RNG IN computeCommunProb IS IN THE PARENT. deparse 158-159:
  #        set.seed(seed.use); permutation <- replicate(nboot, sample.int(nC, size = nC))
  #      above BOTH my.sapply calls (line 160, once per condition; line 214, once per LR pair).
  #      Both bodies only INDEX it: group[permutation[, nE]].
  #   2. NO RNG IS REACHABLE FROM A WORKER. A namespace-wide grep finds RNG only in
  #      computeCommunProb itself and in computeEnrichmentScore / netAnalysis_contribution /
  #      netVisual* / runPCA / runUMAP. computeExpr_LR, computeExpr_complex,
  #      computeExpr_coreceptor, computeExpr_agonist, computeExpr_antagonist, triMean,
  #      thresholdedMean, computeRegionDistance -- everything the two bodies call -- have none.
  #   3. THE TEST CAN FAIL, so passing it means something: future_sapply(1:8, \(i) runif(1))
  #      returns three DIFFERENT vectors under sequential / multicore / multisession, and
  #      formals(future_sapply)$future.seed is NULL, i.e. future actively detects and warns.
  #   4. NO BLAS CHANNEL EITHER: every crossprod in either body has inner dimension k = 1
  #      (matrix(x, nrow = 1), deparse 169-170 and 227) -- a rank-1 outer product, no
  #      accumulation, so no reduction order to permute. Plus OMP/OPENBLAS/MKL_NUM_THREADS=1.
  #   5. END TO END on a real 4,200-cell / 6-group / 104-LR-pair CellChat object at the
  #      runner's own nboot = 100: sequential vs multisession vs multicore give
  #      identical() == TRUE and max |difference| == 0 for BOTH net$prob and net$pval, in all
  #      three pairwise comparisons. (Coverage caveat: those 104 pairs are all simple Secreted
  #      Signaling, so the agonist / antagonist / co-receptor branches and the contact
  #      P.spatial mutation are covered by 1/2/4 above rather than by the measurement.)
  #
  # THE SCARY WARNING IS NOT ABOUT THIS FUNCTION. Every run, under every plan, emits per
  # condition: "UNRELIABLE VALUE: One of the 'future.apply' iterations ... unexpectedly
  # generated random numbers ... results might be invalid". Per-stage instrumentation puts
  # ZERO of those inside computeCommunProb. All 4 come from netAnalysis_computeCentrality,
  # which has its OWN my.sapply called once per condition, and the single RNG consumer is
  # CellChat:::computeCentralityLocal deparse line 12,
  #     centr$eigen <- igraph::eigen_centrality(G)$vector
  # i.e. ARPACK with a random start vector. It shifts centr$eigen by ~3e-14 and NOTHING reads
  # centr$eigen -- cellchat_io.R:89-92 exports only outdeg/indeg/flowbet/info, plot_cellchat.R
  # never mentions it, and netClustering's `eigen` is base R eigen() on a different matrix.
  # The warning is also present in the macOS/multisession GBM logs, so it did not arrive with
  # this port. Do NOT set future.rng.onMisuse = "ignore": it is the authors' own diagnostic and
  # silencing it would hide a real change if CellChat ever grows one. Point readers at D-10.
  #
  #   workers <= 1 takes an even quieter path: CellChat does
  #     my.sapply <- ifelse(nbrOfWorkers() == 1, sapply, future.apply::future_sapply)
  #   (deparse lines 35-36), so a single worker bypasses `future` ENTIRELY -- no export, no
  #   fork, no serialisation. That is the zero-risk fallback if multicore misbehaves.
  #
  # HAZARD: forking a process that already has live OpenMP/BLAS threads can deadlock. Export
  # OMP_NUM_THREADS=1 before starting R when using multicore; 04_cellchat.sh does.
  plan_used <- local({
    if (!is.finite(workers) || workers <= 1) {
      future::plan("sequential")
      return("sequential (workers<=1: CellChat's my.sapply becomes plain sapply; future bypassed)")
    }
    if (identical(plan_arg, "multicore")) {
      if (isTRUE(parallelly::supportsMulticore())) {
        future::plan("multicore", workers = workers)
        return(paste0("multicore, workers=", workers, " (fork; globals are NOT serialised)"))
      }
      future::plan("multisession", workers = workers)
      return(paste0("multisession, workers=", workers,
                    " (multicore requested but parallelly::supportsMulticore() is FALSE)"))
    }
    future::plan(plan_arg, workers = workers)
    paste0(plan_arg, ", workers=", workers)
  })
  log(cond, ": future plan -> ", plan_used)
  cellchat <- identifyOverExpressedGenes(cellchat)                          # C3
  cellchat <- identifyOverExpressedInteractions(cellchat)                   # C4 (variable.both=TRUE)
  log(cond, ": ", nrow(cellchat@data.signaling), " signaling genes on the panel; ",
      length(cellchat@var.features$features), " over-expressed genes; ",
      nrow(cellchat@LR$LRsig), " over-expressed interactions")
  log(cond, ": requested LRIs surviving identifyOverExpressedInteractions: ",
      paste(sprintf("%s=%s", REQUESTED_LR, REQUESTED_LR %in% cellchat@LR$LRsig$interaction_name),
            collapse = ", "))

  # Stage D1 -- resolve scale.distance if distance.use (modeling.R:152-156)
  scale_distance <- NULL
  if (distance_use) {
    if (scale_distance_in == "auto") {
      d.grp <- tryCatch({
        res <- CellChat:::computeRegionDistance(
          coordinates = cellchat@images$coordinates,
          meta = data.frame(group = cellchat@idents, samples = cellchat@meta$samples,
                            row.names = rownames(cellchat@meta)),
          interaction.range = interaction_range, ratio = inp$spatial.factors$ratio,
          tol = inp$spatial.factors$tol, k.min = 10,
          contact.dependent = TRUE, contact.range = contact_range, contact.knn.k = NULL)
        dm <- res$d.spatial; diag(dm) <- NaN; min(dm, na.rm = TRUE)
      }, error = function(e) NA_real_)
      # CellChat wants min(d * scale.distance) in [1,2]; target 1.5
      scale_distance <- signif(1.5 / d.grp, 3)
      log(cond, ": auto scale.distance = ", scale_distance,
          " (min group distance ", round(d.grp, 2), " um -> scaled ",
          round(d.grp * scale_distance, 2), ")")
    } else {
      scale_distance <- as.numeric(scale_distance_in)
    }
  }

  cellchat <- computeCommunProb(cellchat, type = mean_type, trim = trim,
                                distance.use = distance_use,
                                interaction.range = interaction_range,
                                scale.distance = scale_distance,
                                contact.dependent = TRUE, contact.range = contact_range,
                                nboot = nboot, seed.use = as.integer(seed))
  cellchat <- filterCommunication(cellchat, min.cells = min_cells)          # D2
  cellchat <- computeCommunProbPathway(cellchat)                            # D4
  cellchat <- aggregateNet(cellchat)                                        # D5
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")   # D6

  stats <- save_cellchat_quant(cellchat, cond, out_dir, log = log)
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  log(cond, ": done in ", round(wall, 1), " s")

  req <- requested_lr_status(cellchat, cond)
  list(object = cellchat, stats = stats, wall = wall, req = req,
       scale_distance = scale_distance, d_obs = d.obs,
       n_cells = ncol(inp$data.input), n_samples = nlevels(inp$meta$samples),
       reused = FALSE)
}

results <- list()
for (cond in conditions) results[[cond]] <- run_condition(cond)

# ------------------------------------------------------------------ requested-LR status
req_all <- do.call(rbind, lapply(results, function(r) r$req))
dir.create(file.path(out_dir, "quant"), showWarnings = FALSE, recursive = TRUE)
utils::write.csv(req_all, file.path(out_dir, "quant", "requested_lr_status.csv"), row.names = FALSE)
log("requested-LR status:"); print(req_all)

# ------------------------------------------------------------------ manifest
git_sha <- tryCatch(system("git rev-parse HEAD", intern = TRUE)[1], error = function(e) NA)
manifest <- list(
  method = "CellChat", dataset = "GBM", tier = tier,
  cellchat_version = as.character(packageVersion("CellChat")),
  cellchat_git_sha = "75253cd0c9e68410e6e721a6d3a0419a1d7e358f",
  r_version = R.version.string, alarmist_git_sha = git_sha,
  presto_installed = presto_available,
  input_dir = input_dir, out_dir = out_dir,
  db = list(source = "bundled CellChatDB.human (= CellChatDB v2)",
            subset = if (tier == "default") "subsetDB(search='Secreted Signaling', key='annotation')"
                     else "subsetDB(CellChatDB)  # all but Non-protein Signaling",
            n_interactions = nrow(CellChatDB.use$interaction),
            n_pathways = length(unique(CellChatDB.use$interaction$pathway_name)),
            categories = sort(unique(CellChatDB.use$interaction$annotation))),
  parameters = list(type = mean_type, trim = trim, distance.use = distance_use,
                    interaction.range = interaction_range,
                    scale.distance = results[[1]]$scale_distance,
                    contact.dependent = TRUE, contact.range = contact_range,
                    raw.use = TRUE, population.size = FALSE, k.min = 10,
                    do.symmetric = TRUE, nboot = nboot, seed.use = seed,
                    Kh = 0.5, n = 1, filterCommunication.min.cells = min_cells,
                    variable.both = TRUE, workers = workers),
  smoke = list(samples_keep = samples_keep, max_cells = max_cells),
  conditions = lapply(names(results), function(k) {
    r <- results[[k]]
    list(condition = k, n_cells = r$n_cells, n_samples = r$n_samples,
         n_celltypes = r$stats$n_celltypes, n_lr_tested = r$stats$n_lr_tested,
         n_significant_links = r$stats$n_significant, n_pathways = r$stats$n_pathways,
         observed_min_cell_distance_um = r$d_obs, wall_seconds = round(r$wall, 1),
         # TRUE => loaded from objects/<cond>.rds by --reuse-existing, not recomputed in
         # this process. wall_seconds is then NA: the figure belongs to the earlier run.
         reused = isTRUE(r$reused))
  }),
  # gc()'s LAST column is the "max used" figure in Mb. Index it positionally from the right:
  # this R adds a "limit (Mb)" column, so a hardcoded index 6 lands on the raw object count.
  peak_gc_mb = local({ g <- gc(); round(sum(g[, ncol(g)]), 1) }),
  rss_mb_at_exit = tryCatch(round(as.numeric(system(
    sprintf("ps -o rss= -p %d", Sys.getpid()), intern = TRUE)) / 1024, 1),
    error = function(e) NA),
  finished = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)
writeLines(jsonlite::toJSON(manifest, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null"),
           file.path(out_dir, "run_manifest.json"))
log("wrote ", file.path(out_dir, "run_manifest.json"))
log("=== done ===")
