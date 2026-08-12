#!/usr/bin/env Rscript
# audit_db_equivalence.R -- is the BUNDLED CellChatDB.human the same resource as the repo's
# data/LRdatabase/CellChatDBv2.0.human.csv?
#
# This decides whether CellChat's `default` and `cellchatdb2` tiers collapse. The skill says
# to confirm equivalence and SAY SO rather than fake a second run, so this produces the
# evidence: it re-derives the repo CSV's flattening from the bundled DB using exactly the
# mapping CLAUDE.md documents --
#     ligand/receptor <- `complex` subunits joined with "_"   (never the *.symbol columns)
#     pathway <- pathway_name ; signaling_type <- annotation ; version <- version
# -- and then diffs the two key sets.
#
#   source scripts/comparators/cellchat/activate_env.sh
#   Rscript scripts/comparators/cellchat/audit_db_equivalence.R \
#     --csv data/LRdatabase/CellChatDBv2.0.human.csv \
#     --out-dir results/comparators/cellchat/db_audit
suppressWarnings(suppressMessages({ library(CellChat); library(jsonlite) }))
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(f, d = NULL) { i <- which(args == f); if (!length(i)) d else args[i[1] + 1] }
csv_path <- getarg("--csv", "data/LRdatabase/CellChatDBv2.0.human.csv")
out_dir  <- getarg("--out-dir", "results/comparators/cellchat/db_audit")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

DB <- CellChatDB.human
inter <- DB$interaction
cplx  <- DB$complex

# --- flatten the bundled DB the way the repo CSV was exported (CLAUDE.md) ---
expand <- function(name) {
  if (name %in% rownames(cplx)) {
    su <- as.character(cplx[name, ])
    su <- trimws(su); su <- su[nchar(su) > 0 & !is.na(su)]
    paste(su, collapse = "_")
  } else name
}
flat <- data.frame(
  ligand         = vapply(as.character(inter$ligand),   expand, character(1), USE.NAMES = FALSE),
  receptor       = vapply(as.character(inter$receptor), expand, character(1), USE.NAMES = FALSE),
  pathway        = inter$pathway_name,
  signaling_type = inter$annotation,
  version        = inter$version,
  interaction_name = inter$interaction_name,
  stringsAsFactors = FALSE
)
flat$key <- paste0(flat$ligand, "|", flat$receptor)

repo <- utils::read.csv(csv_path)
repo$key <- paste0(repo$ligand, "|", repo$receptor)

only_bundled <- setdiff(unique(flat$key), unique(repo$key))
only_repo    <- setdiff(unique(repo$key), unique(flat$key))
shared       <- intersect(unique(flat$key), unique(repo$key))

# annotation-category agreement on the shared keys
m <- merge(unique(flat[, c("key", "signaling_type", "pathway")]),
           unique(repo[, c("key", "signaling_type", "pathway")]),
           by = "key", suffixes = c("_bundled", "_repo"))
cat_disagree <- m[m$signaling_type_bundled != m$signaling_type_repo, , drop = FALSE]
pw_disagree  <- m[m$pathway_bundled != m$pathway_repo, , drop = FALSE]

utils::write.csv(flat, file.path(out_dir, "bundled_flattened.csv"), row.names = FALSE)
if (length(only_bundled)) utils::write.csv(data.frame(key = only_bundled),
  file.path(out_dir, "only_in_bundled.csv"), row.names = FALSE)
if (length(only_repo)) utils::write.csv(data.frame(key = only_repo),
  file.path(out_dir, "only_in_repo_csv.csv"), row.names = FALSE)
if (nrow(cat_disagree)) utils::write.csv(cat_disagree,
  file.path(out_dir, "signaling_type_disagreements.csv"), row.names = FALSE)
if (nrow(pw_disagree)) utils::write.csv(pw_disagree,
  file.path(out_dir, "pathway_disagreements.csv"), row.names = FALSE)

# per-annotation-category counts, which is what the two tiers actually differ by
tab_bundled <- as.data.frame(table(inter$annotation), stringsAsFactors = FALSE)
colnames(tab_bundled) <- c("annotation", "n_bundled")
tab_repo <- as.data.frame(table(repo$signaling_type), stringsAsFactors = FALSE)
colnames(tab_repo) <- c("annotation", "n_repo_csv")
tab <- merge(tab_bundled, tab_repo, by = "annotation", all = TRUE)
utils::write.csv(tab, file.path(out_dir, "category_counts.csv"), row.names = FALSE)

tier_default <- subsetDB(DB, search = "Secreted Signaling", key = "annotation")
tier_db2     <- subsetDB(DB)

summary <- list(
  cellchat_version = as.character(packageVersion("CellChat")),
  bundled = list(n_rows = nrow(inter), n_unique_keys = length(unique(flat$key)),
                 n_complexes = nrow(cplx), n_cofactors = nrow(DB$cofactor),
                 n_genes_in_geneInfo = nrow(DB$geneInfo)),
  repo_csv = list(path = csv_path, n_rows = nrow(repo), n_unique_keys = length(unique(repo$key))),
  overlap = list(n_shared = length(shared), n_only_bundled = length(only_bundled),
                 n_only_repo = length(only_repo),
                 jaccard = round(length(shared) /
                   length(union(unique(flat$key), unique(repo$key))), 4)),
  disagreements = list(signaling_type = nrow(cat_disagree), pathway = nrow(pw_disagree)),
  category_counts = tab,
  tiers = list(
    default = list(rule = "subsetDB(search='Secreted Signaling', key='annotation')",
                   n_interactions = nrow(tier_default$interaction)),
    cellchatdb2 = list(rule = "subsetDB(CellChatDB)  # all but Non-protein Signaling",
                       n_interactions = nrow(tier_db2$interaction))),
  requested_lr = lapply(c("GRN_SORT1", "ANXA1_FPR1"), function(lr) list(
    interaction_name = lr,
    in_bundled = lr %in% inter$interaction_name,
    in_tier_default = lr %in% tier_default$interaction$interaction_name,
    in_tier_cellchatdb2 = lr %in% tier_db2$interaction$interaction_name,
    annotation = if (lr %in% inter$interaction_name) inter$annotation[match(lr, inter$interaction_name)] else NA,
    pathway = if (lr %in% inter$interaction_name) inter$pathway_name[match(lr, inter$interaction_name)] else NA)),
  # what a `updateCellChatDB` re-import of the flat CSV would COST
  columns_lost_by_reimport = setdiff(colnames(inter), c("ligand", "receptor", "pathway_name",
                                                        "annotation", "version", "interaction_name"))
)
writeLines(jsonlite::toJSON(summary, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null"),
           file.path(out_dir, "db_equivalence.json"))

cat("=== CellChatDB bundled vs repo CSV ===\n")
cat(sprintf("bundled CellChatDB.human : %d rows, %d unique ligand|receptor keys\n",
            nrow(inter), length(unique(flat$key))))
cat(sprintf("repo CSV                 : %d rows, %d unique keys\n", nrow(repo), length(unique(repo$key))))
cat(sprintf("shared %d | only-bundled %d | only-repo %d | Jaccard %.4f\n",
            length(shared), length(only_bundled), length(only_repo), summary$overlap$jaccard))
cat(sprintf("signaling_type disagreements on shared keys: %d ; pathway: %d\n",
            nrow(cat_disagree), nrow(pw_disagree)))
cat("\ncategory counts:\n"); print(tab)
cat(sprintf("\ntier default (Secreted Signaling only): %d interactions\n", nrow(tier_default$interaction)))
cat(sprintf("tier cellchatdb2 (all but Non-protein): %d interactions\n", nrow(tier_db2$interaction)))
cat("\ncolumns a flat-CSV re-import would LOSE:\n  ",
    paste(summary$columns_lost_by_reimport, collapse = ", "), "\n")
cat("\nwrote ", out_dir, "\n")
