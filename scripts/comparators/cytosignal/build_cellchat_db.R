#!/usr/bin/env Rscript
# Build a CytoSignal LR database from alarmist's CellChatDBv2.0.*.csv using the package's
# own formatLRDB(). CytoSignal works in an opaque protein-ID space, so we use GENE SYMBOLS
# as the identifiers (identity g_to_u) — no UniProt mapping needed, and it matches the gene
# symbols in the (already gene-symbol) expression matrix.
#
# signaling_type -> interaction_type: "Cell-Cell Contact" => contact-dependent;
#   everything else (Secreted Signaling, ECM-Receptor, Non-protein Signaling) => diffusion-dependent
#   (mirrors ALARMIST's juxtacrine = Cell-Cell Contact vs secreted split).
#
# Usage: Rscript build_cellchat_db.R <cellchat_csv> <out_rds>
suppressWarnings(suppressMessages(library(cytosignal)))
args <- commandArgs(trailingOnly = TRUE)
csv <- if (length(args) >= 1) args[1] else "data/LRdatabase/CellChatDBv2.0.human.csv"
out <- if (length(args) >= 2) args[2] else "results/comparators/cytosignal/cellchat_db_human.rds"

db <- read.csv(csv, stringsAsFactors = FALSE)
split_genes <- function(s) {                       # split complex on "_" or "," , strip, drop empty
  toks <- trimws(unlist(strsplit(as.character(s), "[_,]")))
  toks[nchar(toks) > 0]
}
lig_list <- lapply(db$ligand,   split_genes)
rec_list <- lapply(db$receptor, split_genes)
itype <- ifelse(db$signaling_type == "Cell-Cell Contact", "contact-dependent", "diffusion-dependent")
cpi <- sprintf("CCI-%05d", seq_len(nrow(db)))

pad <- function(lst, prefix) {                     # list of char vecs -> data.frame of padded columns
  k <- max(lengths(lst))
  m <- t(sapply(lst, function(v) { length(v) <- k; v }))    # NA-pads short ones
  d <- as.data.frame(m, stringsAsFactors = FALSE)
  colnames(d) <- paste0(prefix, seq_len(k)); d
}
interaction_type_df <- data.frame(cpi_cc_id = cpi, interaction_type = itype, stringsAsFactors = FALSE)
ligands_df   <- cbind(cpi_cc_id = cpi, pad(lig_list, "lig"))
receptors_df <- cbind(cpi_cc_id = cpi, pad(rec_list, "rec"))

# identity gene->uniprot over every gene token that appears (ligand or receptor)
all_genes <- sort(unique(c(unlist(lig_list), unlist(rec_list))))
g_to_u <- data.frame(gene_name = all_genes, uniprot = all_genes, hgnc_symbol = all_genes,
                     stringsAsFactors = FALSE)

res <- formatLRDB(interaction_type_df, ligands_df, receptors_df, g_to_u)
# CytoSignal's getLigandNames/getReceptorNames read inter.index cols 4,5 (protein_name_a/b) and fall
# back to cols 2,3 (partner_a/b) when empty. formatLRDB only writes 3 cols -> add empty 4,5 so the
# naming (showIntr return.name, plotSignif, mergeCytoSignal/getIntrNames) works with this custom DB.
ii <- res$intr.index
ii$protein_name_a <- ""; ii$protein_name_b <- ""
ii <- ii[, c("id_cp_interaction", "partner_a", "partner_b", "protein_name_a", "protein_name_b")]
# CytoSignal's getIntrValue (plotSignif) reads the raw db by POSITION: [[2]]=ligand, [[3]]=receptor
# (bundled order is combined,ligands,receptors). formatLRDB returns ligands,receptors,combined, so
# reorder to match — otherwise plotSignif swaps ligand/receptor. (Scoring uses names, so unaffected.)
res$db.diff <- res$db.diff[c("combined", "ligands", "receptors")]
res$db.cont <- res$db.cont[c("combined", "ligands", "receptors")]
out_obj <- list(g_to_u = g_to_u, db.diff = res$db.diff, db.cont = res$db.cont,
                inter.index = ii, source = csv)
saveRDS(out_obj, out)

cat("== CellChat -> CytoSignal DB built ==\n")
cat("rows in CSV:", nrow(db), " unique LR genes:", length(all_genes), "\n")
cat("interactions after dedup: diff =", nlevels(res$db.diff$combined),
    " cont =", nlevels(res$db.cont$combined), " intr.index rows =", nrow(res$intr.index), "\n")
cat("signaling_type -> type:\n"); print(table(db$signaling_type, itype))
cat("saved ->", out, "\n")
