#!/usr/bin/env Rscript
# Match ALARMIST motif-1 top LRIs to the CytoSignal per-interaction grade (high-vs-low) test results.
suppressWarnings(suppressMessages(library(cytosignal)))
d  <- "results/comparators/cytosignal/GBM/cellchatdb2/nebula_grade"
cdb <- readRDS("results/comparators/cytosignal/cellchat_db_human.rds"); ii <- cdb$inter.index
setkey <- function(a) vapply(strsplit(a, "_"), function(x) paste(sort(x), collapse = "|"), character(1))
ii$name    <- paste0(ii$partner_a, "-", ii$partner_b)                     # getIntrNames format
ii$lig_set <- setkey(ii$partner_a)
ii$rec_set <- setkey(sub("-receptor$", "", ii$partner_b))
key2name <- setNames(ii$name, paste(ii$lig_set, ii$rec_set))

lm <- read.csv("results/GBM/bptf/lri_motifs.csv"); m1 <- lm[lm$motif_idx == 1, ]
m1$lig_set <- setkey(m1$ligand); m1$rec_set <- setkey(m1$receptor)
m1u <- m1[order(-m1$factor_lrnorm), ]; m1u <- m1u[!duplicated(paste(m1u$lig_set, m1u$rec_set)), ]
m1u$cs_name <- key2name[paste(m1u$lig_set, m1u$rec_set)]

allg <- read.csv(file.path(d, "nebula_grade_results.csv"))               # runNEBULA grade results (interaction, logFC, p, padj, mode)
m1u <- merge(m1u, allg, by.x = "cs_name", by.y = "interaction", all.x = TRUE, sort = FALSE)
m1u <- m1u[order(-m1u$factor_lrnorm), ]
m1u$alarmist_rank <- seq_len(nrow(m1u))

out <- m1u[, c("alarmist_rank", "ligand", "receptor", "celltype1", "celltype2", "signaling_type",
               "factor_lrnorm", "cs_name", "mode", "logFC", "p", "padj")]
write.csv(out, file.path(d, "motif1_top_lris_grade.csv"), row.names = FALSE)

cat(sprintf("motif-1 unique top LR pairs: %d | matched to a CytoSignal grade test: %d\n",
            nrow(m1u), sum(!is.na(m1u$p))))
cat(sprintf("  ... with grade p<0.05: %d | padj<0.05: %d\n",
            sum(m1u$p < 0.05, na.rm = TRUE), sum(m1u$padj < 0.05, na.rm = TRUE)))
cat("\n== motif-1 top 25 LRIs and their CytoSignal grade high-vs-low result ==\n")
print(head(out[, c("alarmist_rank", "ligand", "receptor", "signaling_type", "factor_lrnorm",
                   "logFC", "p", "padj")], 25), row.names = FALSE, digits = 3)
