#!/usr/bin/env Rscript

# Differential expression analysis using limma-voom
# Replicating the approach from Studniberg et al. (2022)
# GSE181179 - Malaria host response systems

suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(dplyr)
  library(tibble)
})

# Optional: for gene type filtering, install and load org.Hs.eg.db
# If not installed, skip that step.
has_annotation <- requireNamespace("org.Hs.eg.db", quietly = TRUE)
has_annotation_dbi <- requireNamespace("AnnotationDbi", quietly = TRUE)

# ---- 1. Load data ----
counts_file <- "data/processed/count_matrix_integer.txt"
metadata_file <- "data/processed/sample_metadata_aligned.csv"

counts <- read.delim(counts_file, row.names = 1, check.names = FALSE)
metadata <- read.csv(metadata_file, row.names = 1, stringsAsFactors = FALSE)

# Ensure sample order matches
counts <- counts[, rownames(metadata)]
metadata$group <- factor(metadata$group, levels = c("Healthy", "Asymptomatic", "Symptomatic"))

# ---- 2. Pre-filtering (CPM >= 0.5 in at least 5 samples) ----
cpm <- sweep(counts, 2, colSums(counts), "/") * 1e6
keep_cpm <- rowSums(cpm >= 0.5) >= 5
counts <- counts[keep_cpm, ]
cat("After CPM filter:", nrow(counts), "genes\n")

# ---- Optional: Filter out sex-linked, immunoglobulin, hemoglobin genes ----
if (has_annotation && has_annotation_dbi) {
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  # Map Entrez IDs to gene symbols and chromosomes
  gene_info <- select(org.Hs.eg.db, keys = rownames(counts),
                      columns = c("SYMBOL", "CHR"), keytype = "ENTREZID")
  # Identify genes to remove
  remove_genes <- gene_info$ENTREZID[
    grepl("^IG[HKL]", gene_info$SYMBOL) |      # immunoglobulins
    grepl("^HB[AB]", gene_info$SYMBOL) |       # hemoglobin
    gene_info$CHR %in% c("X", "Y")             # sex chromosomes
  ]
  remove_genes <- unique(remove_genes[!is.na(remove_genes)])
  cat("Removing", length(remove_genes), "sex-linked/immunoglobulin/hemoglobin genes\n")
  counts <- counts[!(rownames(counts) %in% remove_genes), ]
} else {
  cat("Annotation package org.Hs.eg.db not available; skipping gene-type filtering.\n")
}

cat("Final gene count:", nrow(counts), "\n")

# ---- 3. Create DGEList and apply TMM normalization ----
dge <- DGEList(counts = counts, group = metadata$group)
dge <- calcNormFactors(dge)  # TMM normalization

# ---- 4. Voom transformation with quantile normalization ----
design <- model.matrix(~ 0 + group, data = metadata)
colnames(design) <- levels(metadata$group)  # "Healthy", "Asymptomatic", "Symptomatic"

v <- voom(dge, design, normalize.method = "quantile", plot = FALSE)

# ---- 5. Fit linear model and compute contrasts ----
fit <- lmFit(v, design)

# Define contrasts for pairwise comparisons
contrast_matrix <- makeContrasts(
  Symptomatic_vs_Healthy = Symptomatic - Healthy,
  Asymptomatic_vs_Healthy = Asymptomatic - Healthy,
  Symptomatic_vs_Asymptomatic = Symptomatic - Asymptomatic,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# ---- 6. Decide tests with global FDR (like original paper) ----
# The original used decideTests with method="global" and an FDR cutoff of 15%.
# We'll use 15% as in the paper, but also report 5% and 10% for comparison.
for (fdr_cut in c(0.05, 0.10, 0.15)) {
  cat("\n--- FDR cutoff:", fdr_cut, "---\n")
  results <- decideTests(fit2, method = "global", p.value = fdr_cut, lfc = 0)
  print(summary(results))
}

# ---- 7. Extract results for each contrast at FDR 0.15 (to compare with original) ----
# We'll create a table of genes with adjusted p-values and logFC for each contrast.
# limma's topTable can be used per contrast, but we want the global FDR.
# The global FDR adjustment is not stored in topTable; we'll compute it manually.
# First, get raw p-values for each contrast.
pvals <- fit2$p.value
# Adjust p-values globally (across all contrasts) using BH.
# The 'global' method in decideTests essentially applies BH across all contrasts jointly.
# We can replicate by stacking p-values and adjusting, then reshaping.
# But for simplicity, we'll just extract results using topTable per contrast with BH (default)
# and note that the global method is more conservative.
# For exact replication, we can use the adjusted p-values from decideTests?
# Actually decideTests returns classification, not adjusted p-values.
# We'll use topTable with adjust.method="BH" and sort by p-value,
# and then filter genes that are significant in the global classification.

# Let's get the global classification matrix at FDR 0.15.
global_class <- decideTests(fit2, method = "global", p.value = 0.15, lfc = 0)
# The row names are genes, columns are contrasts with -1, 0, 1.

# For each contrast, get the significant genes (up or down)
sig_sm_vs_hc <- rownames(global_class)[global_class[, "Symptomatic_vs_Healthy"] != 0]
sig_am_vs_hc <- rownames(global_class)[global_class[, "Asymptomatic_vs_Healthy"] != 0]
sig_sm_vs_am <- rownames(global_class)[global_class[, "Symptomatic_vs_Asymptomatic"] != 0]

cat("\nNumber of significant genes (global FDR<0.15):\n")
cat("Symptomatic vs Healthy:", length(sig_sm_vs_hc), "\n")
cat("Asymptomatic vs Healthy:", length(sig_am_vs_hc), "\n")
cat("Symptomatic vs Asymptomatic:", length(sig_sm_vs_am), "\n")

# ---- 8. Save full results tables (with BH adjusted p-values per contrast) ----
# For reproducibility and later pathway analysis, save the topTable results per contrast.
res_sm_vs_hc <- topTable(fit2, coef = "Symptomatic_vs_Healthy", number = Inf, sort.by = "none")
res_am_vs_hc <- topTable(fit2, coef = "Asymptomatic_vs_Healthy", number = Inf, sort.by = "none")
res_sm_vs_am <- topTable(fit2, coef = "Symptomatic_vs_Asymptomatic", number = Inf, sort.by = "none")

# Add gene ID column
res_sm_vs_hc <- rownames_to_column(res_sm_vs_hc, "GeneID")
res_am_vs_hc <- rownames_to_column(res_am_vs_hc, "GeneID")
res_sm_vs_am <- rownames_to_column(res_sm_vs_am, "GeneID")

write.csv(res_sm_vs_hc, "data/processed/limma_voom_symptomatic_vs_healthy.csv", row.names = FALSE)
write.csv(res_am_vs_hc, "data/processed/limma_voom_asymptomatic_vs_healthy.csv", row.names = FALSE)
write.csv(res_sm_vs_am, "data/processed/limma_voom_symptomatic_vs_asymptomatic.csv", row.names = FALSE)

# ---- 9. Basic plots ----
# MA plot for key contrast
pdf("results/figures/MA_limma_sm_vs_am.pdf", width = 6, height = 5)
limma::plotMD(fit2, coef = "Symptomatic_vs_Asymptomatic", status = decideTests(fit2, method = "global", p.value = 0.15)[, "Symptomatic_vs_Asymptomatic"],
              main = "Symptomatic vs Asymptomatic (limma-voom)")
dev.off()

# Volcano plot (using base R)
pdf("results/figures/Volcano_limma_sm_vs_am.pdf", width = 6, height = 5)
plot(res_sm_vs_am$logFC, -log10(res_sm_vs_am$adj.P.Val),
     xlab = "log2 fold change", ylab = "-log10 adjusted p-value",
     main = "Volcano: Symptomatic vs Asymptomatic",
     pch = 20, cex = 0.4, col = ifelse(res_sm_vs_am$adj.P.Val < 0.15, "red", "grey"))
abline(h = -log10(0.15), lty = 2)
dev.off()

# Save session info
writeLines(capture.output(sessionInfo()), "data/processed/limma_session_info.txt")
cat("\nScript completed successfully.\n")