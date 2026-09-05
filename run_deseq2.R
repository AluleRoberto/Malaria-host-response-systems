#!/usr/bin/env Rscript

# Differential expression analysis with DESeq2
# for GSE181179 (Malaria host response systems)

# Load required libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tibble)
})

# Set working directory to the repository root (adjust if needed)
# If running via Rscript from repo root, no need.
# If running interactively, setwd("/workspaces/Malaria-host-response-systems")

# ---- 1. Load count matrix and metadata ----
counts_file <- "data/processed/count_matrix_integer.txt"
metadata_file <- "data/processed/sample_metadata_aligned.csv"

counts <- read.delim(counts_file, row.names = 1, check.names = FALSE)
metadata <- read.csv(metadata_file, row.names = 1, stringsAsFactors = FALSE)

# Ensure sample order matches
counts <- counts[, rownames(metadata)]

# Convert group to factor with reference level "Healthy"
metadata$group <- factor(metadata$group, levels = c("Healthy", "Asymptomatic", "Symptomatic"))

# Check dimensions
cat("Count matrix dimensions:", dim(counts), "\n")
cat("Metadata dimensions:", dim(metadata), "\n")
stopifnot(all(colnames(counts) == rownames(metadata)))

# ---- 2. Pre-filtering: remove lowly expressed genes ----
# Keep genes with at least 0.5 CPM in at least 5 samples (smallest group size)
# Compute CPM
cpm <- sweep(counts, 2, colSums(counts), "/") * 1e6
keep <- rowSums(cpm >= 0.5) >= 5
counts_filtered <- counts[keep, ]
cat("Genes before filtering:", nrow(counts), "\n")
cat("Genes after filtering:", nrow(counts_filtered), "\n")

# ---- 3. Create DESeqDataSet ----
dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = metadata,
  design = ~ group
)

# ---- 4. Run DESeq2 pipeline ----
dds <- DESeq(dds)
cat("DESeq2 completed.\n")

# ---- 5. Extract results for pairwise comparisons ----
# We use alpha = 0.1 (FDR threshold) to match your earlier choice,
# but also can report at 0.05 and 0.15 later.

res_sm_vs_hc <- results(dds, contrast = c("group", "Symptomatic", "Healthy"), alpha = 0.1)
res_am_vs_hc <- results(dds, contrast = c("group", "Asymptomatic", "Healthy"), alpha = 0.1)
res_sm_vs_am <- results(dds, contrast = c("group", "Symptomatic", "Asymptomatic"), alpha = 0.1)

# Print summary
cat("\n--- Summary of results ---\n")
cat("Symptomatic vs Healthy:\n"); print(summary(res_sm_vs_hc))
cat("\nAsymptomatic vs Healthy:\n"); print(summary(res_am_vs_hc))
cat("\nSymptomatic vs Asymptomatic:\n"); print(summary(res_sm_vs_am))

# ---- 6. Save results as CSV ----
# Convert to data.frame and add gene IDs as column
res_sm_vs_hc_df <- as.data.frame(res_sm_vs_hc) %>% rownames_to_column("GeneID")
res_am_vs_hc_df <- as.data.frame(res_am_vs_hc) %>% rownames_to_column("GeneID")
res_sm_vs_am_df <- as.data.frame(res_sm_vs_am) %>% rownames_to_column("GeneID")

write.csv(res_sm_vs_hc_df, "data/processed/DESeq2_R_symptomatic_vs_healthy.csv", row.names = FALSE)
write.csv(res_am_vs_hc_df, "data/processed/DESeq2_R_asymptomatic_vs_healthy.csv", row.names = FALSE)
write.csv(res_sm_vs_am_df, "data/processed/DESeq2_R_symptomatic_vs_asymptomatic.csv", row.names = FALSE)

# ---- 7. MA plots and volcano plots ----
# MA plot function
plot_ma <- function(res, title) {
  df <- data.frame(baseMean = res$baseMean, lfc = res$log2FoldChange, sig = res$padj < 0.1)
  ggplot(df, aes(x = log10(baseMean), y = lfc)) +
    geom_point(aes(color = sig), size = 0.8, alpha = 0.5) +
    scale_color_manual(values = c("grey", "red")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x = "log10(baseMean)", y = "log2 fold change", title = title) +
    theme_bw() +
    theme(legend.position = "none")
}

# Volcano plot function
plot_volcano <- function(res, title, fdr_cutoff = 0.1) {
  df <- data.frame(lfc = res$log2FoldChange, pval = -log10(res$padj), sig = res$padj < fdr_cutoff)
  ggplot(df, aes(x = lfc, y = pval)) +
    geom_point(aes(color = sig), size = 0.8, alpha = 0.5) +
    scale_color_manual(values = c("grey", "red")) +
    geom_hline(yintercept = -log10(fdr_cutoff), linetype = "dashed") +
    labs(x = "log2 fold change", y = "-log10 adjusted p-value", title = title) +
    theme_bw() +
    theme(legend.position = "none")
}

# Create and save figures
ggsave("results/figures/MA_sm_vs_am_R.png", plot_ma(res_sm_vs_am, "Symptomatic vs Asymptomatic"), width = 6, height = 5, dpi = 150)
ggsave("results/figures/Volcano_sm_vs_am_R.png", plot_volcano(res_sm_vs_am, "Symptomatic vs Asymptomatic"), width = 6, height = 5, dpi = 150)
# Also for other comparisons if desired
ggsave("results/figures/MA_sm_vs_hc_R.png", plot_ma(res_sm_vs_hc, "Symptomatic vs Healthy"), width = 6, height = 5, dpi = 150)
ggsave("results/figures/Volcano_sm_vs_hc_R.png", plot_volcano(res_sm_vs_hc, "Symptomatic vs Healthy"), width = 6, height = 5, dpi = 150)

# ---- 8. Count significant genes at different FDR thresholds ----
sig_counts <- sapply(list(
  "Symptomatic vs Healthy" = res_sm_vs_hc,
  "Asymptomatic vs Healthy" = res_am_vs_hc,
  "Symptomatic vs Asymptomatic" = res_sm_vs_am
), function(res) {
  c(FDR0.05 = sum(res$padj < 0.05, na.rm = TRUE),
    FDR0.10 = sum(res$padj < 0.10, na.rm = TRUE),
    FDR0.15 = sum(res$padj < 0.15, na.rm = TRUE))
})
print(sig_counts)

# ---- Optional: Save session info for reproducibility ----
writeLines(capture.output(sessionInfo()), "data/processed/R_session_info.txt")