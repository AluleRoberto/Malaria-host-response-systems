#!/usr/bin/env Rscript

# Gene annotation and expression plots for the malaria host-response project.
# Install the required packages once in R if needed.
# install.packages(c("ggplot2", "circlize", "BiocManager"))
# BiocManager::install(c("ComplexHeatmap", "mygene"), update = FALSE)
#
# Run with Rscript 03_gene_annotation_and_plots.R.
# The final annotation step requires internet access to MyGene.info.

project_dir <- "/workspaces/Malaria-host-response-systems"
processed_dir <- file.path(project_dir, "data", "processed")
figures_dir <- file.path(project_dir, "results", "figures")

plot_packages <- c("ggplot2", "ComplexHeatmap", "circlize")
missing_packages <- plot_packages[
  !vapply(plot_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Install the missing plotting packages: ",
       paste(missing_packages, collapse = ", "), call. = FALSE)
}

cohort_colours <- c(
  Healthy = "#008000",
  Asymptomatic = "#0000FF",
  Symptomatic = "#FF0000"
)
cohort_labels <- c(
  Healthy = "Healthy controls",
  Asymptomatic = "Asymptomatic malaria",
  Symptomatic = "Symptomatic malaria"
)

# 1. Load precomputed limma-voom results and align samples.
limma_results <- read.csv(
  file.path(processed_dir, "limma_voom_symptomatic_vs_asymptomatic.csv"),
  check.names = FALSE, colClasses = c(GeneID = "character"),
  stringsAsFactors = FALSE, row.names = NULL
)
count_table <- read.delim(
  file.path(processed_dir, "count_matrix_integer.txt"),
  check.names = FALSE, colClasses = "character", row.names = NULL
)
metadata <- read.csv(
  file.path(processed_dir, "sample_metadata_aligned.csv"),
  check.names = FALSE, colClasses = "character", row.names = NULL
)

required_limma_columns <- c("GeneID", "adj.P.Val", "logFC", "AveExpr")
if (!all(required_limma_columns %in% names(limma_results))) {
  stop("The limma-voom table needs GeneID, adj.P.Val, logFC and AveExpr.",
       call. = FALSE)
}
if (!"group" %in% names(metadata)) {
  stop("The sample metadata needs a group column.", call. = FALSE)
}

sample_ids <- metadata[[1L]]
gene_ids <- count_table[[1L]]
limma_ids <- limma_results[["GeneID"]]

if (anyNA(sample_ids) || any(!nzchar(sample_ids)) ||
    anyDuplicated(sample_ids) > 0L) {
  stop("Sample IDs in the metadata must be present and unique.", call. = FALSE)
}
if (anyNA(gene_ids) || anyNA(limma_ids) ||
    any(!nzchar(gene_ids)) || any(!nzchar(limma_ids)) ||
    anyDuplicated(gene_ids) > 0L || anyDuplicated(limma_ids) > 0L) {
  stop("Gene IDs in the count and limma-voom tables must be present and unique.",
       call. = FALSE)
}

count_samples <- names(count_table)[-1L]
missing_samples <- setdiff(sample_ids, count_samples)
if (length(missing_samples) > 0L) {
  stop("Samples missing from the count matrix: ",
       paste(missing_samples, collapse = ", "), call. = FALSE)
}
if (anyDuplicated(count_samples) > 0L) {
  stop("The count matrix has duplicate sample column names.", call. = FALSE)
}

cohort <- metadata[["group"]]
if (anyNA(cohort) || !all(cohort %in% names(cohort_colours))) {
  stop("Metadata group labels must be Healthy, Asymptomatic or Symptomatic.",
       call. = FALSE)
}

counts <- as.matrix(count_table[, sample_ids, drop = FALSE])
suppressWarnings(storage.mode(counts) <- "double")
dimnames(counts) <- list(gene_ids, sample_ids)
if (any(!is.finite(counts)) || any(counts < 0)) {
  stop("The count matrix must contain finite, nonnegative counts.", call. = FALSE)
}

raw_adjusted_p <- limma_results[["adj.P.Val"]]
adjusted_p <- suppressWarnings(as.numeric(raw_adjusted_p))
if (any(!is.na(raw_adjusted_p) & is.na(adjusted_p))) {
  stop("The adj.P.Val column must contain numbers or missing values.",
       call. = FALSE)
}
limma_results[["adj.P.Val"]] <- adjusted_p

# Rank by BH-adjusted P value. Ties follow the original CSV row order.
ranked_rows <- order(adjusted_p, na.last = NA, method = "radix")

# 2. Filter using full library sizes, then calculate log2(CPM + 1)
# using the filtered library sizes, as in the Python notebook.
library_sizes <- colSums(counts)
if (any(library_sizes <= 0)) {
  stop("Each sample must have a positive library size.", call. = FALSE)
}

cpm <- sweep(counts, 2L, library_sizes, "/") * 1e6
keep_genes <- rowSums(cpm >= 0.5) >= 5L
counts_filtered <- counts[keep_genes, , drop = FALSE]
filtered_sizes <- colSums(counts_filtered)
if (nrow(counts_filtered) < 2L || any(filtered_sizes <= 0)) {
  stop("The CPM filter left too few genes or a sample with no counts.",
       call. = FALSE)
}

log2_cpm <- log2(sweep(counts_filtered, 2L, filtered_sizes, "/") * 1e6 + 1)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("Retained %d genes across %d samples.\n", nrow(log2_cpm), ncol(log2_cpm)))

# 3. Poster heatmap. Select the top 50 before matching to filtered counts.
top50_ids <- limma_ids[head(ranked_rows, 50L)]
matched_ids <- top50_ids[top50_ids %in% rownames(log2_cpm)]
if (length(matched_ids) < 2L) {
  stop("Fewer than two ranked genes match the filtered expression matrix.",
       call. = FALSE)
}
top50_expression <- log2_cpm[matched_ids, , drop = FALSE]

cohort_factor <- factor(cohort, levels = names(cohort_colours))
cohort_counts <- as.integer(table(cohort_factor))
cohort_legend <- sprintf("%s (n=%d)", unname(cohort_labels), cohort_counts)

column_annotation <- ComplexHeatmap::HeatmapAnnotation(
  Cohort = cohort_factor,
  col = list(Cohort = cohort_colours),
  annotation_legend_param = list(
    Cohort = list(at = names(cohort_colours), labels = cohort_legend, nrow = 1L)
  )
)

expression_range <- range(top50_expression)
if (diff(expression_range) == 0) {
  expression_range <- expression_range + c(-0.5, 0.5)
}
expression_colours <- circlize::colorRamp2(
  seq(expression_range[1L], expression_range[2L], length.out = 5L),
  c("#FFF5F0", "#FCBBA1", "#FB6A4A", "#CB181D", "#67000D")
)

poster_heatmap <- ComplexHeatmap::Heatmap(
  top50_expression,
  name = "log2(CPM + 1)",
  col = expression_colours,
  top_annotation = column_annotation,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_columns = "euclidean",
  clustering_method_rows = "average",
  clustering_method_columns = "average",
  row_dend_reorder = FALSE,
  column_dend_reorder = FALSE,
  row_title = "Entrez Gene ID",
  column_title = "PBMC samples",
  column_title_side = "bottom",
  row_names_gp = grid::gpar(fontsize = 8),
  column_names_gp = grid::gpar(fontsize = 8),
  column_names_rot = 90
)

grDevices::png(
  file.path(figures_dir, "heatmap_top50_sm_vs_am.png"),
  width = 11, height = 12, units = "in", res = 300, bg = "white"
)
invisible(tryCatch(
  ComplexHeatmap::draw(
    poster_heatmap,
    column_title = sprintf(
      "%d top-ranked genes (symptomatic vs asymptomatic)", length(matched_ids)
    ),
    heatmap_legend_side = "right",
    annotation_legend_side = "bottom"
  ),
  finally = grDevices::dev.off()
))
cat(sprintf("Heatmap displayed %d of %d selected genes.\n",
            length(matched_ids), length(top50_ids)))

# 4. PCA. Standardise genes across samples using population SD.
sample_expression <- t(log2_cpm)
centred <- sweep(sample_expression, 2L, colMeans(sample_expression), "-")
population_sd <- sqrt(colMeans(centred^2))
population_sd[population_sd == 0] <- 1
scaled_expression <- sweep(centred, 2L, population_sd, "/")

pca <- stats::prcomp(scaled_expression, center = TRUE, scale. = FALSE, rank. = 2L)
explained <- 100 * pca$sdev[1:2]^2 / sum(pca$sdev^2)
pca_data <- data.frame(
  PC1 = pca$x[, 1L],
  PC2 = pca$x[, 2L],
  group = factor(cohort, levels = unique(cohort)),
  row.names = sample_ids
)
pca_plot <- ggplot2::ggplot(pca_data, ggplot2::aes(PC1, PC2, colour = group)) +
  ggplot2::geom_point(size = 3) +
  ggplot2::scale_colour_brewer(palette = "Set1", name = "Group") +
  ggplot2::labs(
    title = sprintf("PCA of filtered log2-CPM (%d genes)", nrow(log2_cpm)),
    x = sprintf("PC1 (%.1f%%)", explained[1L]),
    y = sprintf("PC2 (%.1f%%)", explained[2L])
  ) +
  ggplot2::theme_minimal()
ggplot2::ggsave(
  file.path(figures_dir, "PCA_limma_filtered.png"),
  plot = pca_plot, width = 8, height = 6, units = "in", dpi = 150, bg = "white"
)

# 5. Euclidean sample distances from the log2-CPM matrix.
distance_matrix <- as.matrix(stats::dist(sample_expression, method = "euclidean"))
n_samples <- length(sample_ids)
distance_data <- data.frame(
  sample_x = factor(rep(sample_ids, each = n_samples), levels = sample_ids),
  sample_y = factor(rep(sample_ids, times = n_samples), levels = rev(sample_ids)),
  distance = as.vector(distance_matrix)
)
distance_plot <- ggplot2::ggplot(
  distance_data, ggplot2::aes(sample_x, sample_y, fill = distance)
) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_viridis_c(name = "Euclidean distance") +
  ggplot2::coord_fixed() +
  ggplot2::labs(title = "Sample distance (filtered log2-CPM)", x = NULL, y = NULL) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggplot2::ggsave(
  file.path(figures_dir, "sample_distance_heatmap.png"),
  plot = distance_plot, width = 10, height = 8, units = "in", dpi = 150, bg = "white"
)

# 6. Annotate human Entrez IDs through MyGene.info after saving the plots.
if (!requireNamespace("mygene", quietly = TRUE)) {
  stop("Install mygene with BiocManager::install('mygene') for gene annotation.",
       call. = FALSE)
}

mygene_client <- mygene::MyGene(base.url = "https://mygene.info/v3", verbose = FALSE)
batches <- split(limma_ids, ceiling(seq_along(limma_ids) / 1000L))
gene_records <- unlist(
  lapply(batches, function(ids) {
    mygene::queryMany(
      ids, scopes = "entrezgene", fields = "symbol",
      species = "human", return.as = "records", mygene = mygene_client
    )
  }),
  recursive = FALSE, use.names = FALSE
)

mapped_records <- Filter(function(record) {
  length(record$query) == 1L && length(record$symbol) == 1L &&
    !is.na(record$symbol) && nzchar(record$symbol)
}, gene_records)

query_ids <- vapply(mapped_records, function(record) as.character(record$query),
                   character(1))
symbols <- vapply(mapped_records, function(record) as.character(record$symbol),
                 character(1))
last_hit <- !duplicated(query_ids, fromLast = TRUE)
symbol_lookup <- setNames(symbols[last_hit], query_ids[last_hit])

annotated_results <- limma_results
annotated_results[["Symbol"]] <- unname(symbol_lookup[limma_ids])
write.csv(
  annotated_results,
  file.path(processed_dir, "limma_voom_symptomatic_vs_asymptomatic_with_symbols.csv"),
  row.names = FALSE, na = ""
)

cat(sprintf("\nGenes with adjusted P < 0.05: %d\n",
            sum(adjusted_p < 0.05, na.rm = TRUE)))
top20_columns <- c("GeneID", "Symbol", "logFC", "AveExpr", "adj.P.Val")
print(annotated_results[head(ranked_rows, 20L), top20_columns, drop = FALSE],
      row.names = FALSE)

genes_of_interest <- c(
  "IFNG", "CXCL9", "SLAMF7", "LAG3", "CTLA4",
  "BTG1", "IRF2BP2", "NR1D1", "SIRT1", "SMAD7"
)
candidate_rows <- annotated_results[["Symbol"]] %in% genes_of_interest
cat("\nGenes of interest\n")
print(annotated_results[
  candidate_rows, c("Symbol", "logFC", "adj.P.Val", "GeneID"), drop = FALSE
], row.names = FALSE)
