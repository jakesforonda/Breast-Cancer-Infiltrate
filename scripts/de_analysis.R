# Load the necessary libraries
library(DESeq2)
library(dplyr)
library(pheatmap)
library(httpgd)
library(apeglm)
library(ggplot2)

# Read in merged_unstranded_counts
data_path <- file.path("results", "merged_unstranded_counts.csv")
counts_data <- read.csv(
  data_path, header = TRUE, row.names = 1, check.names = FALSE
)

# Read in metadata
metadata_path <- file.path(
  "metadata", "clinical.cart.2025-01-29", "clinical.tsv"
)

metadata <- read.csv(
  metadata_path, header = TRUE, check.names = FALSE, sep = "\t"
)

# Filter metadata to only include samples in counts_data
metadata <- metadata[, c("case_submitter_id", "ajcc_pathologic_stage")]
metadata <- metadata[!duplicated(metadata), ]
metadata <- metadata[metadata$case_submitter_id %in% colnames(counts_data), ]
rownames(metadata) <- metadata$case_submitter_id
metadata$case_submitter_id <- NULL

# Reorder metadata to match counts_data
reorder_idx <- match(rownames(metadata), colnames(counts_data))
reordered_counts_data <- counts_data[, reorder_idx]

# Clean factor levels in metadata
metadata$ajcc_pathologic_stage <- trimws(
  metadata$ajcc_pathologic_stage
)
metadata$ajcc_pathologic_stage <- gsub(
  "[^A-Za-z0-9_.]", "_", metadata$ajcc_pathologic_stage
)
metadata$ajcc_pathologic_stage <- factor(
  metadata$ajcc_pathologic_stage
)

# Create DESeqDataSet object
dds_obj <- DESeqDataSetFromMatrix(
  countData = reordered_counts_data,
  colData = metadata,
  design = ~ ajcc_pathologic_stage
)

# Perform DESeq2 normalization
dds_obj <- estimateSizeFactors(dds_obj)
dds_obj_norm <- counts(dds_obj, normalized = TRUE)

# Perform DESeq2 analysis and extract results
vsd_obj <- vst(dds_obj, blind = TRUE)
vsd_mat <- assay(vsd_obj)
vsd_cor <- cor(vsd_mat)

# Plot heatmap of correlation matrix
output_file <- file.path("results", "correlation_matrix.png")
png(output_file, width = 1200, height = 800)
pheatmap(
  vsd_cor,
  annotation = select(metadata, ajcc_pathologic_stage),
  main = "Variance Stabilizing Transformation Correlation Heatmap"
)
dev.off()

# Perform PCA analysis
output_file <- file.path("results", "PCA_plot.png")
png(output_file, width = 1200, height = 800)
print(plotPCA(vsd_obj, intgroup = "ajcc_pathologic_stage"), main = "PCA Plot")
dev.off()

# Perform DESeq2 analysis
dds_obj <- DESeq(dds_obj)
plotDispEsts(dds_obj, main = "Dispersion Estimates")

# Extract results
res_stage_2a_vs_1 <- results(
  dds_obj, contrast = c("ajcc_pathologic_stage", "Stage_IIA", "Stage_I")
)
res_stage_3c_vs_1 <- results(
  dds_obj, contrast = c("ajcc_pathologic_stage", "Stage_IIIC", "Stage_I")
)
res_stage_3c_vs_2a <- results(
  dds_obj, contrast = c("ajcc_pathologic_stage", "Stage_IIIC", "Stage_IIA")
)

# Perform shrinkage
res_stage_2a_vs_1 <- lfcShrink(
  dds_obj, coef = "ajcc_pathologic_stage_Stage_IIA_vs_Stage_I", type = "apeglm"
)

res_stage_3c_vs_1 <- lfcShrink(
  dds_obj, coef = "ajcc_pathologic_stage_Stage_IIIC_vs_Stage_I", type = "apeglm"
)

res_stage_3c_vs_2a <- lfcShrink(
  dds_obj,
  contrast = c("ajcc_pathologic_stage", "Stage_IIIC", "Stage_IIA"),
  type = "normal"
)

# Extract significant genes
sig_genes_2a_vs_i <- res_stage_2a_vs_1[
  which(res_stage_2a_vs_1$padj < 0.05),
]
sig_genes_3c_vs_1 <- res_stage_3c_vs_1[
  which(res_stage_3c_vs_1$padj < 0.05),
]
sig_genes_3c_vs_2a <- res_stage_3c_vs_2a[
  which(res_stage_3c_vs_2a$padj < 0.05),
]

# Write results to CSV
output_file1 <- file.path("results", "DE_results_Stage_IIA_vs_I.csv")
output_file2 <- file.path("results", "DE_results_Stage_IIIC_vs_I.csv")
output_file3 <- file.path("results", "DE_results_Stage_IIIC_vs_IIA.csv")
write.csv(sig_genes_2a_vs_i, output_file1)
write.csv(sig_genes_3c_vs_1, output_file2)
write.csv(sig_genes_3c_vs_2a, output_file3)

# Plot MA and Volcano plots

output_file <- file.path("results", "MA_plot_IIA_vs_I.png")
png(output_file, width = 1200, height = 800)
plotMA(res_stage_2a_vs_1, ylim = c(-10, 10),
  main = "Stage IIA vs Stage I"
)
dev.off()

output_file <- file.path("results", "MA_plot_IIIC_vs_I.png")
png(output_file, width = 1200, height = 800)
plotMA(res_stage_3c_vs_1, ylim = c(-10, 10),
  main = "Stage IIIC vs Stage I"
)
dev.off()

output_file <- file.path("results", "MA_plot_IIIC_vs_IIA.png")
png(output_file, width = 1200, height = 800)
plotMA(res_stage_3c_vs_2a, ylim = c(-10, 10),
  main = "Stage IIIC vs Stage IIA"
)
dev.off()

# Function to plot volcano plot
plot_volcano <- function(res, title) {
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df$significant <- res_df$padj < 0.05

  ggplot(res_df,
    aes(x = .data$log2FoldChange,
        y = -log10(.data$padj), color = .data$significant)
  ) +
    geom_point(alpha = 0.7) +
    scale_color_manual(values = c("black", "red")) +
    theme_minimal() +
    labs(title = title, x = "Log2 Fold Change", y = "-Log10 Adjusted P-value")
}

# Plot volcano plots
output_file <- file.path("results", "vol_plot_IIA_vs_I.png")
png(output_file, width = 1200, height = 800)
print(plot_volcano(res_stage_2a_vs_1, "Stage IIA vs Stage I"))
dev.off()

output_file <- file.path("results", "vol_plot_IIIC_vs_I.png")
png(output_file, width = 1200, height = 800)
print(plot_volcano(res_stage_3c_vs_1, "Stage IIIC vs Stage I"))
dev.off()

output_file <- file.path("results", "vol_plot_IIIC_vs_IIA.png")
png(output_file, width = 1200, height = 800)
print(plot_volcano(res_stage_3c_vs_2a, "Stage IIIC vs Stage IIA"))
dev.off()

# Plot heatmap of top 20 significant genes
output_file <- file.path("results", "20_genes_IIA_vs_I.png")
png(output_file, width = 1200, height = 800)
top_genes_2a_vs_1 <- head(order(res_stage_2a_vs_1$padj), 20)
pheatmap(assay(vsd_obj)[top_genes_2a_vs_1, ],
  cluster_rows = TRUE, cluster_cols = TRUE,
  annotation_col = metadata,
  main = "Top 20 Significant Genes Heatmap (Stage IIA vs Stage I)"
)
dev.off()

output_file <- file.path("results", "20_genes_IIIC_vs_I.png")
png(output_file, width = 1200, height = 800)
top_genes_3c_vs_1 <- head(order(res_stage_3c_vs_1$padj), 20)
pheatmap(assay(vsd_obj)[top_genes_3c_vs_1, ],
  cluster_rows = TRUE, cluster_cols = TRUE,
  annotation_col = metadata,
  main = "Top 20 Significant Genes Heatmap (Stage IIIC vs Stage I)"
)
dev.off()

output_file <- file.path("results", "20_genes_IIIC_vs_IIA.png")
png(output_file, width = 1200, height = 800)
top_genes_3c_vs_2a <- head(order(res_stage_3c_vs_2a$padj), 20)
pheatmap(assay(vsd_obj)[top_genes_3c_vs_2a, ],
  cluster_rows = TRUE, cluster_cols = TRUE,
  annotation_col = metadata,
  main = "Top 20 Significant Genes Heatmap (Stage IIIC vs Stage IIA)"
)
dev.off()

# Plot log2 fold change for significant genes
sig_genes <- res_stage_2a_vs_1[which(res_stage_2a_vs_1$padj < 0.05), ]
output_file <- file.path("results", "sig_genes_IIA_vs_I.png")
png(output_file, width = 1200, height = 800)
p <- ggplot(sig_genes,
  aes(x = reorder(rownames(sig_genes), log2FoldChange), y = log2FoldChange)
) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Log2 Fold Change for Significant Genes (Stage IIA vs Stage I)",
    x = "Gene", y = "Log2 Fold Change"
  )
print(p)
dev.off()

sig_genes <- res_stage_3c_vs_1[which(res_stage_3c_vs_1$padj < 0.05), ]
output_file <- file.path("results", "sig_genes_IIIC_vs_I.png")
png(output_file, width = 1200, height = 800)
p <- ggplot(sig_genes,
  aes(x = reorder(rownames(sig_genes), log2FoldChange), y = log2FoldChange)
) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Log2 Fold Change for Significant Genes (Stage IIIC vs Stage I)",
    x = "Gene", y = "Log2 Fold Change"
  )
print(p)
dev.off()

sig_genes <- res_stage_3c_vs_2a[which(res_stage_3c_vs_2a$padj < 0.05), ]
output_file <- file.path("results", "sig_genes_IIIC_vs_IIA.png")
png(output_file, width = 1200, height = 800)
p <- ggplot(sig_genes,
  aes(x = reorder(rownames(sig_genes), log2FoldChange), y = log2FoldChange)
) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    title = "Log2 Fold Change for Significant Genes (Stage IIIC vs Stage IIA)",
    x = "Gene", y = "Log2 Fold Change"
  )
print(p)
dev.off()
