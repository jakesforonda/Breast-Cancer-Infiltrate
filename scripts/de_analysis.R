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
metadata <- metadata[
  metadata$ajcc_pathologic_stage %in% c("Stage I", "Stage IV"),
]
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
print(
  plotPCA(vsd_obj, intgroup = "ajcc_pathologic_stage") + ggtitle("PCA Plot")
)
dev.off()

# Perform DESeq2 analysis
dds_obj <- DESeq(dds_obj)
plotDispEsts(dds_obj, main = "Dispersion Estimates")

# Extract results
res_stage_1_vs_4 <- results(
  dds_obj, contrast = c("ajcc_pathologic_stage", "Stage_I", "Stage_IV")
)

# Perform shrinkage
res_stage_1_vs_4 <- lfcShrink(
  dds_obj,
  contrast = c("ajcc_pathologic_stage", "Stage_I", "Stage_IV"), type = "normal"
)

# Extract significant genes
sig_genes_1_vs_4 <- res_stage_1_vs_4[
  which(res_stage_1_vs_4$padj < 0.05),
]

# Write results to CSV
output_file1 <- file.path("results", "DE_results_Stage_I_vs_IV.csv")
write.csv(sig_genes_1_vs_4, output_file1)

# Plot MA and Volcano plots
output_file <- file.path("results", "MA_plot_I_vs_IV.png")
png(output_file, width = 1200, height = 800)
plotMA(res_stage_1_vs_4, ylim = c(-2, 2),
  main = "Stage I vs Stage IV MA Plot"
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
output_file <- file.path("results", "vol_plot_I_vs_IV.png")
png(output_file, width = 1200, height = 800)
print(plot_volcano(res_stage_1_vs_4, "Stage I vs Stage IV Volcano Plot"))
dev.off()

# Ensure metadata annotation is in the correct format
annotation_col <- data.frame(
  ajcc_pathologic_stage = metadata$ajcc_pathologic_stage
)
rownames(annotation_col) <- rownames(metadata)

# Define output file
output_file <- file.path("results", "20_genes_I_vs_IV.png")
png(output_file, width = 1200, height = 800)

# Extract top significant genes
top_genes_1_vs_4 <- rownames(res_stage_1_vs_4)[
  order(res_stage_1_vs_4$padj, na.last = NA)
][1:min(20, sum(!is.na(res_stage_1_vs_4$padj)))]

# Check if any valid genes exist before plotting
if (length(top_genes_1_vs_4) > 0) {
  pheatmap(
    assay(vsd_obj)[top_genes_1_vs_4, ],
    cluster_rows = TRUE, cluster_cols = TRUE,
    annotation_col = annotation_col,
    main = "Top 20 Significant Genes Heatmap (Stage I vs Stage IV)"
  )
} else {
  message("No significant genes found for heatmap.")
}

dev.off()

# Extract the 20 genes with the largest absolute log2 fold changes
sig_genes <- as.data.frame(sig_genes_1_vs_4)
sig_genes <- sig_genes[!is.na(sig_genes$log2FoldChange), ]
sig_genes$gene <- rownames(sig_genes)  # Store gene names as a column

sorted_indices <- order(
  abs(sig_genes$log2FoldChange), decreasing = TRUE
)
sig_genes_sorted <- sig_genes[sorted_indices, ]

top_genes <- sig_genes_sorted[seq_len(min(20, nrow(sig_genes_sorted))), ]

# Define output file
output_file <- file.path("results", "sig_genes_I_vs_IV.png")
png(output_file, width = 1200, height = 800)

# Create bar plot
p <- ggplot(top_genes, aes(
  x = reorder(gene, log2FoldChange),  # Use the new 'gene' column
  y = log2FoldChange
)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    title = "Top 20 Largest Log2 Fold Changes (Stage I vs Stage IV)",
    x = "Gene", y = "Log2 Fold Change"
  )

print(p)
dev.off()

print("DE analysis complete")