# Load the necessary libraries
library(DESeq2)
library(dplyr)
library(pheatmap)
library(httpgd)

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

hgd()

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
pheatmap(
  vsd_cor,
  annotation = select(metadata, ajcc_pathologic_stage),
  main = "Variance Stabilizing Transformation Correlation Heatmap"
)

# Perform PCA analysis
print(plotPCA(vsd_obj, intgroup = "ajcc_pathologic_stage"))

hgd_view()

# Perform DESeq2 analysis
dds_obj <- DESeq(dds_obj)
plotDispEsts(dds_obj)