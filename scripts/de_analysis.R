# Load installed DESeq2 package
library(DESeq2)

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

# Create DESeqDataSet object
dds_obj <- DESeqDataSetFromMatrix(
  countData = reordered_counts_data,
  colData = metadata,
  design = ~ ajcc_pathologic_stage
)
print(head(reordered_counts_data))
print(head(metadata))