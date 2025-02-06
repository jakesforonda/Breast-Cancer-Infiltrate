# Load data.table package
library(data.table)

# Make a list of files to access for iterating -> remove head() after testing
files <- list.files(path = "data/", pattern = "*.tsv", full.names = FALSE)

# Read in metadata file
metadata <- fread("metadata/gdc_sample_sheet.2025-01-29.tsv")
metadata <- metadata[, c("File Name", "Case ID")]


# Make empty lists to collect data.tables for analyses
de_list <- list()
infil_list <- list()


print("Processing TSV Files...")

# Read every file and process columns for analyses
for (file in files) {
  # Full path to file
  full_path <- paste("data/", file, sep = "")

  # Read file as a data.table
  data <- fread(full_path)

  # Only keep rows that have protein_coding genes
  data <- data[gene_type == "protein_coding"]

  # Process infiltrate data
  infil_tbl <- data[, .(gene_name, tpm_unstranded)]
  infil_tbl[, gene_name := sub("\\..*", "", gene_name)]
  infil_tbl <- infil_tbl[
    , .(tpm_unstranded = as.integer(round(sum(tpm_unstranded)))),
    by = gene_name
  ]
  setnames(infil_tbl, "tpm_unstranded", file)

  # Process DE data
  de_tbl <- data[, .(gene_name, unstranded)]
  de_tbl[, gene_name := sub("\\..*", "", gene_name)]
  de_tbl <- de_tbl[
    , .(unstranded = as.integer(round(sum(unstranded)))), by = gene_name
  ]
  setnames(de_tbl, "unstranded", file)

  # Store data.tables in lists
  infil_list[[file]] <- infil_tbl
  de_list[[file]] <- de_tbl
}

print("Merging Data...")

# Use Reduce to merge data.tables
merged_de_tbl <- Reduce(
  function(x, y) merge(x, y, by = "gene_name", all = TRUE), de_list
)
merged_infil_tbl <- Reduce(
  function(x, y) merge(x, y, by = "gene_name", all = TRUE), infil_list
)

# Preserve gene names to use as row names in matrices
gene_names <- merged_de_tbl$gene_name

# Remove gene_name column to prep data for conversion
merged_de_tbl <- merged_de_tbl[, -1, with = FALSE]
merged_infil_tbl <- merged_infil_tbl[, -1, with = FALSE]

# Convert data.tables into numeric matrices
merged_de_matrix <- as.matrix(merged_de_tbl)
merged_infil_matrix <- as.matrix(merged_infil_tbl)

# Preserve gene_names as row names in matrices
rownames(merged_de_matrix) <- gene_names
rownames(merged_infil_matrix) <- gene_names

# Replace filenames with metadata Case_ID
de_colnames <- metadata[
  match(colnames(merged_de_matrix), metadata$"File Name"), "Case ID"
]
infil_colnames <- metadata[
  match(colnames(merged_infil_matrix), metadata$"File Name"), "Case ID"
]
colnames(merged_de_matrix) <- de_colnames$"Case ID"
colnames(merged_infil_matrix) <- infil_colnames$"Case ID"

print("Writing Merged Data...")

# Label file paths for output
output_file1 <- file.path("results", "merged_tpm_counts.csv")
output_file2 <- file.path("results", "merged_unstranded_counts.csv")

# Write processed data for output
write.csv(merged_de_matrix, output_file1, row.names = TRUE)
write.csv(merged_infil_matrix, output_file2, row.names = TRUE)

print("Processing Complete!")