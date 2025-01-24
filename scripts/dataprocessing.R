# Install the data.table package to handle large data
install.packages("data.table")

# Load installed package
library(data.table)

# Make a list of files to access for iterating
files <- list.files(path = "data/", pattern = "*.tsv", full.names = TRUE)

# Make empty lists to collect data.tables for analyses
de_list <- list()
infil_list <- list()

print("Processing TSV Files...")

# Read every file and process columns for analyses
for (file in files) {

  # Take only sample name from filename
  sample_name <- sub("\\..*", "", basename(file))

  # Read file as a data.table
  data <- fread(file)

  # Only keep rows that have protein_coding genes
  data <- data[gene_type == "protein_coding"]

  # Process infiltrate data
  infil_tbl <- data[, .(gene_name, tpm_unstranded)]

  # Clean the gene names by removing transcript versions
  infil_tbl[, gene_name := sub("\\..*", "", gene_name)]

  # Remove duplicate genes
  infil_tbl <- infil_tbl[!duplicated(gene_name)]

  # Set count column to sample name
  setnames(infil_tbl, "tpm_unstranded", sample_name)

  # Process DE data
  de_tbl <- data[, .(gene_name, unstranded)]

  # Clean the gene names by removing transcript versions
  de_tbl[, gene_name := sub("\\..*", "", gene_name)]

  # Remove duplicate genes
  de_tbl <- de_tbl[!duplicated(gene_name)]

  # Set count column to sample name
  setnames(de_tbl, "unstranded", sample_name)

  # Store data.tables in lists
  infil_list[[sample_name]] <- infil_tbl
  de_list[[sample_name]] <- de_tbl
}

# Sort the data.tables by gene_name to speed up merging
setorder(de_list[[1]], gene_name)
setorder(infil_list[[1]], gene_name)

print("Merging Data...")

# Use Reduce to merge data.tables efficiently
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

print("Writing Merged Data...")

# Label file paths for output
output_file1 <- file.path("results", "merged_tpm_counts.csv")
output_file2 <- file.path("results", "merged_unstranded_counts.csv")

# Write processed data for output
write.csv(merged_de_matrix, output_file1, row.names = TRUE)
write.csv(merged_infil_matrix, output_file2, row.names = TRUE)

print("Processing Complete!")