# Need to install the readr package to read RNA Seq TSV files
install.packages("readr")

# Install packages to help manipulate data structures
install.packages("dplyr")
install.packages("tidyr")
install.packages("purrr")

# Load installed packages
library(dplyr)
library(tidyr)
library(purrr)
library(readr)

# Make a list of files to access for iterating
files <- list.files(path = "test/", pattern = "*.tsv", full.names = TRUE)

# Make empty lists to collect tibbles for analyses
de_list <- list()
infil_list <- list()

# Read every file and process columns for analyses
for (file in files) {

  # Take only sample name from filename
  sample_name <- sub("\\..*", "", basename(file))

  # Read file as a tibble
  data <- read_tsv(file)

  # Keep only gene_name and tpm_unstranded columns for infiltrate
  infil_tbl <- data %>%

    # Remove duplicate values in gene_names
    distinct(gene_name, .keep_all = TRUE) %>%

    # Collect gene_name and tpm_usntranded columns
    select(gene_name, tpm_unstranded) %>%

    # Rename the column name as the sample name
    rename(!!sample_name := tpm_unstranded)

  # Keep only gene_name and tpm_unstranded columns for DE
  de_tbl <- data %>%

    # Remove duplicate values in gene_names
    distinct(gene_name, .keep_all = TRUE) %>%

    # Collect gene_name and usntranded columns
    select(gene_name, unstranded) %>%

    # Rename the column name as the sample name
    rename(!!sample_name := unstranded)

  # Store tibbles in list
  infil_list[[sample_name]] <- infil_tbl
  de_list[[sample_name]] <- de_tbl
}

# Merge the collected tibbles of each sample
merged_de_tbl <- reduce(de_list, full_join, by = "gene_name")
merged_infil_tbl <- reduce(infil_list, full_join, by = "gene_name")

# Preserve gene names to use as row names in matrices
gene_names <- merged_de_tbl$gene_name

# Remove gene_name column to prep data for conversion
merged_de_tbl <- merged_de_tbl[, -1]
merged_infil_tbl <- merged_infil_tbl[, -1]

# Convert tibbles into numeric matrices
merged_de_matrix <- apply(as.matrix(merged_de_tbl), 2, as.numeric)
merged_infil_matrix <- apply(as.matrix(merged_infil_tbl), 2, as.numeric)

# Preserve gene_names as row names in matrices
rownames(merged_de_matrix) <- gene_names
rownames(merged_infil_matrix) <- gene_names

# Save the merged data into results folder

# Label file paths for output
output_file1 <- file.path("results", "merged_tpm_counts.csv")
output_file2 <- file.path("results", "merged_unstranded_counts.csv")

# Write processed data for output
write.csv(merged_de_matrix, output_file1, row.names = TRUE)
write.csv(merged_infil_matrix, output_file2, row.names = TRUE)