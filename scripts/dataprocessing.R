# Need to install the readr package to read RNA Seq TSV files
install.packages("readr")

# Load installed readr package
library(readr)

# Make a list of files to access for iterating
files <- list.files(path = "test/", pattern = "*.tsv", full.names = TRUE)

# Install packages to help manipulate data structures
install.packages("dplyr")
install.packages("tidyr")
install.packages("purrr")

# Load installed packages
library(dplyr)
library(tidyr)
library(purrr)

# Make an empty list to collect data frames
df_list <- list()

# Read every file into a data frame
for (file in files) {

  # Take only sample name from filename
  sample_name <- sub("\\..*", "", basename(file))

  # Read data
  data <- read_tsv(file)

  # Keep only gene_name and tpm_unstranded columns to use EPIC on
  df <- data %>%
    select(gene_name, tpm_unstranded) %>%
    rename(!!sample_name := tpm_unstranded)

  # Store in list
  df_list[[sample_name]] <- df
}

# Merge the collected data frames of each sample
merged_df <- reduce(df_list, full_join, by = "gene_name")

# Save the merged data into results folder
output_folder <- "results"
output_file <- file.path(output_folder, "merged_gene_counts.csv")
write.csv(merged_df, output_file, row.names = FALSE)
