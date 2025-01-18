# Need to install the readr package to read RNA Seq TSV files
install.packages("readr")

# Load installed readr package
library(readr)

# Make a list of files to access for iterating
files <- list.files(path = "test/", pattern = "*.tsv", full.names = TRUE)

# Install dplyr package to help manipulate data structures
install.packages("dplyr")

# Load installed dplyr package
library(dyplr)

# Make an empty dataframe for merging
merged_df <- data.frame()

# Read every file into a data frame
for (file in files) {
  df <- read_tsv(file)

  # Keep only gene_name and tpm_unstranded columns to use EPIC on
  df1 <- df[, c("gene_name", "tpm_unstranded")]
  print(df1)
  merged_df <- merge(df1, merged_df, by = "gene_name", all = TRUE)
}
print(merged_df)