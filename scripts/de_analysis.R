# Load installed DESeq2 package
library(DESeq2)

# Read in merged_unstranded_counts
data_path <- file.path("results", "merged_unstranded_counts.csv")
merged_data <- read.csv(
  data_path, header = TRUE, row.names = 1, check.names = FALSE
)
