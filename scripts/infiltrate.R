# Install EPIC


# Load in merged_data
merged_path <- file.path("results", "merged_gene_counts.csv")
data <- read.csv(merged_path)
print(data)