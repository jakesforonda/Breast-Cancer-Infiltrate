# Load installed pheatmap package
library(pheatmap)

# Load MCPCounter
library(MCPcounter)

print("Running infiltrate.R")
# Assign file path to variable
data_path <- file.path("results", "merged_tpm_counts.csv")

# Read merged TPM_unstranded counts
merged_data <- read.csv(
  data_path, header = TRUE, row.names = 1, check.names = FALSE
)

# Input the matrix with MCPcounter
results <- MCPcounter.estimate(merged_data, featuresType = "HUGO_symbols")

# Log-transform the results (add 1 to avoid log(0))
results_log <- log2(results + 1)

# Output the log-transformed results
output_file <- file.path("results", "infil_counts.csv")
write.csv(results_log, output_file, row.names = TRUE)

# Generate png for heatmap
output_file <- file.path("results", "infil_heatmap.png")
png(output_file, width = 1200, height = 800)

# Plot the heatmap with customizations
pheatmap(as.matrix(results_log),
  # Cluster Rows
  cluster_rows = TRUE,
  # Cluster Columns
  cluster_cols = TRUE,
  # Title the graph
  main = "MCPcounter Cell Type Estimates (Log Transformed)",
  # Remove numbers to reduce clutter
  display_numbers = FALSE,
  # Remove sample_names to reduce clutter
  show_colnames = FALSE
)
print("infiltrate.R complete")
# Close png file
dev.off()