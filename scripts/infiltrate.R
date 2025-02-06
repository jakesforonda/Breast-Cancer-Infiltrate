# Load necessary packages
library(pheatmap)
library(MCPcounter)
library(httpgd)

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

# Save the heatmap as PNG
output_file <- file.path("results", "infil_heatmap.png")
png(output_file, width = 1200, height = 800)
pheatmap(as.matrix(results_log),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "MCPcounter Cell Type Estimates (Log Transformed)",
  display_numbers = FALSE,
  show_colnames = FALSE
)
dev.off()

# Start httpgd graphics device
hgd()
hgd_view()

# Plot the heatmap with customizations
pheatmap(as.matrix(results_log),
  cluster_rows = TRUE,   # Cluster Rows
  cluster_cols = TRUE,   # Cluster Columns
  main = "MCPcounter Cell Type Estimates (Log Transformed)", # Title
  display_numbers = FALSE, # Remove numbers to reduce clutter
  show_colnames = FALSE  # Remove sample names to reduce clutter
)

print("infiltrate.R complete")
