# Install MCPcounter for immune deconvolution per GitHub instructions
install.packages(c("devtools", "curl"))
library(devtools)
install_github("ebecht/MCPcounter", ref = "master", subdir = "Source")

# Load installed MCPcounter package
library(MCPcounter)

# Load the merged data
data_path <- file.path("results", "merged_tpm_counts.csv")
merged_data <- read.csv(
  data_path, header = TRUE, row.names = 1, check.names = FALSE
)

# Input the matrix with MCPcounter
results <- MCPcounter.estimate(merged_data, featuresType = "HUGO_symbols")

# Output the results
output_file <- file.path("results", "infil_counts.csv")
write.csv(results, output_file, row.names = TRUE)

# Install pheatmap to being analyzing results
install.packages("pheatmap")

# Load installed pheatmap package
library(pheatmap)

# Generate png for heatmap
output_file <- file.path("results", "infil_heatmap.png")
png(output_file, width = 800, height = 600)

# Plot the heatmap
pheatmap(as.matrix(results), cluster_rows = TRUE, cluster_cols = TRUE,
         main = "MCPcounter Cell Type Estimates")

# Close png file
dev.off()
