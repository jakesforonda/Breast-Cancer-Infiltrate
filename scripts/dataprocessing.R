# Need to install the readr package to access RNA Seq TSV files
install.packages("readr")

# Load installed readr package
library(readr)

# Use file indexing to better access files
files <- list.files(path = "test/", pattern = "*.tsv", full.names = TRUE)

# Read file into data frame

df1 <- read_tsv(files[1])
print(df1)
