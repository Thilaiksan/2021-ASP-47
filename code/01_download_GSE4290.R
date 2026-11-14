
#################################################
# Project: Glioblastoma Bioinformatics Analysis
# Dataset: GSE4290
# Source: NCBI Gene Expression Omnibus (GEO)
#################################################

# Install GEOquery if necessary
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("GEOquery", quietly = TRUE)) {
  BiocManager::install("GEOquery")
}

library(GEOquery)

# Set project directory
project_path <- "C:/Users/hp/Desktop/Glioblastoma_Bioinformatics_Project"
setwd(project_path)

# Create data folder if necessary
dir.create("data", showWarnings = FALSE)

# Download GSE4290
gse_list <- getGEO(
  "GSE4290",
  GSEMatrix = TRUE
)

# Save downloaded GEO object
saveRDS(
  gse_list,
  "data/GSE4290_GEO_Data.rds"
)

# Display dataset information
cat("Dataset downloaded successfully!\n")
cat("Number of platform objects:", length(gse_list), "\n")

