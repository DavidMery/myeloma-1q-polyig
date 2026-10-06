################################################################################
# PREPARE PROJECT ENVIRONMENT
#
# Installs and loads packages required for the analysis.
#
# Run from the project root:
#   source("scripts/00_prepare-environment.R")
################################################################################

# CRAN
install.packages("data.table")

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("Biobase", ask = FALSE, update = FALSE)

# Load packages
library(data.table)
library(Biobase)

