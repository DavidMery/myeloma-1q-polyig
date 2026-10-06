################################################################################
# STEP 00 — PREPARE R ENVIRONMENT
#
# Purpose:
#   Install all CRAN and Bioconductor packages required by the complete
#   reproducible analysis pipeline.
#
# This script:
#   1. Checks whether each required package is already installed.
#   2. Installs only packages that are missing.
#   3. Installs CRAN packages from CRAN.
#   4. Installs Bioconductor packages using BiocManager.
#   5. Verifies that all required packages are available before continuing.
#
# Intended use:
#   Run once before Steps 01–33, or source automatically from the master
#   pipeline.
#
# Run from the project root with:
#
#   source("scripts/00_prepare-environment.R")
#
# Notes:
#   - Base R packages such as stats, utils, grid, methods, parallel, and
#     grDevices do not need to be installed.
#   - This script is compatible with the synthetic reproducibility workflow
#     used in GitHub and Code Ocean.
################################################################################


################################################################################
# 1. CRAN PACKAGES
################################################################################

CRAN_PACKAGES <- c(
  "BiocManager",
  "caret",
  "circlize",
  "data.table",
  "dplyr",
  "ggplot2",
  "ggpubr",
  "ggrepel",
  "glmnet",
  "gridExtra",
  "lubridate",
  "matrixStats",
  "readr",
  "survival",
  "survminer",
  "tibble"
)


################################################################################
# 2. BIOCONDUCTOR PACKAGES
################################################################################

BIOC_PACKAGES <- c(
  "Biobase",
  "ComplexHeatmap",
  "ConsensusClusterPlus",
  "limma"
)


################################################################################
# 3. INSTALL MISSING CRAN PACKAGES
################################################################################

missing_cran <- CRAN_PACKAGES[
  !vapply(
    CRAN_PACKAGES,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_cran) > 0) {
  
  cat(
    "\nInstalling missing CRAN packages:\n",
    paste(missing_cran, collapse = ", "),
    "\n\n"
  )
  
  install.packages(
    missing_cran,
    repos = "https://cloud.r-project.org",
    dependencies = TRUE
  )
  
} else {
  
  cat(
    "\nAll required CRAN packages are already installed.\n"
  )
}


################################################################################
# 4. VERIFY BIOCMANAGER
################################################################################

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  
  install.packages(
    "BiocManager",
    repos = "https://cloud.r-project.org"
  )
}


################################################################################
# 5. INSTALL MISSING BIOCONDUCTOR PACKAGES
################################################################################

missing_bioc <- BIOC_PACKAGES[
  !vapply(
    BIOC_PACKAGES,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_bioc) > 0) {
  
  cat(
    "\nInstalling missing Bioconductor packages:\n",
    paste(missing_bioc, collapse = ", "),
    "\n\n"
  )
  
  BiocManager::install(
    missing_bioc,
    ask = FALSE,
    update = FALSE
  )
  
} else {
  
  cat(
    "\nAll required Bioconductor packages are already installed.\n"
  )
}


################################################################################
# 6. VERIFY ALL REQUIRED PACKAGES
################################################################################

ALL_PACKAGES <- c(
  CRAN_PACKAGES,
  BIOC_PACKAGES
)

package_available <- vapply(
  ALL_PACKAGES,
  requireNamespace,
  quietly = TRUE,
  FUN.VALUE = logical(1)
)

if (!all(package_available)) {
  
  missing_packages <- ALL_PACKAGES[
    !package_available
  ]
  
  stop(
    "\nThe following required package(s) could not be installed:\n",
    paste(missing_packages, collapse = ", "),
    "\n"
  )
}


################################################################################
# 7. REPORT ENVIRONMENT
################################################################################

cat(
  "\n",
  "============================================================\n",
  "R ENVIRONMENT SUCCESSFULLY PREPARED\n",
  "============================================================\n",
  "\n",
  "R version:\n",
  R.version.string,
  "\n\n",
  "CRAN packages:\n",
  paste(CRAN_PACKAGES, collapse = ", "),
  "\n\n",
  "Bioconductor packages:\n",
  paste(BIOC_PACKAGES, collapse = ", "),
  "\n",
  "============================================================\n"
)