################################################################################
# STEP 00 — VERIFY R ENVIRONMENT
#
# Purpose:
#   Verify that all CRAN and Bioconductor packages required by the complete
#   reproducible analysis pipeline are available in the runtime environment.
#
# Notes:
#   - Packages are installed through the Code Ocean environment configuration.
#   - This script does not install or update packages.
#   - The pipeline stops immediately if a required package is unavailable.
################################################################################


################################################################################
# 1. REQUIRED PACKAGES
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

BIOC_PACKAGES <- c(
  "Biobase",
  "ComplexHeatmap",
  "ConsensusClusterPlus",
  "limma"
)

ALL_PACKAGES <- c(
  CRAN_PACKAGES,
  BIOC_PACKAGES
)


################################################################################
# 2. CHECK PACKAGE AVAILABILITY
################################################################################

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
    "\nRequired package(s) are missing from the R environment:\n",
    paste(missing_packages, collapse = ", "),
    "\n\n",
    "Install the missing package(s) in the Code Ocean Environment configuration.\n"
  )
}


################################################################################
# 3. REPORT ENVIRONMENT
################################################################################

cat(
  "\n",
  "============================================================\n",
  "R ENVIRONMENT VERIFIED\n",
  "============================================================\n",
  "\n",
  "R version:\n",
  R.version.string,
  "\n\n",
  "CRAN packages available:\n",
  paste(CRAN_PACKAGES, collapse = ", "),
  "\n\n",
  "Bioconductor packages available:\n",
  paste(BIOC_PACKAGES, collapse = ", "),
  "\n",
  "============================================================\n"
)