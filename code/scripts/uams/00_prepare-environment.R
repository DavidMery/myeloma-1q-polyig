################################################################################
# STEP 00 — PREPARE / VERIFY R ENVIRONMENT
#
# Purpose:
#   Verify that all CRAN and Bioconductor packages required by the complete
#   reproducible analysis pipeline are available.
#
# Runtime behavior:
#
#   ENVIRONMENT_MODE = "LOCAL"
#     - Intended for local VS Code / RStudio use.
#     - Missing CRAN packages are installed automatically.
#     - Missing Bioconductor packages are installed automatically.
#
#   ENVIRONMENT_MODE = "CODE_OCEAN"
#     - Intended for Code Ocean.
#     - This script does NOT install or update packages.
#     - The pipeline stops if a required package is unavailable.
#
# Important:
#   ENVIRONMENT_MODE is defined in:
#
#     code/run_pipeline.R
#
#   Do not define a second environment toggle in this script.
################################################################################


################################################################################
# 1. VERIFY ENVIRONMENT MODE
################################################################################

if (!exists(
  "ENVIRONMENT_MODE"
)) {
  stop(
    "\nENVIRONMENT_MODE is not defined.\n",
    "Run this script through code/run_pipeline.R.\n"
  )
}


ENVIRONMENT_MODE <- toupper(
  trimws(
    ENVIRONMENT_MODE
  )
)


if (!ENVIRONMENT_MODE %in% c(
  "LOCAL",
  "CODE_OCEAN"
)) {
  stop(
    "ENVIRONMENT_MODE must be either 'LOCAL' or 'CODE_OCEAN'."
  )
}


################################################################################
# 2. REQUIRED PACKAGES
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


ALL_PACKAGES <- unique(
  c(
    CRAN_PACKAGES,
    BIOC_PACKAGES
  )
)


################################################################################
# 3. HELPER — CHECK PACKAGE AVAILABILITY
################################################################################

check_package_availability <- function(
    packages
) {
  
  vapply(
    packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
}


################################################################################
# 4. CHECK CURRENT PACKAGE AVAILABILITY
################################################################################

package_available <- check_package_availability(
  ALL_PACKAGES
)


missing_packages <- ALL_PACKAGES[
  !package_available
]


################################################################################
# 5. REPORT INITIAL STATUS
################################################################################

cat(
  "\n",
  "============================================================\n",
  "STEP 00 — R ENVIRONMENT CHECK\n",
  "============================================================\n",
  "\n",
  "Environment mode: ",
  ENVIRONMENT_MODE,
  "\n",
  "R version: ",
  R.version.string,
  "\n",
  sep = ""
)


if (length(
  missing_packages
) == 0) {
  
  cat(
    "\nAll required packages are already available.\n"
  )
  
} else {
  
  cat(
    "\nMissing package(s):\n",
    paste(
      missing_packages,
      collapse = ", "
    ),
    "\n",
    sep = ""
  )
}


################################################################################
# 6. CODE OCEAN MODE
#
# Do NOT install or update packages.
################################################################################

if (
  ENVIRONMENT_MODE == "CODE_OCEAN" &&
  length(missing_packages) > 0
) {
  
  stop(
    "\nRequired package(s) are missing from the Code Ocean environment:\n",
    paste(
      missing_packages,
      collapse = ", "
    ),
    "\n\n",
    "Add the missing package(s) to the Code Ocean environment configuration ",
    "before running the pipeline.\n"
  )
}


################################################################################
# 7. LOCAL MODE — INSTALL MISSING PACKAGES
################################################################################

if (
  ENVIRONMENT_MODE == "LOCAL" &&
  length(missing_packages) > 0
) {
  
  cat(
    "\n",
    "============================================================\n",
    "LOCAL MODE — INSTALLING MISSING PACKAGES\n",
    "============================================================\n",
    sep = ""
  )
  
  
  ##############################################################################
  # 7A. INSTALL BiocManager FIRST IF NEEDED
  ##############################################################################
  
  if (!requireNamespace(
    "BiocManager",
    quietly = TRUE
  )) {
    
    cat(
      "\nInstalling BiocManager...\n"
    )
    
    install.packages(
      "BiocManager",
      repos = "https://cloud.r-project.org"
    )
  }
  
  
  if (!requireNamespace(
    "BiocManager",
    quietly = TRUE
  )) {
    
    stop(
      "\nBiocManager could not be installed.\n",
      "Bioconductor packages cannot be installed without BiocManager.\n"
    )
  }
  
  
  ##############################################################################
  # 7B. INSTALL MISSING CRAN PACKAGES
  ##############################################################################
  
  current_package_available <- check_package_availability(
    ALL_PACKAGES
  )
  
  
  missing_cran <- CRAN_PACKAGES[
    !current_package_available[
      match(
        CRAN_PACKAGES,
        ALL_PACKAGES
      )
    ]
  ]
  
  
  missing_cran <- setdiff(
    missing_cran,
    "BiocManager"
  )
  
  
  if (length(
    missing_cran
  ) > 0) {
    
    cat(
      "\nInstalling missing CRAN packages:\n",
      paste(
        missing_cran,
        collapse = ", "
      ),
      "\n\n",
      sep = ""
    )
    
    
    install.packages(
      missing_cran,
      repos = "https://cloud.r-project.org"
    )
  }
  
  
  ##############################################################################
  # 7C. INSTALL MISSING BIOCONDUCTOR PACKAGES
  ##############################################################################
  
  current_package_available <- check_package_availability(
    ALL_PACKAGES
  )
  
  
  missing_bioc <- BIOC_PACKAGES[
    !current_package_available[
      match(
        BIOC_PACKAGES,
        ALL_PACKAGES
      )
    ]
  ]
  
  
  if (length(
    missing_bioc
  ) > 0) {
    
    cat(
      "\nInstalling missing Bioconductor packages:\n",
      paste(
        missing_bioc,
        collapse = ", "
      ),
      "\n\n",
      sep = ""
    )
    
    
    BiocManager::install(
      missing_bioc,
      ask = FALSE,
      update = FALSE
    )
  }
}


################################################################################
# 8. VERIFY AGAIN AFTER OPTIONAL INSTALLATION
################################################################################

package_available <- check_package_availability(
  ALL_PACKAGES
)


missing_packages <- ALL_PACKAGES[
  !package_available
]


if (length(
  missing_packages
) > 0) {
  
  if (
    ENVIRONMENT_MODE == "LOCAL"
  ) {
    
    stop(
      "\nRequired package(s) are still missing from the local R environment:\n",
      paste(
        missing_packages,
        collapse = ", "
      ),
      "\n\n",
      "One or more package installations failed.\n",
      "Review the installation messages above before running the pipeline again.\n"
    )
    
  } else {
    
    stop(
      "\nRequired package(s) are missing from the Code Ocean environment:\n",
      paste(
        missing_packages,
        collapse = ", "
      ),
      "\n"
    )
  }
}


################################################################################
# 9. REPORT PACKAGE VERSIONS
################################################################################

package_versions <- vapply(
  ALL_PACKAGES,
  function(package_name) {
    
    as.character(
      utils::packageVersion(
        package_name
      )
    )
    
  },
  FUN.VALUE = character(1)
)


################################################################################
# 10. FINAL ENVIRONMENT REPORT
################################################################################

cat(
  "\n",
  "============================================================\n",
  "R ENVIRONMENT VERIFIED\n",
  "============================================================\n",
  "\n",
  "Environment mode:\n",
  ENVIRONMENT_MODE,
  "\n\n",
  "R version:\n",
  R.version.string,
  "\n\n",
  "CRAN packages available:\n",
  paste(
    CRAN_PACKAGES,
    collapse = ", "
  ),
  "\n\n",
  "Bioconductor packages available:\n",
  paste(
    BIOC_PACKAGES,
    collapse = ", "
  ),
  "\n\n",
  "Package versions:\n",
  sep = ""
)


for (package_name in ALL_PACKAGES) {
  
  cat(
    "  ",
    package_name,
    ": ",
    package_versions[
      package_name
    ],
    "\n",
    sep = ""
  )
}


cat(
  "\n",
  "============================================================\n",
  "STEP 00 COMPLETE\n",
  "============================================================\n",
  sep = ""
)