################################################################################
# MASTER R PIPELINE
#
# Purpose:
#   Run the complete reproducible myeloma analysis pipeline in a single
#   R session.
#
#   The master pipeline controls:
#     1. Whether EXAMPLE or protected UAMS data are used.
#     2. Whether packages may be installed locally or only verified in
#        Code Ocean.
#     3. Where data, models, and results are read/written.
#     4. Execution of Steps 00–33, including Step 25b.
#
# IMPORTANT:
#   Run this script from the project root:
#
#     myeloma-1q-polyig/
#
#   Example:
#
#     Rscript code/run_pipeline.R
#
################################################################################


################################################################################
# USER SETTINGS — CHANGE ONLY THIS SECTION
################################################################################

# -------------------------------------------------------------------------------
# PIPELINE_MODE
#
# "EXAMPLE"
#   Input data:
#     GitHub/myeloma-1q-polyig/data/simulated/uams
#
#   Models:
#     GitHub/myeloma-1q-polyig/models/simulated/uams
#
#   Results:
#     GitHub/myeloma-1q-polyig/results/simulated/uams
#
#
# "UAMS"
#   Input data:
#     Box/PROJECTS/myeloma-1q-polyig/data
#
#   Models:
#     Box/PROJECTS/myeloma-1q-polyig/models
#
#   Results:
#     Box/PROJECTS/myeloma-1q-polyig/results
# -------------------------------------------------------------------------------

PIPELINE_MODE <- "EXAMPLE"


# -------------------------------------------------------------------------------
# ENVIRONMENT_MODE
#
# "LOCAL"
#   Intended for local VS Code / RStudio use.
#   Step 00 installs packages that are missing.
#
# "CODE_OCEAN"
#   Intended for Code Ocean.
#   Step 00 does NOT install or update packages.
#   It only verifies that the capsule environment is complete.
# -------------------------------------------------------------------------------

ENVIRONMENT_MODE <- "LOCAL"


################################################################################
# 1. VERIFY PROJECT ROOT
################################################################################

if (
  !dir.exists("code") ||
  !dir.exists(
    file.path(
      "code",
      "scripts",
      "uams"
    )
  )
) {
  stop(
    "\nPipeline must be run from the project root.\n",
    "\nExpected to find:\n",
    "  code/\n",
    "  code/scripts/uams/\n",
    "\nCurrent working directory:\n",
    normalizePath(
      ".",
      winslash = "/",
      mustWork = TRUE
    ),
    "\n"
  )
}


################################################################################
# 2. DEFINE PROJECT ROOT
################################################################################

PROJECT_ROOT <- normalizePath(
  ".",
  winslash = "/",
  mustWork = TRUE
)

SCRIPT_ROOT <- file.path(
  PROJECT_ROOT,
  "code",
  "scripts",
  "uams"
)


################################################################################
# 3. VALIDATE USER SETTINGS
################################################################################

PIPELINE_MODE <- toupper(
  trimws(
    PIPELINE_MODE
  )
)

ENVIRONMENT_MODE <- toupper(
  trimws(
    ENVIRONMENT_MODE
  )
)


if (!PIPELINE_MODE %in% c(
  "EXAMPLE",
  "UAMS"
)) {
  stop(
    "PIPELINE_MODE must be either 'EXAMPLE' or 'UAMS'."
  )
}


if (!ENVIRONMENT_MODE %in% c(
  "LOCAL",
  "CODE_OCEAN"
)) {
  stop(
    "ENVIRONMENT_MODE must be either 'LOCAL' or 'CODE_OCEAN'."
  )
}


################################################################################
# 4. DEFINE DATA / MODEL / RESULTS ROOTS
################################################################################

if (PIPELINE_MODE == "EXAMPLE") {
  
  DATA_ROOT <- file.path(
    PROJECT_ROOT,
    "data",
    "simulated",
    "uams"
  )
  
  MODEL_ROOT <- file.path(
    PROJECT_ROOT,
    "models",
    "simulated",
    "uams"
  )
  
  RESULTS_ROOT <- file.path(
    PROJECT_ROOT,
    "results",
    "simulated",
    "uams"
  )
  
} else if (PIPELINE_MODE == "UAMS") {
  
  USER_PROFILE <- Sys.getenv(
    "USERPROFILE"
  )
  
  if (
    is.na(USER_PROFILE) ||
    USER_PROFILE == ""
  ) {
    stop(
      "USERPROFILE environment variable is not available. ",
      "Cannot construct the protected Box project path."
    )
  }
  
  UAMS_ROOT <- file.path(
    USER_PROFILE,
    "Box",
    "PROJECTS",
    "myeloma-1q-polyig"
  )
  
  DATA_ROOT <- file.path(
    UAMS_ROOT,
    "data"
  )
  
  MODEL_ROOT <- file.path(
    UAMS_ROOT,
    "models"
  )
  
  RESULTS_ROOT <- file.path(
    UAMS_ROOT,
    "results"
  )
}


################################################################################
# 5. VERIFY / CREATE ROOT DIRECTORIES
################################################################################

if (!dir.exists(DATA_ROOT)) {
  stop(
    "\nDATA_ROOT does not exist:\n",
    DATA_ROOT,
    "\n"
  )
}


if (!dir.exists(MODEL_ROOT)) {
  
  created <- dir.create(
    MODEL_ROOT,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  if (
    !created &&
    !dir.exists(MODEL_ROOT)
  ) {
    stop(
      "\nCould not create MODEL_ROOT:\n",
      MODEL_ROOT,
      "\n"
    )
  }
}


if (!dir.exists(RESULTS_ROOT)) {
  
  created <- dir.create(
    RESULTS_ROOT,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  if (
    !created &&
    !dir.exists(RESULTS_ROOT)
  ) {
    stop(
      "\nCould not create RESULTS_ROOT:\n",
      RESULTS_ROOT,
      "\n"
    )
  }
}


################################################################################
# 6. REPORT ACTIVE CONFIGURATION
################################################################################

cat(
  "\n",
  "######################################################################\n",
  "# ACTIVE PIPELINE CONFIGURATION\n",
  "######################################################################\n",
  "#\n",
  "# PIPELINE MODE:\n",
  "#   ",
  PIPELINE_MODE,
  "\n",
  "#\n",
  "# ENVIRONMENT MODE:\n",
  "#   ",
  ENVIRONMENT_MODE,
  "\n",
  "#\n",
  "# PROJECT ROOT:\n",
  "#   ",
  PROJECT_ROOT,
  "\n",
  "#\n",
  "# CODE:\n",
  "#   ",
  SCRIPT_ROOT,
  "\n",
  "#\n",
  "# INPUT DATA:\n",
  "#   ",
  DATA_ROOT,
  "\n",
  "#\n",
  "# MODELS:\n",
  "#   ",
  MODEL_ROOT,
  "\n",
  "#\n",
  "# RESULTS:\n",
  "#   ",
  RESULTS_ROOT,
  "\n",
  "#\n",
  "######################################################################\n",
  "\n",
  sep = ""
)


################################################################################
# 7. PREPARE / VERIFY R ENVIRONMENT
################################################################################

ENVIRONMENT_SCRIPT <- file.path(
  SCRIPT_ROOT,
  "00_prepare-environment.R"
)


if (!file.exists(ENVIRONMENT_SCRIPT)) {
  stop(
    "Environment script not found:\n",
    ENVIRONMENT_SCRIPT
  )
}


source(
  ENVIRONMENT_SCRIPT,
  local = FALSE
)


################################################################################
# 8. LOAD REQUIRED PACKAGES ONCE
################################################################################

suppressPackageStartupMessages({
  
  library(Biobase)
  library(caret)
  library(circlize)
  library(ComplexHeatmap)
  library(ConsensusClusterPlus)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(ggrepel)
  library(glmnet)
  library(gridExtra)
  library(limma)
  library(lubridate)
  library(matrixStats)
  library(readr)
  library(survival)
  library(survminer)
  library(tibble)
  
})


################################################################################
# 9. PIPELINE STEPS
################################################################################

PIPELINE_SCRIPTS <- c(
  
  "01_check-missing-values.R",
  
  "02_prepare-uams-data.R",
  
  "03_train-gep1q-elastic-net-cv.R",
  
  "04_evaluate-gep1q-penalized-heldout-confusion-matrix.R",
  
  "05_apply-gep1q-lasso-score.R",
  
  "06_add-fish-data.R",
  
  "07_train-and-apply-polyig-score.R",
  
  "08_classify-pc-me-genes.R",
  
  "09_fig1a-gep70-km.R",
  
  "10_fig1b-cox-volcano.R",
  
  "11_fig1c-consensus-clustering-heatmap.R",
  
  "12_fig1d_cluster_enrichment.R",
  
  "13_fish1q-ttp-km.R",
  
  "14_unigh-ttp-km.R",
  
  "15_fish1q-unigh-ttp-km.R",
  
  "16_gep1q-ttp-km.R",
  
  "17_polyig-ttp-km.R",
  
  "18_fig1e-gep1q-polyig-km.R",
  
  "19_add-uams-binary-pdata-columns.R",
  
  "20_add-uams-composite-rnas-groups.R",
  
  "21_match-rnas-groups-to-rnbx-by-patid.R",
  
  "22_rnbx-mgus-vs-gep70-limma.R",
  
  "23_fig2b_Microenvironment_Heatmap.R",
  
  "24_fig2a_Microenvironment_Heatmap.R",
  
  "25_fig3a-gep1q-polyig-km-tt2.R",
  
  "25b_compare_TT_Thal_NoThal.R",
  
  "26_fig3b-gep1q-polyig-km-tt3plus.R",
  
  "27_fig4A-D_4_6_8yr_full_followup.R",
  
  "28_fig5a-gep1q-polyig-km.R",
  
  "29_fig5b-fish1q-unigh-ttp-km.R",
  
  "30_overlay-gep-and-clinical-ttr-curves.R",
  
  "31_clinical-gep-confusion-matrix.R",
  
  "32_univariate-cox.R",
  
  "33_multivariate-cox.R"
)


################################################################################
# 10. VERIFY ALL PIPELINE SCRIPTS EXIST
################################################################################

PIPELINE_PATHS <- file.path(
  SCRIPT_ROOT,
  PIPELINE_SCRIPTS
)


missing_scripts <- PIPELINE_SCRIPTS[
  !file.exists(
    PIPELINE_PATHS
  )
]


if (length(missing_scripts) > 0) {
  
  stop(
    "\nThe following pipeline script(s) are missing:\n",
    paste(
      missing_scripts,
      collapse = "\n"
    ),
    "\n"
  )
}


################################################################################
# 11. RUN PIPELINE
################################################################################

PIPELINE_START_TIME <- Sys.time()


for (i in seq_along(
  PIPELINE_SCRIPTS
)) {
  
  script <- PIPELINE_SCRIPTS[i]
  
  script_path <- file.path(
    SCRIPT_ROOT,
    script
  )
  
  
  cat(
    "\n\n",
    "============================================================\n",
    "RUNNING STEP ",
    i,
    " OF ",
    length(PIPELINE_SCRIPTS),
    "\n",
    script,
    "\n",
    "============================================================\n",
    sep = ""
  )
  
  
  STEP_START_TIME <- Sys.time()
  
  
  source(
    script_path,
    local = FALSE
  )
  
  
  STEP_END_TIME <- Sys.time()
  
  STEP_DURATION <- difftime(
    STEP_END_TIME,
    STEP_START_TIME,
    units = "mins"
  )
  
  
  cat(
    "\n",
    "COMPLETED: ",
    script,
    "\n",
    "Elapsed time: ",
    round(
      as.numeric(
        STEP_DURATION
      ),
      2
    ),
    " minutes\n",
    sep = ""
  )
}


################################################################################
# 12. PIPELINE COMPLETE
################################################################################

PIPELINE_END_TIME <- Sys.time()

PIPELINE_DURATION <- difftime(
  PIPELINE_END_TIME,
  PIPELINE_START_TIME,
  units = "mins"
)


cat(
  "\n\n",
  "============================================================\n",
  "COMPLETE PIPELINE FINISHED SUCCESSFULLY\n",
  "============================================================\n",
  "\n",
  "Pipeline mode: ",
  PIPELINE_MODE,
  "\n",
  "Environment mode: ",
  ENVIRONMENT_MODE,
  "\n",
  "Total pipeline time: ",
  round(
    as.numeric(
      PIPELINE_DURATION
    ),
    2
  ),
  " minutes\n",
  "\n",
  "Results directory:\n",
  RESULTS_ROOT,
  "\n",
  "============================================================\n",
  sep = ""
)