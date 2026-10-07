################################################################################
# MASTER R PIPELINE
#
# Runs the complete synthetic reproducibility workflow in a single R session.
################################################################################


################################################################################
# 1. VERIFY PROJECT ROOT
################################################################################

if (!dir.exists("data") ||
    !dir.exists("code") ||
    !dir.exists(file.path("code", "scripts"))) {

  stop(
    "Pipeline must be run from the capsule/project root."
  )
}


################################################################################
# 2. PREPARE ENVIRONMENT
################################################################################

source(
  file.path(
    "code",
    "scripts",
    "00_prepare-environment.R"
  )
)


################################################################################
# 3. LOAD REQUIRED PACKAGES ONCE
################################################################################

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


################################################################################
# 4. PIPELINE STEPS
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
  "26_fig3b-gep1q-polyig-km-tt3plus.R",
  "27_Fig4A-D_4_6_8yr_Full_followup.R",
  "28_fig5a-gep1q-polyig-km.R",
  "29_fig5b-fish1q-unigh-ttp-km.R",
  "30_overlay-gep-and-clinical-ttr-curves.R",
  "31_clinical-gep-confusion-matrix.R",
  "32_univariate-cox.R",
  "33_multivariate-cox.R"
)


################################################################################
# 5. RUN PIPELINE
################################################################################

for (script in PIPELINE_SCRIPTS) {

  script_path <- file.path(
    "code",
    "scripts",
    script
  )

  cat(
    "\n\n",
    "============================================================\n",
    "RUNNING: ",
    script,
    "\n",
    "============================================================\n",
    sep = ""
  )

  if (!file.exists(script_path)) {

    stop(
      "Pipeline script not found: ",
      script_path
    )
  }

  source(
    script_path,
    local = FALSE
  )

  cat(
    "\nCOMPLETED: ",
    script,
    "\n",
    sep = ""
  )
}


################################################################################
# 6. COMPLETE
################################################################################

cat(
  "\n\n",
  "============================================================\n",
  "COMPLETE PIPELINE FINISHED SUCCESSFULLY\n",
  "============================================================\n"
)