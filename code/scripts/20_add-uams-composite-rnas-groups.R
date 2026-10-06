################################################################################
# STEP 20B — ADD UAMS COMPOSITE RNAS GROUP COLUMNS
#
# Creates:
#   1. GEP_RNAS_Groups_MGUS_Lasso
#   2. Clinical_RNAS_Groups_MGUS
#
# MGUS is assigned first.
# For NDMM, GEP70 > 0.66 overrides the other variables.
#
# GEP groups use:
#   GEP1q_binary_lasso × PolyIG_binary
#
# Clinical groups use:
#   FISH1q_binary × IG_ARD_binary
#
# The Step 19B ExpressionSet is not overwritten.
################################################################################

suppressPackageStartupMessages({
  library(Biobase)
  library(dplyr)
  library(tibble)
  library(readr)
})


################################################################################
# 1. CONFIGURATION
################################################################################

ESET_PATH <- file.path(
  "data",
  "processed",
  "uams",
  "ESET_uams_all_gep1q_fish_polyig_pcme_binaries_lasso.rds"
)

ESET_OUT_PATH <- file.path(
  "data",
  "processed",
  "uams",
  "ESET_uams_all_gep1q_fish_polyig_pcme_binaries_lasso_groups.rds"
)

QC_DIR <- file.path(
  "results",
  "qc"
)

SUMMARY_FILE <- file.path(
  QC_DIR,
  "uams_composite_rnas_groups_summary_step20b.csv"
)

VALUES_FILE <- file.path(
  QC_DIR,
  "uams_composite_rnas_groups_values_step20b.csv"
)


################################################################################
# 2. DEFINITIONS
################################################################################

NDMM <- "RNAS_CD138_NDMM"
MGUS <- "RNAS_CD138_MGUS"

GEP_LASSO_GROUP <- "GEP_RNAS_Groups_MGUS_Lasso"
CLINICAL_GROUP <- "Clinical_RNAS_Groups_MGUS"

GEP_LEVELS <- c(
  "MGUS",
  "GEP70>0.66",
  "GEP1q-CN2/PolyIG-high",
  "GEP1q-CN2/PolyIG-low",
  "GEP1q-CN>=3/PolyIG-high",
  "GEP1q-CN>=3/PolyIG-low"
)

CLINICAL_LEVELS <- c(
  "MGUS",
  "GEP70>0.66",
  "FISH1q-CN2/IG_ARDhigh",
  "FISH1q-CN2/IG_ARDlow",
  "FISH1q-CN>=3/IG_ARDhigh",
  "FISH1q-CN>=3/IG_ARDlow"
)


################################################################################
# 3. DIRECTORIES + INPUT
################################################################################

dir.create(
  dirname(ESET_OUT_PATH),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  QC_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(ESET_PATH)) {
  stop(
    "Cannot find Step 19B ExpressionSet:\n",
    ESET_PATH
  )
}

if (
  normalizePath(
    ESET_PATH,
    winslash = "/",
    mustWork = FALSE
  ) ==
  normalizePath(
    ESET_OUT_PATH,
    winslash = "/",
    mustWork = FALSE
  )
) {
  stop(
    "ESET_OUT_PATH must differ from ESET_PATH."
  )
}


################################################################################
# 4. LOAD EXPRESSIONSET
################################################################################

eset <- readRDS(
  ESET_PATH
)

if (!methods::is(
  eset,
  "ExpressionSet"
)) {
  stop(
    "Input file is not an ExpressionSet."
  )
}

phenotype <- as.data.frame(
  Biobase::pData(eset),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

phenotype_row_names <- rownames(
  phenotype
)


################################################################################
# 5. REQUIRED COLUMNS
################################################################################

REQUIRED_COLUMNS <- c(
  "sample_group",
  "GEP1q_binary_lasso",
  "PolyIG_binary",
  "GEP70_binary",
  "FISH1q_binary",
  "IG_ARD_binary",
  "MGUS_Binary"
)

missing_columns <- setdiff(
  REQUIRED_COLUMNS,
  colnames(phenotype)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required Step 19B column(s): ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}


################################################################################
# 6. VALIDATE BINARY VARIABLES
################################################################################

BINARY_COLUMNS <- c(
  "GEP1q_binary_lasso",
  "PolyIG_binary",
  "GEP70_binary",
  "FISH1q_binary",
  "IG_ARD_binary",
  "MGUS_Binary"
)

for (column in BINARY_COLUMNS) {
  
  invalid_values <- setdiff(
    unique(
      phenotype[[column]][
        !is.na(
          phenotype[[column]]
        )
      ]
    ),
    c(
      0L,
      1L
    )
  )
  
  if (length(invalid_values) > 0) {
    stop(
      column,
      " contains values other than 0, 1, or NA: ",
      paste(
        invalid_values,
        collapse = ", "
      )
    )
  }
}


################################################################################
# 7. DEFINE SAMPLE GROUPS
################################################################################

is_ndmm <-
  !is.na(phenotype$sample_group) &
  phenotype$sample_group == NDMM

is_mgus <-
  !is.na(phenotype$sample_group) &
  phenotype$sample_group == MGUS


################################################################################
# 8. CHECK MGUS BINARY
################################################################################

if (
  any(
    is_mgus &
    (
      is.na(
        phenotype$MGUS_Binary
      ) |
      phenotype$MGUS_Binary != 1L
    )
  )
) {
  stop(
    "MGUS_Binary is inconsistent with sample_group."
  )
}


################################################################################
# 9. CREATE LASSO GEP GROUP
################################################################################

gep_lasso_group <- rep(
  NA_character_,
  nrow(phenotype)
)

# MGUS first
gep_lasso_group[
  is_mgus
] <- "MGUS"

# NDMM
i <- which(
  is_ndmm
)

gep_lasso_group[i] <- dplyr::case_when(
  
  phenotype$GEP70_binary[i] == 1L ~
    "GEP70>0.66",
  
  phenotype$GEP70_binary[i] == 0L &
    phenotype$GEP1q_binary_lasso[i] == 0L &
    phenotype$PolyIG_binary[i] == 1L ~
    "GEP1q-CN2/PolyIG-high",
  
  phenotype$GEP70_binary[i] == 0L &
    phenotype$GEP1q_binary_lasso[i] == 0L &
    phenotype$PolyIG_binary[i] == 0L ~
    "GEP1q-CN2/PolyIG-low",
  
  phenotype$GEP70_binary[i] == 0L &
    phenotype$GEP1q_binary_lasso[i] == 1L &
    phenotype$PolyIG_binary[i] == 1L ~
    "GEP1q-CN>=3/PolyIG-high",
  
  phenotype$GEP70_binary[i] == 0L &
    phenotype$GEP1q_binary_lasso[i] == 1L &
    phenotype$PolyIG_binary[i] == 0L ~
    "GEP1q-CN>=3/PolyIG-low",
  
  TRUE ~ NA_character_
)


################################################################################
# 10. CREATE CLINICAL / FISH GROUP
################################################################################

clinical_group <- rep(
  NA_character_,
  nrow(phenotype)
)

# MGUS first
clinical_group[
  is_mgus
] <- "MGUS"

# NDMM
clinical_group[i] <- dplyr::case_when(
  
  phenotype$GEP70_binary[i] == 1L ~
    "GEP70>0.66",
  
  phenotype$GEP70_binary[i] == 0L &
    phenotype$FISH1q_binary[i] == 0L &
    phenotype$IG_ARD_binary[i] == 1L ~
    "FISH1q-CN2/IG_ARDhigh",
  
  phenotype$GEP70_binary[i] == 0L &
    phenotype$FISH1q_binary[i] == 0L &
    phenotype$IG_ARD_binary[i] == 0L ~
    "FISH1q-CN2/IG_ARDlow",
  
  phenotype$GEP70_binary[i] == 0L &
    phenotype$FISH1q_binary[i] == 1L &
    phenotype$IG_ARD_binary[i] == 1L ~
    "FISH1q-CN>=3/IG_ARDhigh",
  
  phenotype$GEP70_binary[i] == 0L &
    phenotype$FISH1q_binary[i] == 1L &
    phenotype$IG_ARD_binary[i] == 0L ~
    "FISH1q-CN>=3/IG_ARDlow",
  
  TRUE ~ NA_character_
)


################################################################################
# 11. ADD GROUP COLUMNS
################################################################################

phenotype[[GEP_LASSO_GROUP]] <- factor(
  gep_lasso_group,
  levels = GEP_LEVELS
)

phenotype[[CLINICAL_GROUP]] <- factor(
  clinical_group,
  levels = CLINICAL_LEVELS
)

rownames(phenotype) <-
  phenotype_row_names

Biobase::pData(eset) <-
  phenotype


################################################################################
# 12. VALIDATE GROUP VALUES
################################################################################

check_values <- function(
    x,
    allowed,
    name
) {
  
  observed <- unique(
    as.character(
      x[
        !is.na(x)
      ]
    )
  )
  
  bad <- setdiff(
    observed,
    allowed
  )
  
  if (length(bad) > 0) {
    stop(
      name,
      " contains unexpected values: ",
      paste(
        bad,
        collapse = ", "
      )
    )
  }
}

check_values(
  phenotype[[GEP_LASSO_GROUP]],
  GEP_LEVELS,
  GEP_LASSO_GROUP
)

check_values(
  phenotype[[CLINICAL_GROUP]],
  CLINICAL_LEVELS,
  CLINICAL_GROUP
)


################################################################################
# 13. VALIDATE SAMPLE SCOPE
################################################################################

# Only NDMM and MGUS should receive composite assignments.
non_target <-
  !is_ndmm &
  !is_mgus

if (
  any(
    non_target &
    !is.na(
      phenotype[[GEP_LASSO_GROUP]]
    )
  )
) {
  stop(
    "Non-NDMM/non-MGUS samples received GEP LASSO composite assignments."
  )
}

if (
  any(
    non_target &
    !is.na(
      phenotype[[CLINICAL_GROUP]]
    )
  )
) {
  stop(
    "Non-NDMM/non-MGUS samples received clinical composite assignments."
  )
}


################################################################################
# 14. VALIDATE GEP70-HIGH OVERRIDE
################################################################################

gep70_high <-
  is_ndmm &
  phenotype$GEP70_binary == 1L

for (column in c(
  GEP_LASSO_GROUP,
  CLINICAL_GROUP
)) {
  
  if (
    any(
      gep70_high &
      as.character(
        phenotype[[column]]
      ) != "GEP70>0.66"
    )
  ) {
    stop(
      "GEP70-high override failed for ",
      column
    )
  }
}


################################################################################
# 15. VALIDATE SAMPLE NAMES
################################################################################

if (
  !identical(
    rownames(
      Biobase::pData(eset)
    ),
    Biobase::sampleNames(eset)
  )
) {
  stop(
    "pData row names do not match ExpressionSet sample names."
  )
}


################################################################################
# 16. QC SUMMARY
################################################################################

make_summary <- function(
    column,
    levels
) {
  
  tibble(
    Composite_Column = column,
    Group = levels,
    
    n = vapply(
      levels,
      function(x) {
        
        sum(
          as.character(
            phenotype[[column]]
          ) == x,
          na.rm = TRUE
        )
      },
      integer(1)
    )
  )
}

group_summary <- bind_rows(
  make_summary(
    GEP_LASSO_GROUP,
    GEP_LEVELS
  ),
  make_summary(
    CLINICAL_GROUP,
    CLINICAL_LEVELS
  )
)


################################################################################
# 17. QC VALUES
################################################################################

values_table <- tibble(
  PHENOTYPE_ROWNAME =
    phenotype_row_names,
  
  sample_group =
    phenotype$sample_group,
  
  GEP1q_binary_lasso =
    phenotype$GEP1q_binary_lasso,
  
  PolyIG_binary =
    phenotype$PolyIG_binary,
  
  GEP70_binary =
    phenotype$GEP70_binary,
  
  FISH1q_binary =
    phenotype$FISH1q_binary,
  
  IG_ARD_binary =
    phenotype$IG_ARD_binary,
  
  MGUS_Binary =
    phenotype$MGUS_Binary,
  
  GEP_RNAS_Groups_MGUS_Lasso =
    as.character(
      phenotype[[GEP_LASSO_GROUP]]
    ),
  
  Clinical_RNAS_Groups_MGUS =
    as.character(
      phenotype[[CLINICAL_GROUP]]
    )
)


################################################################################
# 18. SAVE QC TABLES
################################################################################

write_csv(
  group_summary,
  SUMMARY_FILE
)

write_csv(
  values_table,
  VALUES_FILE
)


################################################################################
# 19. SAVE NEW EXPRESSIONSET
################################################################################

saveRDS(
  eset,
  ESET_OUT_PATH
)

if (!file.exists(ESET_OUT_PATH)) {
  stop(
    "Output ExpressionSet was not created."
  )
}


################################################################################
# 20. SUMMARY
################################################################################

message(
  "\nStep 20B complete.",
  "\nInput = ", ESET_PATH,
  "\nOutput = ", ESET_OUT_PATH,
  "\nNDMM = ", sum(is_ndmm),
  "\nMGUS = ", sum(is_mgus),
  "\n\nGEP LASSO groups:"
)

print(
  table(
    phenotype[[GEP_LASSO_GROUP]],
    useNA = "ifany"
  )
)

message(
  "\nClinical groups:"
)

print(
  table(
    phenotype[[CLINICAL_GROUP]],
    useNA = "ifany"
  )
)