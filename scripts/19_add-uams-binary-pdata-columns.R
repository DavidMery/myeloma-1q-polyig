################################################################################
# STEP 19B — ADD UAMS BINARY pData COLUMNS
#
# Creates six phenotype columns:
#
# Expression-derived binaries:
#   GEP1q_binary_lasso : 0 = GEP1qcopy_lasso == "1q2"
#                        1 = GEP1qcopy_lasso == "1q3plus"
#
#   PolyIG_binary      : 0 = PolyIG_Score <= 11
#                        1 = PolyIG_Score > 11
#
#   GEP70_binary       : 0 = GEP70 <= 0.66
#                        1 = GEP70 > 0.66
#
# NDMM-only clinical binaries:
#   FISH1q_binary      : 0 = FISH_1q_num == 2
#                        1 = FISH_1q_num >= 3
#
#   IG_ARD_binary      : 0 = IG_ARD <= 0.29
#                        1 = IG_ARD > 0.29
#
# MGUS binary:
#   MGUS_Binary        : 1 = sample_group == "RNAS_CD138_MGUS"
#                        0 = all other nonmissing sample_group values
#                       NA = missing sample_group
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
  "ESET_uams_all_gep1q_fish_polyig_pcme.rds"
)

ESET_OUT_PATH <- file.path(
  "data",
  "processed",
  "uams",
  "ESET_uams_all_gep1q_fish_polyig_pcme_binaries_lasso.rds"
)

QC_DIR <- file.path(
  "results",
  "qc"
)

SUMMARY_FILE <- file.path(
  QC_DIR,
  "uams_binary_pdata_summary_step19b.csv"
)

VALUES_FILE <- file.path(
  QC_DIR,
  "uams_binary_pdata_values_step19b.csv"
)

EXPRESSION_BY_GROUP_FILE <- file.path(
  QC_DIR,
  "uams_expression_binaries_by_sample_group_step19b.csv"
)

OVERWRITE_EXISTING_BINARY_COLUMNS <- TRUE


################################################################################
# 2. DEFINITIONS
################################################################################

NDMM_GROUP <- "RNAS_CD138_NDMM"
MGUS_GROUP <- "RNAS_CD138_MGUS"

EXPRESSION_GROUPS <- c(
  "RNAS_CD138_MGUS",
  "RNAS_CD138_NDMM"
)

POLYIG_CUTOFF <- 11
GEP70_CUTOFF <- 0.66
IGARD_CUTOFF <- 0.29

OUT_GEP1Q_LASSO <- "GEP1q_binary_lasso"
OUT_POLYIG <- "PolyIG_binary"
OUT_GEP70 <- "GEP70_binary"
OUT_FISH1Q <- "FISH1q_binary"
OUT_IGARD <- "IG_ARD_binary"
OUT_MGUS <- "MGUS_Binary"

BINARY_COLUMNS <- c(
  OUT_GEP1Q_LASSO,
  OUT_POLYIG,
  OUT_GEP70,
  OUT_FISH1Q,
  OUT_IGARD,
  OUT_MGUS
)


################################################################################
# 3. HELPERS
################################################################################

make_dir <- function(path) {
  
  if (!dir.exists(path)) {
    dir.create(
      path,
      recursive = TRUE,
      showWarnings = FALSE
    )
  }
  
  if (!dir.exists(path)) {
    stop(
      "Could not create directory: ",
      path
    )
  }
}


threshold_binary <- function(x, cutoff) {
  
  ifelse(
    is.finite(x),
    as.integer(x > cutoff),
    NA_integer_
  )
}


check_binary <- function(x, name) {
  
  bad <- unique(
    x[
      !is.na(x) &
        !x %in% c(0L, 1L)
    ]
  )
  
  if (length(bad) > 0) {
    stop(
      name,
      " contains invalid values: ",
      paste(
        bad,
        collapse = ", "
      )
    )
  }
}


################################################################################
# 4. DIRECTORIES AND INPUT
################################################################################

if (!file.exists(ESET_PATH)) {
  stop(
    "Cannot find ExpressionSet:\n",
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

make_dir(
  dirname(ESET_OUT_PATH)
)

make_dir(
  QC_DIR
)

eset <- readRDS(
  ESET_PATH
)

if (!methods::is(eset, "ExpressionSet")) {
  stop(
    "Input file is not a Biobase ExpressionSet."
  )
}

phenotype <- as.data.frame(
  Biobase::pData(eset),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

original_rownames <- rownames(
  phenotype
)


################################################################################
# 5. REQUIRED SOURCE COLUMNS
################################################################################

REQUIRED_COLUMNS <- c(
  "sample_group",
  "GEP1qcopy_lasso",
  "PolyIG_Score",
  "GEP70",
  "FISH_1q_num",
  "IG_ARD"
)

missing_columns <- setdiff(
  REQUIRED_COLUMNS,
  colnames(phenotype)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required column(s): ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}


################################################################################
# 6. VALIDATE GEP1q COPY VALUES
################################################################################

invalid_gep1q_copy <- setdiff(
  unique(
    phenotype$GEP1qcopy_lasso[
      !is.na(
        phenotype$GEP1qcopy_lasso
      )
    ]
  ),
  c(
    "1q2",
    "1q3plus"
  )
)

if (length(invalid_gep1q_copy) > 0) {
  stop(
    "GEP1qcopy_lasso contains unexpected value(s): ",
    paste(
      invalid_gep1q_copy,
      collapse = ", "
    )
  )
}


################################################################################
# 7. EXISTING OUTPUT COLUMNS
################################################################################

existing <- intersect(
  BINARY_COLUMNS,
  colnames(phenotype)
)

if (
  length(existing) > 0 &&
  !OVERWRITE_EXISTING_BINARY_COLUMNS
) {
  stop(
    "These output columns already exist: ",
    paste(
      existing,
      collapse = ", "
    )
  )
}


################################################################################
# 8. DEFINE SAMPLE SCOPES
################################################################################

is_expression_scope <-
  !is.na(phenotype$sample_group) &
  phenotype$sample_group %in% EXPRESSION_GROUPS

is_ndmm <-
  !is.na(phenotype$sample_group) &
  phenotype$sample_group == NDMM_GROUP


################################################################################
# 9. CREATE EMPTY BINARY VARIABLES
################################################################################

gep1q_binary_lasso <- rep(
  NA_integer_,
  nrow(phenotype)
)

polyig_binary <- rep(
  NA_integer_,
  nrow(phenotype)
)

gep70_binary <- rep(
  NA_integer_,
  nrow(phenotype)
)

fish1q_binary <- rep(
  NA_integer_,
  nrow(phenotype)
)

igard_binary <- rep(
  NA_integer_,
  nrow(phenotype)
)

mgus_binary <- rep(
  NA_integer_,
  nrow(phenotype)
)


################################################################################
# 10. EXPRESSION-DERIVED BINARIES
################################################################################

# GEP1q LASSO copy class
gep1q_binary_lasso[
  is_expression_scope &
    phenotype$GEP1qcopy_lasso == "1q2"
] <- 0L

gep1q_binary_lasso[
  is_expression_scope &
    phenotype$GEP1qcopy_lasso == "1q3plus"
] <- 1L


# PolyIG
polyig_binary[is_expression_scope] <-
  threshold_binary(
    phenotype$PolyIG_Score,
    POLYIG_CUTOFF
  )[is_expression_scope]


# GEP70
gep70_binary[is_expression_scope] <-
  threshold_binary(
    phenotype$GEP70,
    GEP70_CUTOFF
  )[is_expression_scope]


################################################################################
# 11. NDMM-ONLY CLINICAL BINARIES
################################################################################

# FISH 1q
fish1q_binary[
  is_ndmm &
    is.finite(phenotype$FISH_1q_num) &
    phenotype$FISH_1q_num == 2
] <- 0L

fish1q_binary[
  is_ndmm &
    is.finite(phenotype$FISH_1q_num) &
    phenotype$FISH_1q_num >= 3
] <- 1L


# IG_ARD
igard_binary[is_ndmm] <-
  threshold_binary(
    phenotype$IG_ARD,
    IGARD_CUTOFF
  )[is_ndmm]


################################################################################
# 12. MGUS BINARY
################################################################################

mgus_binary[
  !is.na(
    phenotype$sample_group
  )
] <- as.integer(
  phenotype$sample_group[
    !is.na(
      phenotype$sample_group
    )
  ] == MGUS_GROUP
)


################################################################################
# 13. ADD BINARY COLUMNS TO pData
################################################################################

phenotype[[OUT_GEP1Q_LASSO]] <-
  gep1q_binary_lasso

phenotype[[OUT_POLYIG]] <-
  polyig_binary

phenotype[[OUT_GEP70]] <-
  gep70_binary

phenotype[[OUT_FISH1Q]] <-
  fish1q_binary

phenotype[[OUT_IGARD]] <-
  igard_binary

phenotype[[OUT_MGUS]] <-
  mgus_binary

rownames(phenotype) <-
  original_rownames

Biobase::pData(eset) <-
  phenotype


################################################################################
# 14. VALIDATE BINARY COLUMNS
################################################################################

for (x in BINARY_COLUMNS) {
  
  check_binary(
    phenotype[[x]],
    x
  )
}


################################################################################
# 15. VALIDATE EXPRESSION-DERIVED SAMPLE SCOPE
################################################################################

expression_columns <- c(
  OUT_GEP1Q_LASSO,
  OUT_POLYIG,
  OUT_GEP70
)

for (x in expression_columns) {
  
  if (
    any(
      !is_expression_scope &
      !is.na(
        phenotype[[x]]
      )
    )
  ) {
    stop(
      x,
      " contains values outside the specified RNAS CD138 groups."
    )
  }
}


################################################################################
# 16. VALIDATE NDMM-ONLY SAMPLE SCOPE
################################################################################

for (x in c(
  OUT_FISH1Q,
  OUT_IGARD
)) {
  
  if (
    any(
      !is_ndmm &
      !is.na(
        phenotype[[x]]
      )
    )
  ) {
    stop(
      x,
      " contains values outside NDMM."
    )
  }
}


################################################################################
# 17. VALIDATE SAMPLE NAMES
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
    "pData row names no longer match ExpressionSet sample names."
  )
}


################################################################################
# 18. QC SUMMARY
################################################################################

summary_table <- tibble(
  Column = BINARY_COLUMNS,
  
  n_zero = vapply(
    BINARY_COLUMNS,
    function(x) {
      sum(
        phenotype[[x]] == 0L,
        na.rm = TRUE
      )
    },
    integer(1)
  ),
  
  n_one = vapply(
    BINARY_COLUMNS,
    function(x) {
      sum(
        phenotype[[x]] == 1L,
        na.rm = TRUE
      )
    },
    integer(1)
  ),
  
  n_missing = vapply(
    BINARY_COLUMNS,
    function(x) {
      sum(
        is.na(
          phenotype[[x]]
        )
      )
    },
    integer(1)
  )
)


################################################################################
# 19. QC VALUES TABLE
################################################################################

values_table <- tibble(
  PHENOTYPE_ROWNAME =
    original_rownames,
  
  sample_group =
    phenotype$sample_group,
  
  GEP1qcopy_lasso =
    phenotype$GEP1qcopy_lasso,
  
  GEP1q_binary_lasso =
    phenotype[[OUT_GEP1Q_LASSO]],
  
  PolyIG_Score =
    phenotype$PolyIG_Score,
  
  PolyIG_binary =
    phenotype[[OUT_POLYIG]],
  
  GEP70 =
    phenotype$GEP70,
  
  GEP70_binary =
    phenotype[[OUT_GEP70]],
  
  FISH_1q_num =
    phenotype$FISH_1q_num,
  
  FISH1q_binary =
    phenotype[[OUT_FISH1Q]],
  
  IG_ARD =
    phenotype$IG_ARD,
  
  IG_ARD_binary =
    phenotype[[OUT_IGARD]],
  
  MGUS_Binary =
    phenotype[[OUT_MGUS]]
)


################################################################################
# 20. EXPRESSION BINARIES BY SAMPLE GROUP
################################################################################

expression_by_group <- values_table |>
  filter(
    sample_group %in% EXPRESSION_GROUPS
  ) |>
  group_by(
    sample_group
  ) |>
  summarise(
    n = n(),
    
    GEP1q_lasso_0 = sum(
      GEP1q_binary_lasso == 0L,
      na.rm = TRUE
    ),
    
    GEP1q_lasso_1 = sum(
      GEP1q_binary_lasso == 1L,
      na.rm = TRUE
    ),
    
    PolyIG_0 = sum(
      PolyIG_binary == 0L,
      na.rm = TRUE
    ),
    
    PolyIG_1 = sum(
      PolyIG_binary == 1L,
      na.rm = TRUE
    ),
    
    GEP70_0 = sum(
      GEP70_binary == 0L,
      na.rm = TRUE
    ),
    
    GEP70_1 = sum(
      GEP70_binary == 1L,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )


################################################################################
# 21. SAVE QC TABLES
################################################################################

write_csv(
  summary_table,
  SUMMARY_FILE
)

write_csv(
  values_table,
  VALUES_FILE
)

write_csv(
  expression_by_group,
  EXPRESSION_BY_GROUP_FILE
)


################################################################################
# 22. SAVE EXPRESSIONSET
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
# 23. COMPLETE
################################################################################

message(
  "\nStep 19B complete.",
  "\nOutput = ", ESET_OUT_PATH,
  "\nQC summary = ", SUMMARY_FILE,
  "\nQC values = ", VALUES_FILE,
  "\nExpression-group QC = ", EXPRESSION_BY_GROUP_FILE,
  "\n\nBinary counts:"
)

print(
  summary_table
)