################################################################################
# STEP 21B — MATCH RNAS LASSO GEP GROUPS TO RNBX BY PATID
#
# Purpose:
#   1. Match RNAS CD138 samples to RNBX whole-bone-marrow samples by PATID.
#   2. Match NDMM only to NDMM and MGUS only to MGUS.
#   3. Require exact one-to-one RNAS:RNBX matches.
#   4. Transfer GEP_RNAS_Groups_MGUS_Lasso to the matched RNBX row.
#   5. Save match QC tables and the updated ExpressionSet.
#
# Matching:
#   RNAS_CD138_NDMM -> RNBX_WB_NDMM
#   RNAS_CD138_MGUS -> RNBX_WB_MGUS
#
# The Step 20B ExpressionSet is not overwritten.
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
  "ESET_uams_all_gep1q_fish_polyig_pcme_binaries_lasso_groups.rds"
)

ESET_OUT_PATH <- file.path(
  "data",
  "processed",
  "uams",
  "ESET_uams_all_gep1q_fish_polyig_pcme_binaries_lasso_groups_rnbx-matched.rds"
)

QC_DIR <- file.path(
  "results",
  "qc"
)

PAIR_AUDIT_FILE <- file.path(
  QC_DIR,
  "uams_rnas_rnbx_pair_audit_step21b_lasso.csv"
)

MATCHED_PAIRS_FILE <- file.path(
  QC_DIR,
  "uams_rnas_rnbx_matched_pairs_step21b_lasso.csv"
)

MATCH_SUMMARY_FILE <- file.path(
  QC_DIR,
  "uams_rnas_rnbx_match_summary_step21b_lasso.csv"
)

OVERWRITE_EXISTING_MATCH_COLUMNS <- TRUE


################################################################################
# 2. DEFINITIONS
################################################################################

RNAS_GROUP <- "GEP_RNAS_Groups_MGUS_Lasso"

RNAS_NDMM <- "RNAS_CD138_NDMM"
RNBX_NDMM <- "RNBX_WB_NDMM"

RNAS_MGUS <- "RNAS_CD138_MGUS"
RNBX_MGUS <- "RNBX_WB_MGUS"

OUT_GROUP <- "GEP_RNBX_Groups_MGUS_Lasso"
OUT_STATUS <- "RNAS_to_RNBX_PATID_match_MGUS_Lasso"
OUT_TYPE <- "RNAS_to_RNBX_match_type_Lasso"

OUT_COLUMNS <- c(
  OUT_GROUP,
  OUT_STATUS,
  OUT_TYPE
)

GROUP_LEVELS <- c(
  "MGUS",
  "GEP70>0.66",
  "GEP1q-CN2/PolyIG-high",
  "GEP1q-CN2/PolyIG-low",
  "GEP1q-CN>=3/PolyIG-high",
  "GEP1q-CN>=3/PolyIG-low"
)

STATUS_LEVELS <- c(
  "matched_one_to_one",
  "no_matching_RNAS_PATID_and_type",
  "missing_PATID",
  "ambiguous_multiple_RNAS",
  "ambiguous_multiple_RNBX",
  "ambiguous_multiple_RNAS_and_RNBX"
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
    "Cannot find input ExpressionSet:\n",
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

row_names <- rownames(
  phenotype
)


################################################################################
# 5. REQUIRED COLUMNS
################################################################################

REQUIRED_COLUMNS <- c(
  "patid",
  "sample_group",
  "GEP_RNAS_Groups_MGUS_Lasso"
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
# 6. EXISTING OUTPUT COLUMNS
################################################################################

existing <- intersect(
  OUT_COLUMNS,
  colnames(phenotype)
)

if (
  length(existing) > 0 &&
  !OVERWRITE_EXISTING_MATCH_COLUMNS
) {
  stop(
    "Output columns already exist: ",
    paste(
      existing,
      collapse = ", "
    )
  )
}


################################################################################
# 7. DEFINE SOURCE VALUES
################################################################################

patid <- phenotype$patid

sample_group <- phenotype$sample_group

rnas_group <- as.character(
  phenotype$GEP_RNAS_Groups_MGUS_Lasso
)


################################################################################
# 8. DEFINE RNAS SOURCES AND RNBX TARGETS
################################################################################

is_rnas_ndmm <-
  sample_group == RNAS_NDMM &
  !is.na(rnas_group)

is_rnbx_ndmm <-
  sample_group == RNBX_NDMM

is_rnas_mgus <-
  sample_group == RNAS_MGUS &
  !is.na(rnas_group)

is_rnbx_mgus <-
  sample_group == RNBX_MGUS

is_rnas <-
  is_rnas_ndmm |
  is_rnas_mgus

is_rnbx <-
  is_rnbx_ndmm |
  is_rnbx_mgus


################################################################################
# 9. VALIDATE RNAS GROUPS
################################################################################

bad_groups <- setdiff(
  unique(
    rnas_group[
      is_rnas
    ]
  ),
  GROUP_LEVELS
)

if (length(bad_groups) > 0) {
  stop(
    "Unexpected RNAS group value(s): ",
    paste(
      bad_groups,
      collapse = ", "
    )
  )
}

if (
  any(
    is_rnas_mgus &
    rnas_group != "MGUS",
    na.rm = TRUE
  )
) {
  stop(
    "RNAS MGUS rows must have group MGUS."
  )
}

if (
  any(
    is_rnas_ndmm &
    rnas_group == "MGUS",
    na.rm = TRUE
  )
) {
  stop(
    "RNAS NDMM rows cannot have group MGUS."
  )
}


################################################################################
# 10. BUILD SOURCE TABLE
################################################################################

source_table <- tibble(
  SOURCE_ROW = which(
    is_rnas
  ),
  
  PATID = patid[
    is_rnas
  ],
  
  TYPE = ifelse(
    is_rnas_ndmm[
      is_rnas
    ],
    "NDMM",
    "MGUS"
  ),
  
  RNAS_GROUP = rnas_group[
    is_rnas
  ]
) |>
  mutate(
    KEY = ifelse(
      is.na(PATID),
      NA_character_,
      paste(
        TYPE,
        PATID,
        sep = "||"
      )
    )
  )


################################################################################
# 11. BUILD TARGET TABLE
################################################################################

target_table <- tibble(
  TARGET_ROW = which(
    is_rnbx
  ),
  
  PATID = patid[
    is_rnbx
  ],
  
  TYPE = ifelse(
    is_rnbx_ndmm[
      is_rnbx
    ],
    "NDMM",
    "MGUS"
  )
) |>
  mutate(
    KEY = ifelse(
      is.na(PATID),
      NA_character_,
      paste(
        TYPE,
        PATID,
        sep = "||"
      )
    )
  )


################################################################################
# 12. PAIR AUDIT
################################################################################

source_counts <- source_table |>
  filter(
    !is.na(KEY)
  ) |>
  group_by(
    KEY,
    PATID,
    TYPE
  ) |>
  summarise(
    N_RNAS = n(),
    .groups = "drop"
  )

target_counts <- target_table |>
  filter(
    !is.na(KEY)
  ) |>
  group_by(
    KEY,
    PATID,
    TYPE
  ) |>
  summarise(
    N_RNBX = n(),
    .groups = "drop"
  )

pair_audit <- full_join(
  source_counts,
  target_counts,
  by = c(
    "KEY",
    "PATID",
    "TYPE"
  )
) |>
  mutate(
    N_RNAS = coalesce(
      N_RNAS,
      0L
    ),
    
    N_RNBX = coalesce(
      N_RNBX,
      0L
    ),
    
    PAIR_STATUS = case_when(
      
      N_RNAS == 1L &
        N_RNBX == 1L ~
        "matched_one_to_one",
      
      N_RNAS == 0L &
        N_RNBX >= 1L ~
        "no_matching_RNAS_PATID_and_type",
      
      N_RNAS > 1L &
        N_RNBX == 1L ~
        "ambiguous_multiple_RNAS",
      
      N_RNAS == 1L &
        N_RNBX > 1L ~
        "ambiguous_multiple_RNBX",
      
      N_RNAS > 1L &
        N_RNBX > 1L ~
        "ambiguous_multiple_RNAS_and_RNBX",
      
      N_RNAS >= 1L &
        N_RNBX == 0L ~
        "RNAS_without_RNBX_target",
      
      TRUE ~
        "other"
    )
  )

write_csv(
  pair_audit,
  PAIR_AUDIT_FILE
)


################################################################################
# 13. CREATE EXACT ONE-TO-ONE MATCHES
################################################################################

matched_pairs <- source_table |>
  inner_join(
    pair_audit |>
      filter(
        PAIR_STATUS == "matched_one_to_one"
      ) |>
      select(
        KEY
      ),
    by = "KEY"
  ) |>
  inner_join(
    target_table,
    by = c(
      "KEY",
      "PATID",
      "TYPE"
    )
  ) |>
  select(
    TYPE,
    PATID,
    SOURCE_ROW,
    RNAS_GROUP,
    TARGET_ROW
  ) |>
  arrange(
    TYPE,
    PATID
  )

if (
  anyDuplicated(
    matched_pairs$SOURCE_ROW
  )
) {
  stop(
    "An RNAS row was assigned to more than one pair."
  )
}

if (
  anyDuplicated(
    matched_pairs$TARGET_ROW
  )
) {
  stop(
    "An RNBX row was assigned to more than one pair."
  )
}

write_csv(
  matched_pairs,
  MATCHED_PAIRS_FILE
)


################################################################################
# 14. CREATE RNBX OUTPUT VARIABLES
################################################################################

rnbx_group <- rep(
  NA_character_,
  nrow(phenotype)
)

match_status <- rep(
  NA_character_,
  nrow(phenotype)
)

match_type <- rep(
  NA_character_,
  nrow(phenotype)
)


################################################################################
# 15. ASSIGN MATCH TYPE
################################################################################

target_rows <- target_table$TARGET_ROW

match_type[
  target_rows
] <- target_table$TYPE


################################################################################
# 16. ASSIGN MATCH STATUS
################################################################################

missing_patid <- is.na(
  target_table$PATID
)

match_status[
  target_rows[
    missing_patid
  ]
] <- "missing_PATID"


target_status <- target_table |>
  left_join(
    pair_audit |>
      select(
        KEY,
        PAIR_STATUS
      ),
    by = "KEY"
  )

has_patid <- !is.na(
  target_status$PATID
)

match_status[
  target_status$TARGET_ROW[
    has_patid
  ]
] <- target_status$PAIR_STATUS[
  has_patid
]


################################################################################
# 17. TRANSFER RNAS GROUP TO MATCHED RNBX ROW
################################################################################

rnbx_group[
  matched_pairs$TARGET_ROW
] <- matched_pairs$RNAS_GROUP


################################################################################
# 18. ADD OUTPUT COLUMNS
################################################################################

phenotype[[OUT_GROUP]] <- factor(
  rnbx_group,
  levels = GROUP_LEVELS
)

phenotype[[OUT_STATUS]] <- factor(
  match_status,
  levels = STATUS_LEVELS
)

phenotype[[OUT_TYPE]] <- factor(
  match_type,
  levels = c(
    "NDMM",
    "MGUS"
  )
)

rownames(phenotype) <-
  row_names

Biobase::pData(eset) <-
  phenotype


################################################################################
# 19. VALIDATE RNBX SCOPE
################################################################################

if (
  any(
    !is_rnbx &
    !is.na(
      phenotype[[OUT_GROUP]]
    )
  )
) {
  stop(
    "A non-RNBX row received an RNBX group."
  )
}


################################################################################
# 20. VALIDATE UNMATCHED TARGETS
################################################################################

if (
  any(
    is_rnbx &
    as.character(
      phenotype[[OUT_STATUS]]
    ) != "matched_one_to_one" &
    !is.na(
      phenotype[[OUT_GROUP]]
    ),
    na.rm = TRUE
  )
) {
  stop(
    "An unmatched or ambiguous RNBX row received a group."
  )
}


################################################################################
# 21. VALIDATE MATCHED TARGETS
################################################################################

if (
  any(
    is_rnbx &
    as.character(
      phenotype[[OUT_STATUS]]
    ) == "matched_one_to_one" &
    is.na(
      phenotype[[OUT_GROUP]]
    ),
    na.rm = TRUE
  )
) {
  stop(
    "A matched RNBX row did not receive a group."
  )
}


################################################################################
# 22. VALIDATE MGUS MATCHES
################################################################################

if (
  any(
    is_rnbx_mgus &
    !is.na(
      phenotype[[OUT_GROUP]]
    ) &
    as.character(
      phenotype[[OUT_GROUP]]
    ) != "MGUS",
    na.rm = TRUE
  )
) {
  stop(
    "An RNBX MGUS row received a non-MGUS group."
  )
}


################################################################################
# 23. VALIDATE SAMPLE NAMES
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
# 24. MATCH SUMMARY
################################################################################

match_summary <- tibble(
  Step = "21b",
  
  Match_Type = c(
    "NDMM",
    "MGUS"
  ),
  
  Source_sample_group = c(
    RNAS_NDMM,
    RNAS_MGUS
  ),
  
  Target_sample_group = c(
    RNBX_NDMM,
    RNBX_MGUS
  ),
  
  Eligible_RNAS_sources = c(
    sum(
      is_rnas_ndmm
    ),
    sum(
      is_rnas_mgus
    )
  ),
  
  Eligible_RNBX_targets = c(
    sum(
      is_rnbx_ndmm
    ),
    sum(
      is_rnbx_mgus
    )
  ),
  
  Matched_one_to_one = c(
    sum(
      matched_pairs$TYPE == "NDMM"
    ),
    sum(
      matched_pairs$TYPE == "MGUS"
    )
  )
)

write_csv(
  match_summary,
  MATCH_SUMMARY_FILE
)


################################################################################
# 25. SAVE
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
# 26. COMPLETE
################################################################################

message(
  "\nStep 21B complete.",
  "\nMatched NDMM = ",
  sum(
    matched_pairs$TYPE == "NDMM"
  ),
  "\nMatched MGUS = ",
  sum(
    matched_pairs$TYPE == "MGUS"
  ),
  "\nTotal matched = ",
  nrow(
    matched_pairs
  ),
  "\nOutput = ",
  ESET_OUT_PATH
)