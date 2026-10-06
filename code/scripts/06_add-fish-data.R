################################################################################
# STEP 06 — ADD FISH 1q DATA
#
# Purpose:
#   1. Load the GEP1q-scored UAMS ExpressionSet.
#   2. Load expanded FISH 1q data.
#   3. Derive categorical and numeric FISH 1q copy-number fields.
#   4. Resolve duplicate FISH records.
#   5. Add FISH annotations to phenotype data.
#   6. Save the FISH-annotated ExpressionSet.
################################################################################

library(dplyr)
library(lubridate)


################################################################################
# 1. Files
################################################################################

ESET_FILE <- file.path(
  "data",
  "processed",
  "uams",
  "ESET_uams_all_gep1q_lasso.rds"
)

FISH_FILE <- file.path(
  "data",
  "raw",
  "uams",
  "fish-1q-1p-gep-bl-overlap.csv"
)

OUTPUT_FILE <- file.path(
  "data",
  "processed",
  "uams",
  "ESET_uams_all_gep1q_fish.rds"
)


################################################################################
# 2. Helper function for CHIPID
################################################################################

normalize_chipid <- function(x) {
  
  x <- toupper(
    trimws(
      as.character(x)
    )
  )
  
  x <- sub(
    "\\.0$",
    "",
    x
  )
  
  x <- sub(
    "^CHIPID",
    "",
    x
  )
  
  numeric_only <- grepl(
    "^[0-9]+$",
    x
  )
  
  x[numeric_only] <- paste0(
    "CHIPID",
    x[numeric_only]
  )
  
  x
}


################################################################################
# 3. Load ExpressionSet
################################################################################

eset <- readRDS(
  ESET_FILE
)

phenotype <- pData(
  eset
)


################################################################################
# 4. Load expanded FISH data
################################################################################

fish <- fread(
  FISH_FILE,
  data.table = FALSE
)

fish$chipid <- normalize_chipid(
  fish$CHIPID
)

fish$.row <- seq_len(
  nrow(fish)
)


################################################################################
# 5. Parse FISH date
################################################################################

fish$FISH_Date <- as.Date(
  parse_date_time(
    fish$FISH_Date,
    orders = c(
      "ymd",
      "mdy",
      "dmy"
    )
  )
)


################################################################################
# 6. Convert FISH fields to numeric
################################################################################

numeric_columns <- c(
  "FISH_1q_2copies",
  "FISH_1q_3copies",
  "FISH_1q_4plus",
  
  "1q21_Percent_Cells_Copy_0",
  "1q21_Percent_Cells_Copy_1",
  "1q21_Percent_Cells_Copy_2",
  "1q21_Percent_Cells_Copy_3",
  "1q21_Percent_Cells_Copy_4",
  "1q21_Percent_Cells_Copy_5",
  "1q21_Percent_Cells_Copy_6"
)

for (column in numeric_columns) {
  
  if (column %in% names(fish)) {
    
    fish[[column]] <- as.numeric(
      fish[[column]]
    )
  }
}


################################################################################
# 7. Derive FISH 1q copy-number class
#
# Worst-state precedence:
#   4plus > 3 > 2
################################################################################

fish$FISH_1q <- case_when(
  fish$FISH_1q_4plus == 1 ~ "4plus",
  fish$FISH_1q_3copies == 1 ~ "3",
  fish$FISH_1q_2copies == 1 ~ "2",
  TRUE ~ NA_character_
)

fish$FISH_1q_num <- case_when(
  fish$FISH_1q == "2" ~ 2,
  fish$FISH_1q == "3" ~ 3,
  fish$FISH_1q == "4plus" ~ 4,
  TRUE ~ NA_real_
)


################################################################################
# 8. Resolve duplicate FISH records
#
# Priority:
#   1. Baseline record
#   2. Earliest FISH date
#   3. Original row order
################################################################################

fish$Baseline <- as.numeric(
  fish$Baseline
)

fish <- fish |>
  arrange(
    desc(
      coalesce(
        Baseline,
        0
      )
    ),
    FISH_Date,
    .row
  ) |>
  group_by(
    chipid
  ) |>
  slice_head(
    n = 1
  ) |>
  ungroup()


################################################################################
# 9. Select FISH fields to add
################################################################################

fish_columns <- c(
  "chipid",
  "PATID",
  "Baseline",
  "TT",
  "TT_Simple",
  "Old_FISHCHIP",
  "FISH_Date",
  
  "FISH_1q_2copies",
  "FISH_1q_3copies",
  "FISH_1q_4plus",
  "FISH_1q",
  "FISH_1q_num",
  
  "1q21_Percent_Cells_Copy_0",
  "1q21_Percent_Cells_Copy_1",
  "1q21_Percent_Cells_Copy_2",
  "1q21_Percent_Cells_Copy_3",
  "1q21_Percent_Cells_Copy_4",
  "1q21_Percent_Cells_Copy_5",
  "1q21_Percent_Cells_Copy_6"
)

fish <- fish[
  ,
  intersect(
    fish_columns,
    names(fish)
  ),
  drop = FALSE
]


################################################################################
# 10. Join FISH data to phenotype data
################################################################################

phenotype$chipid <- normalize_chipid(
  phenotype$chipid
)

sample_order <- rownames(
  phenotype
)

fish$FISH_PATID <- fish$PATID
fish$PATID <- NULL

phenotype <- phenotype |>
  tibble::rownames_to_column(
    "sample_id"
  ) |>
  left_join(
    fish,
    by = "chipid"
  )

phenotype <- phenotype[
  match(
    sample_order,
    phenotype$sample_id
  ),
  ,
  drop = FALSE
]

rownames(
  phenotype
) <- phenotype$sample_id

phenotype$sample_id <- NULL


################################################################################
# 11. Update ExpressionSet
################################################################################

phenoData(
  eset
) <- AnnotatedDataFrame(
  phenotype
)


################################################################################
# 12. Save FISH-annotated ExpressionSet
################################################################################

saveRDS(
  eset,
  OUTPUT_FILE
)


################################################################################
# 13. Complete
################################################################################

cat(
  "ExpressionSet samples:",
  ncol(eset),
  "\n"
)

cat(
  "Samples with FISH 1q:",
  sum(
    !is.na(
      phenotype$FISH_1q
    )
  ),
  "\n"
)

cat(
  "FISH 1q2:",
  sum(
    phenotype$FISH_1q == "2",
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "FISH 1q3:",
  sum(
    phenotype$FISH_1q == "3",
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "FISH 1q4plus:",
  sum(
    phenotype$FISH_1q == "4plus",
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Saved:",
  OUTPUT_FILE,
  "\n"
)