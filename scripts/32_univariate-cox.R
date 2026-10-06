################################################################################
# STEP 32 — UNIVARIATE COX REGRESSION
#
# Endpoint:
#   TTP from first transplant
#
# Cohort:
#   - baseline RNAS_CD138_NDMM
#   - TT3a, TT3b, TT4_S-TT3
#   - GEP70 <= 0.66
#
# Clinical predictors:
#   - names from Clinical_variables.csv
#   - values from large UAMS clinical table
#
# Molecular predictors:
#   - PolyIG_Score
#   - IG_ARD
#   - FISH_1q_num
#   - GEP1qcopy_lasso
#   - dated FISH binaries
#
# Output:
#   HR, 95% CI, p value, BH FDR, N, events, type, coding
################################################################################

suppressPackageStartupMessages({
  library(Biobase)
  library(dplyr)
  library(readr)
  library(survival)
})

################################################################################
# 1. FILES / SETTINGS
################################################################################

ESET_PATH <- file.path(
  "data", "processed", "uams",
  "ESET_uams_all_gep1q_fish_polyig_pcme_binaries_lasso_groups_rnbx-matched.rds"
)

FACTORS_PATH <- file.path(
  "data", "raw", "uams", "clinical",
  "Clinical_variables.csv"
)

CLINICAL_PATH <- file.path(
  "data", "raw", "uams", "clinical",
  paste0(
    "xlsjs_all_090624_n154395_v68_2026-10-05.csv"
  )
)

OUT_DIR <- file.path(
  "results", "table1", "Univariate"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# Exact pData columns
PATID_COL <- "patid"
SAMPLE_GROUP_COL <- "sample_group"
TT_COL <- "TT_On_for_Tx1"
GEP70_COL <- "GEP70"

TIME_COL <- "YearsTTP"
EVENT_COL <- "CensTTP"

POLYIG_COL <- "PolyIG_Score"
UNIGH_COL <- "IG_ARD"
FISH1Q_COL <- "FISH_1q_num"
GEP1Q_COL <- "GEP1qcopy_lasso"

TT_KEEP <- c(
  "TT3a",
  "TT3b",
  "TT4_S-TT3"
)

GEP70_CUTOFF <- 0.66
POLYIG_CUTOFF <- 11
UNIGH_CUTOFF <- 0.29

MIN_N <- 20
MIN_EVENTS <- 5

FISH_BINARY_COLS <- c(
  "FISH_1q21_2026-07-28",
  "FISH_1pdel_2026-07-29",
  "FISH_17pdel_2026-07-28",
  "FISH_13qdel_2026-07-28"
)

################################################################################
# 2. HELPERS
################################################################################

num <- function(x) {
  suppressWarnings(
    as.numeric(
      as.character(x)
    )
  )
}

clean_id <- function(x) {
  x <- trimws(
    as.character(x)
  )
  
  x[
    is.na(x) |
      x == ""
  ] <- NA_character_
  
  x
}

################################################################################
# 3. LOAD CLINICAL VARIABLE LIST
################################################################################

factors_raw <- read.csv(
  FACTORS_PATH,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

CLINICAL_VARS <- as.vector(
  t(
    factors_raw
  )
)

CLINICAL_VARS <- trimws(
  as.character(
    CLINICAL_VARS
  )
)

CLINICAL_VARS <- unique(
  CLINICAL_VARS[
    !is.na(CLINICAL_VARS) &
      CLINICAL_VARS != ""
  ]
)

cat(
  "Clinical predictors requested:",
  length(CLINICAL_VARS),
  "\n"
)

################################################################################
# 4. LOAD CLINICAL TABLE
################################################################################

clinical <- read.csv(
  CLINICAL_PATH,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

bad_names <- which(
  is.na(
    names(clinical)
  ) |
    trimws(
      names(clinical)
    ) == ""
)

if (length(bad_names)) {
  names(clinical)[bad_names] <- paste0(
    "Unnamed_",
    bad_names
  )
}

names(clinical) <- make.unique(
  names(clinical)
)

if (!"PATID" %in% names(clinical)) {
  stop(
    "Required clinical-table column 'PATID' was not found."
  )
}

clinical <- clinical |>
  mutate(
    PATID = clean_id(
      PATID
    )
  )

missing_clinical_vars <- setdiff(
  CLINICAL_VARS,
  names(clinical)
)

if (length(missing_clinical_vars)) {
  cat(
    "\nClinical variables absent from clinical table:\n"
  )
  print(
    missing_clinical_vars
  )
}

keep_clinical <- intersect(
  c(
    "PATID",
    CLINICAL_VARS
  ),
  names(clinical)
)

clinical <- clinical |>
  select(
    all_of(
      keep_clinical
    )
  ) |>
  filter(
    !is.na(PATID)
  ) |>
  distinct(
    PATID,
    .keep_all = TRUE
  )

################################################################################
# 5. LOAD EXPRESSIONSET / pData
################################################################################

if (!file.exists(ESET_PATH)) {
  stop(
    "ESET not found:\n",
    ESET_PATH
  )
}

eset <- readRDS(
  ESET_PATH
)

if (!methods::is(eset, "ExpressionSet")) {
  stop(
    "ESET_PATH is not an ExpressionSet."
  )
}

ph <- as.data.frame(
  Biobase::pData(eset),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

################################################################################
# 6. CHECK EXACT REQUIRED pData COLUMNS
################################################################################

required_pdata_columns <- c(
  PATID_COL,
  SAMPLE_GROUP_COL,
  TT_COL,
  GEP70_COL,
  TIME_COL,
  EVENT_COL,
  POLYIG_COL,
  UNIGH_COL,
  FISH1Q_COL,
  GEP1Q_COL
)

missing_pdata_columns <- setdiff(
  required_pdata_columns,
  names(ph)
)

if (length(missing_pdata_columns)) {
  stop(
    "Missing required pData columns: ",
    paste(
      missing_pdata_columns,
      collapse = ", "
    )
  )
}

################################################################################
# 7. BUILD ANALYSIS COHORT
################################################################################

df <- ph |>
  mutate(
    PATID = clean_id(
      .data[[PATID_COL]]
    ),
    
    sample_group_analysis =
      .data[[SAMPLE_GROUP_COL]],
    
    TT_analysis =
      .data[[TT_COL]],
    
    GEP70_analysis =
      num(
        .data[[GEP70_COL]]
      ),
    
    time =
      num(
        .data[[TIME_COL]]
      ),
    
    event =
      num(
        .data[[EVENT_COL]]
      )
  ) |>
  filter(
    sample_group_analysis ==
      "RNAS_CD138_NDMM",
    
    TT_analysis %in%
      TT_KEEP,
    
    is.finite(
      GEP70_analysis
    ),
    
    GEP70_analysis <=
      GEP70_CUTOFF,
    
    is.finite(
      time
    ),
    
    time >= 0,
    
    event %in%
      c(0, 1),
    
    !is.na(
      PATID
    )
  ) |>
  distinct(
    PATID,
    .keep_all = TRUE
  )

################################################################################
# 8. LEFT-JOIN CLINICAL VARIABLES
#
# pData defines the analysis cohort.
# Missing clinical-table matches remain NA and are handled per Cox model.
################################################################################

clinical_collisions <- intersect(
  setdiff(
    names(df),
    "PATID"
  ),
  setdiff(
    names(clinical),
    "PATID"
  )
)

if (length(clinical_collisions)) {
  clinical <- clinical |>
    select(
      -all_of(
        clinical_collisions
      )
    )
}

df <- left_join(
  df,
  clinical,
  by = "PATID"
)

if (
  nrow(df) !=
  n_distinct(
    df$PATID
  )
) {
  stop(
    "Duplicated PATIDs after clinical-table join."
  )
}

################################################################################
# 9. DERIVE MOLECULAR PREDICTORS
################################################################################

df <- df |>
  mutate(
    PolyIG_continuous =
      num(
        .data[[POLYIG_COL]]
      ),
    
    PolyIG_high_bin =
      if_else(
        is.na(
          PolyIG_continuous
        ),
        NA_integer_,
        as.integer(
          PolyIG_continuous >
            POLYIG_CUTOFF
        )
      ),
    
    IG_ARD_continuous =
      num(
        .data[[UNIGH_COL]]
      ),
    
    unIGH_high_bin =
      if_else(
        is.na(
          IG_ARD_continuous
        ),
        NA_integer_,
        as.integer(
          IG_ARD_continuous >
            UNIGH_CUTOFF
        )
      ),
    
    FISH1q_3plus_bin =
      case_when(
        num(
          .data[[FISH1Q_COL]]
        ) == 2 ~ 0,
        
        num(
          .data[[FISH1Q_COL]]
        ) >= 3 ~ 1,
        
        TRUE ~ NA_real_
      ),
    
    GEP1q_3plus_bin =
      case_when(
        .data[[GEP1Q_COL]] ==
          "1q2" ~ 0,
        
        .data[[GEP1Q_COL]] ==
          "1q3plus" ~ 1,
        
        TRUE ~ NA_real_
      )
  )

################################################################################
# 10. DATED FISH BINARIES
################################################################################

available_fish <- intersect(
  FISH_BINARY_COLS,
  names(df)
)

for (v in available_fish) {
  df[[v]] <- num(
    df[[v]]
  )
}

################################################################################
# 11. VARIABLES TO TEST
################################################################################

MOLECULAR_VARS <- c(
  "PolyIG_continuous",
  "PolyIG_high_bin",
  "IG_ARD_continuous",
  "unIGH_high_bin",
  "FISH1q_3plus_bin",
  "GEP1q_3plus_bin",
  available_fish
)

VARS_TO_TEST <- unique(
  c(
    intersect(
      CLINICAL_VARS,
      names(df)
    ),
    intersect(
      MOLECULAR_VARS,
      names(df)
    )
  )
)

cat(
  "\nTotal predictors to test:",
  length(VARS_TO_TEST),
  "\n"
)

################################################################################
# 12. VARIABLE CODING
################################################################################

CODING_MAP <- c(
  PolyIG_continuous =
    "Continuous PolyIG_Score; HR per 1-unit increase",
  
  PolyIG_high_bin =
    paste0(
      "1 = PolyIG_Score > ",
      POLYIG_CUTOFF,
      "; 0 = PolyIG_Score <= ",
      POLYIG_CUTOFF
    ),
  
  IG_ARD_continuous =
    "Continuous IG_ARD; HR per 1-unit increase",
  
  unIGH_high_bin =
    paste0(
      "1 = IG_ARD > ",
      UNIGH_CUTOFF,
      "; 0 = IG_ARD <= ",
      UNIGH_CUTOFF
    ),
  
  FISH1q_3plus_bin =
    "1 = FISH_1q_num >=3 copies; 0 = FISH_1q_num = 2 copies",
  
  GEP1q_3plus_bin =
    "1 = GEP1qcopy_lasso = 1q3plus; 0 = GEP1qcopy_lasso = 1q2",
  
  "FISH_1q21_2026-07-28" =
    "1 = 1q21 >=3 copies; 0 = 2 copies",
  
  "FISH_1pdel_2026-07-29" =
    "1 = 1p deletion; 0 = no 1p deletion",
  
  "FISH_17pdel_2026-07-28" =
    "1 = 17p deletion; 0 = no 17p deletion",
  
  "FISH_13qdel_2026-07-28" =
    "1 = 13q deletion; 0 = no 13q deletion"
)

################################################################################
# 13. UNIVARIATE COX FUNCTION
################################################################################

fit_univariate_cox <- function(var) {
  
  d <- df |>
    select(
      time,
      event,
      x = all_of(var)
    ) |>
    filter(
      !is.na(time),
      !is.na(event),
      !is.na(x)
    )
  
  n <- nrow(d)
  events <- sum(
    d$event == 1
  )
  
  empty <- function(
    type = NA_character_,
    note = ""
  ) {
    tibble(
      Variable = var,
      Type = type,
      N = n,
      Events = events,
      Term = NA_character_,
      HR = NA_real_,
      CI_low = NA_real_,
      CI_high = NA_real_,
      p = NA_real_,
      Note = note
    )
  }
  
  if (n < MIN_N) {
    return(
      empty(
        note = paste0(
          "N < ",
          MIN_N
        )
      )
    )
  }
  
  if (events < MIN_EVENTS) {
    return(
      empty(
        note = paste0(
          "Events < ",
          MIN_EVENTS
        )
      )
    )
  }
  
  x_num <- num(
    d$x
  )
  
  nonmissing_numeric <- sum(
    !is.na(
      x_num
    )
  )
  
  unique_numeric <- sort(
    unique(
      x_num[
        !is.na(
          x_num
        )
      ]
    )
  )
  
  ##########################################################################
  # Binary 0/1
  ##########################################################################
  
  if (
    nonmissing_numeric == n &&
    all(
      unique_numeric %in%
      c(0, 1)
    )
  ) {
    
    if (
      length(
        unique_numeric
      ) < 2
    ) {
      return(
        empty(
          "binary",
          "Only one level"
        )
      )
    }
    
    d$x <- factor(
      x_num,
      levels = c(
        0,
        1
      )
    )
    
    fit <- coxph(
      Surv(
        time,
        event
      ) ~ x,
      data = d
    )
    
    s <- summary(
      fit
    )
    
    return(
      tibble(
        Variable = var,
        Type = "binary",
        N = n,
        Events = events,
        Term = "1 vs 0",
        HR = s$coef[
          1,
          "exp(coef)"
        ],
        CI_low = s$conf.int[
          1,
          "lower .95"
        ],
        CI_high = s$conf.int[
          1,
          "upper .95"
        ],
        p = s$coef[
          1,
          "Pr(>|z|)"
        ],
        Note = ""
      )
    )
  }
  
  ##########################################################################
  # Continuous numeric
  ##########################################################################
  
  if (
    nonmissing_numeric /
    n >= 0.90
  ) {
    
    d$x <- x_num
    
    d <- d |>
      filter(
        !is.na(x)
      )
    
    if (
      n_distinct(
        d$x
      ) < 2
    ) {
      return(
        empty(
          "continuous",
          "No variation"
        )
      )
    }
    
    fit <- coxph(
      Surv(
        time,
        event
      ) ~ x,
      data = d
    )
    
    s <- summary(
      fit
    )
    
    return(
      tibble(
        Variable = var,
        Type = "continuous",
        N = nrow(d),
        Events = sum(
          d$event == 1
        ),
        Term =
          "Per 1-unit increase",
        HR = s$coef[
          1,
          "exp(coef)"
        ],
        CI_low = s$conf.int[
          1,
          "lower .95"
        ],
        CI_high = s$conf.int[
          1,
          "upper .95"
        ],
        p = s$coef[
          1,
          "Pr(>|z|)"
        ],
        Note = ""
      )
    )
  }
  
  ##########################################################################
  # Categorical
  ##########################################################################
  
  d$x <- droplevels(
    factor(
      as.character(
        d$x
      )
    )
  )
  
  if (
    nlevels(
      d$x
    ) < 2
  ) {
    return(
      empty(
        "categorical",
        "Only one level"
      )
    )
  }
  
  fit <- coxph(
    Surv(
      time,
      event
    ) ~ x,
    data = d
  )
  
  s <- summary(
    fit
  )
  
  tibble(
    Variable = var,
    Type = paste0(
      "categorical_",
      nlevels(
        d$x
      ),
      "_level"
    ),
    N = nrow(d),
    Events = sum(
      d$event == 1
    ),
    Term =
      "Global Wald test",
    HR = NA_real_,
    CI_low = NA_real_,
    CI_high = NA_real_,
    p = unname(
      s$waldtest[
        "pvalue"
      ]
    ),
    Note =
      "Global p value; individual level HRs not shown"
  )
}

################################################################################
# 14. RUN ALL UNIVARIATE COX MODELS
################################################################################

results <- bind_rows(
  lapply(
    VARS_TO_TEST,
    fit_univariate_cox
  )
) |>
  mutate(
    FDR_BH =
      p.adjust(
        p,
        method = "BH"
      ),
    
    Coding =
      unname(
        CODING_MAP[
          Variable
        ]
      ),
    
    Coding =
      if_else(
        is.na(Coding),
        "See source clinical variable",
        Coding
      ),
    
    HR_95CI =
      if_else(
        is.finite(HR),
        sprintf(
          "%.2f (%.2f–%.2f)",
          HR,
          CI_low,
          CI_high
        ),
        NA_character_
      )
  ) |>
  select(
    Variable,
    Coding,
    Type,
    N,
    Events,
    Term,
    HR,
    CI_low,
    CI_high,
    HR_95CI,
    p,
    FDR_BH,
    Note
  )

################################################################################
# 15. SAVE RESULTS
################################################################################

RESULT_FILE <- file.path(
  OUT_DIR,
  "Step32_UnivariateCox_TT3a_TT3b_TT4S_GEP70standard_TTP.csv"
)

readr::write_excel_csv(
  results,
  RESULT_FILE
)

################################################################################
# 16. SAVE ANALYSIS COHORT
################################################################################

COHORT_FILE <- file.path(
  OUT_DIR,
  "Step32_UnivariateCox_analysis-cohort_TTP.csv"
)

readr::write_csv(
  df,
  COHORT_FILE
)

################################################################################
# 17. FINAL SUMMARY
################################################################################

cat(
  "\n============================================================\n",
  "STEP 32 — UNIVARIATE COX COMPLETE\n",
  "============================================================\n",
  "Endpoint: TTP\n",
  "Analysis cohort N: ", nrow(df), "\n",
  "TTP events: ", sum(df$event == 1, na.rm = TRUE), "\n",
  "Predictors tested: ", length(VARS_TO_TEST), "\n",
  "Results: ", RESULT_FILE, "\n",
  "Analysis cohort: ", COHORT_FILE, "\n",
  sep = ""
)

cat(
  "\nTop results by p value:\n"
)

print(
  results |>
    filter(
      is.finite(p)
    ) |>
    arrange(
      p
    ) |>
    select(
      Variable,
      Type,
      N,
      Events,
      HR_95CI,
      p,
      FDR_BH
    ) |>
    head(
      20
    )
)