################################################################################
# STEP 25B — COMPARE TT2 THALIDOMIDE VS NO THALIDOMIDE
# 1q2 / PolyIG-High patients only
#
# Purpose:
#   Compare the favorable GEP-derived group between:
#     1. TT2_NoThal
#     2. TT2_Thal
#
# Cohort:
#   - Baseline CD138-selected NDMM
#   - GEP70-standard-risk: GEP70 <= 0.66
#   - Frozen penalized GEP1q class: GEP1qcopy_lasso == "1q2"
#   - PolyIG High: PolyIG > 11
#   - TT2_NoThal or TT2_Thal
#
# Analysis:
#   - Kaplan–Meier cumulative relapse/progression probability: 1 - S(t)
#   - Cox proportional-hazards model
#       Reference: TT2_NoThal
#       Reported HR: TT2_Thal vs TT2_NoThal
#   - Two-sided Cox Wald p-value
#   - Two-sided log-rank p-value
#
# Primary outputs:
#   results/figures/Fig3/Fig3A_TT2_Thal_vs_NoThal_1q2_PolyIGHigh/
#
# Suggested script filename:
#   code/scripts/uams/25b_compare_TT_Thal_NoThal.R
################################################################################

################################################################################
# 1. REQUIRED PACKAGES
################################################################################

REQUIRED_PACKAGES <- c(
  "Biobase",
  "survival",
  "survminer",
  "dplyr",
  "tibble",
  "readr",
  "ggplot2",
  "ggpubr"
)

missing_packages <- REQUIRED_PACKAGES[
  !vapply(
    REQUIRED_PACKAGES,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    "The following required package(s) are not installed: ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them in code/scripts/uams/00_prepare-environment.R before running this script."
  )
}

suppressPackageStartupMessages({
  library(Biobase)
  library(survival)
  library(survminer)
  library(dplyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(ggpubr)
})

################################################################################
# 2. CONFIGURATION
################################################################################

set.seed(123)

################################################################################
# 2A. VERIFY PIPELINE ROOTS
################################################################################

if (!exists("DATA_ROOT")) {
  stop(
    "DATA_ROOT is not defined. ",
    "Run this script through code/run_pipeline.R."
  )
}

if (!exists("RESULTS_ROOT")) {
  stop(
    "RESULTS_ROOT is not defined. ",
    "Run this script through code/run_pipeline.R."
  )
}

if (!dir.exists(DATA_ROOT)) {
  stop(
    "DATA_ROOT does not exist: ",
    DATA_ROOT
  )
}

dir.create(
  RESULTS_ROOT,
  recursive = TRUE,
  showWarnings = FALSE
)


STEP_ID <- "25b"
FIGURE_ID <- "Fig3A_TT2_Thal_vs_NoThal_1q2_PolyIGHigh"

ESET_PATH <- file.path(
  DATA_ROOT,
  "processed",
  "ESET_uams_all_gep1q_fish_polyig_pcme.rds"
)

FIGURE_DIR <- file.path(
  RESULTS_ROOT,
  "figures",
  "Fig3",
  FIGURE_ID
)

SOURCE_DATA_DIR <- file.path(
  FIGURE_DIR,
  "source-data"
)

ENDPOINT <- "TTP"

# Frozen GEP-derived classification settings.
GEP1Q_COLUMN_CANDIDATES <- c(
  "GEP1qcopy_lasso"
)

POLYIG_COLUMN_CANDIDATES <- c(
  "DEM_PolyIG",
  "PolyIG_Score_Log2",
  "PolyIG_score",
  "PolyIG"
)

GEP70_COLUMN_CANDIDATES <- c(
  "GEP70",
  "GEP70new"
)

TT_PROTOCOL_COLUMN_CANDIDATES <- c(
  "TT_On_for_Tx1",
  "TT_Guido"
)

POLYIG_CUT <- 11
GEP70_THRESHOLD <- 0.66

TT_REFERENCE <- "TT2_NoThal"
TT_COMPARATOR <- "TT2_Thal"
TT_KEEP <- c(
  TT_REFERENCE,
  TT_COMPARATOR
)

# TRUE means the selected endpoint field uses:
#   1 = relapse/progression/death event
#   0 = censored
EVENT_VALUE_ONE_IS_EVENT <- TRUE

# Plot settings.
BREAK_TIME_BY_YEARS <- 2
PLOT_WIDTH_IN <- 9
PLOT_HEIGHT_IN <- 8
PNG_DPI <- 300

SHOW_CONFIDENCE_INTERVAL <- FALSE
SHOW_CENSOR_MARKS <- TRUE

# Both curves remain blue, but different shades and line types make them readable.
PROTOCOL_COLORS <- c(
  "TT2_NoThal" = "steelblue2",
  "TT2_Thal" = "steelblue4"
)

PROTOCOL_LINE_TYPES <- c(
  "TT2_NoThal" = "solid",
  "TT2_Thal" = "dashed"
)

################################################################################
# 3. HELPER FUNCTIONS
################################################################################

create_dir_or_stop <- function(path) {
  if (!dir.exists(path)) {
    created <- dir.create(
      path,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    if (!created && !dir.exists(path)) {
      stop(
        "Could not create output directory:\n",
        normalizePath(
          path,
          winslash = "/",
          mustWork = FALSE
        )
      )
    }
  }
  
  invisible(path)
}

clean_text <- function(x) {
  output <- trimws(
    as.character(x)
  )
  
  output[
    is.na(output) |
      output %in% c(
        "",
        "NA",
        "N/A",
        "<NA>",
        "NaN"
      )
  ] <- NA_character_
  
  output
}

to_numeric_safe <- function(x) {
  suppressWarnings(
    as.numeric(
      as.character(x)
    )
  )
}

find_first_column <- function(
    data,
    candidates,
    label,
    required = TRUE
) {
  available <- colnames(data)
  
  hit <- match(
    tolower(candidates),
    tolower(available)
  )
  
  hit <- hit[!is.na(hit)]
  
  if (length(hit) == 0) {
    if (!required) {
      return(NULL)
    }
    
    stop(
      "Could not find a ",
      label,
      " column.\nCandidates: ",
      paste(candidates, collapse = ", "),
      "\nAvailable columns:\n",
      paste(available, collapse = ", ")
    )
  }
  
  available[hit[1]]
}

resolve_survival_columns <- function(
    data,
    endpoint
) {
  endpoint <- toupper(endpoint)
  
  if (endpoint == "TTP") {
    time_candidates <- c(
      "YearsTTP"
    )
    
    event_candidates <- c(
      "CensTTP"
    )
  } else if (endpoint == "TTR") {
    time_candidates <- c(
      "YearsTTR_PostTx1_20260714",
      "YearsTTR_PostTx1",
      "TimeToRel_PostTx1WorkINProgress",
      "yearsTTR_post_tx1"
    )
    
    event_candidates <- c(
      "CensTTR_PostTx1_20260714",
      "CensTTR_PostTx1",
      "CensTimeToRel_PostTx1WorkINProgress",
      "censTTR_post_tx1"
    )
  } else if (endpoint == "PFS") {
    time_candidates <- c(
      "yearsPFS_post_tx1",
      "YearsPFS_PostTx1",
      "YearsPFS_PostTx1_20260714"
    )
    
    event_candidates <- c(
      "censPFS_post_tx1",
      "CensPFS_PostTx1",
      "CensPFS_PostTx1_20260714"
    )
  } else if (endpoint == "OS") {
    time_candidates <- c(
      "yearsOS_post_tx1",
      "YearsOS_PostTx1",
      "YearsOS_PostTx1_20260714"
    )
    
    event_candidates <- c(
      "censOS_post_tx1",
      "CensOS_PostTx1",
      "CensOS_PostTx1_20260714"
    )
  } else {
    stop(
      "ENDPOINT must be 'TTP', 'TTR', 'PFS', or 'OS'."
    )
  }
  
  list(
    time = find_first_column(
      data = data,
      candidates = time_candidates,
      label = paste0(endpoint, " time"),
      required = TRUE
    ),
    event = find_first_column(
      data = data,
      candidates = event_candidates,
      label = paste0(endpoint, " event"),
      required = TRUE
    )
  )
}

normalize_gep1q_class <- function(x) {
  raw <- clean_text(x)
  
  normalized <- tolower(raw)
  normalized <- gsub("\\s+", "", normalized)
  normalized <- gsub(
    "≥",
    ">=",
    normalized,
    fixed = TRUE
  )
  
  output <- rep(
    NA_character_,
    length(normalized)
  )
  
  output[
    normalized %in% c(
      "1q2",
      "1q=2",
      "1q2copy",
      "1q2copies",
      "diploid",
      "2"
    )
  ] <- "1q2"
  
  output[
    normalized %in% c(
      "1q3plus",
      "1q3+",
      "1q>=3",
      "1q>2",
      "1qgain",
      "gain",
      "gained",
      "3plus",
      "3+",
      "3",
      "4"
    )
  ] <- "1q3plus"
  
  output
}

fmt_p <- function(p) {
  if (
    length(p) != 1 ||
    !is.finite(p)
  ) {
    return("NA")
  }
  
  if (p < 0.001) {
    return(
      format(
        p,
        scientific = TRUE,
        digits = 2
      )
    )
  }
  
  formatC(
    p,
    format = "f",
    digits = 3
  )
}

################################################################################
# 4. VALIDATE INPUTS AND CREATE OUTPUT DIRECTORIES
################################################################################

if (!file.exists(ESET_PATH)) {
  stop(
    "Cannot find the processed UAMS ExpressionSet:\n",
    ESET_PATH
  )
}

if (
  length(EVENT_VALUE_ONE_IS_EVENT) != 1 ||
  is.na(EVENT_VALUE_ONE_IS_EVENT) ||
  !is.logical(EVENT_VALUE_ONE_IS_EVENT)
) {
  stop(
    "EVENT_VALUE_ONE_IS_EVENT must be exactly TRUE or FALSE."
  )
}

if (
  length(POLYIG_CUT) != 1 ||
  !is.finite(POLYIG_CUT)
) {
  stop(
    "POLYIG_CUT must be one finite numeric value."
  )
}

create_dir_or_stop(FIGURE_DIR)
create_dir_or_stop(SOURCE_DATA_DIR)

################################################################################
# 5. LOAD EXPRESSIONSET
################################################################################

eset <- readRDS(ESET_PATH)

if (!methods::is(
  eset,
  "ExpressionSet"
)) {
  stop(
    "ESET_PATH did not load an ExpressionSet."
  )
}

phenotype <- as.data.frame(
  Biobase::pData(eset),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

phenotype$PHENOTYPE_ROWNAME <- rownames(phenotype)

df <- tibble::as_tibble(phenotype)

message(
  "Loaded processed UAMS ExpressionSet: ",
  nrow(df),
  " phenotype rows."
)

################################################################################
# 6. RESOLVE REQUIRED COLUMNS
################################################################################

sample_id_column <- find_first_column(
  data = df,
  candidates = c(
    "CHIPID",
    "chipid",
    "SampleID",
    "sample_id",
    "PHENOTYPE_ROWNAME"
  ),
  label = "sample identifier",
  required = TRUE
)

sample_group_column <- find_first_column(
  data = df,
  candidates = c(
    "sample_group",
    "SAMPLE_GROUP"
  ),
  label = "sample_group",
  required = FALSE
)

group_column <- find_first_column(
  data = df,
  candidates = c(
    "GROUP",
    "GROUPHD"
  ),
  label = "NDMM group",
  required = is.null(sample_group_column)
)

gep1q_column <- find_first_column(
  data = df,
  candidates = GEP1Q_COLUMN_CANDIDATES,
  label = "GEP1qcopy_lasso chromosome 1q class",
  required = TRUE
)

if (!identical(
  gep1q_column,
  "GEP1qcopy_lasso"
)) {
  stop(
    "This script must use GEP1qcopy_lasso exactly. Resolved column: ",
    gep1q_column
  )
}

polyig_column <- find_first_column(
  data = df,
  candidates = POLYIG_COLUMN_CANDIDATES,
  label = "PolyIG",
  required = TRUE
)

gep70_column <- find_first_column(
  data = df,
  candidates = GEP70_COLUMN_CANDIDATES,
  label = "GEP70",
  required = TRUE
)

tt_protocol_column <- find_first_column(
  data = df,
  candidates = TT_PROTOCOL_COLUMN_CANDIDATES,
  label = "Total Therapy protocol",
  required = TRUE
)

survival_columns <- resolve_survival_columns(
  data = df,
  endpoint = ENDPOINT
)

time_column <- survival_columns$time
event_column <- survival_columns$event

message(
  "Resolved columns:",
  "\n  Sample ID: ",
  sample_id_column,
  "\n  Sample group: ",
  ifelse(
    is.null(sample_group_column),
    group_column,
    sample_group_column
  ),
  "\n  TT protocol: ",
  tt_protocol_column,
  "\n  GEP1qcopy_lasso: ",
  gep1q_column,
  "\n  PolyIG: ",
  polyig_column,
  "\n  GEP70: ",
  gep70_column,
  "\n  Time: ",
  time_column,
  "\n  Event: ",
  event_column
)

################################################################################
# 7. BUILD THE COMPARISON COHORT
################################################################################

df <- df |>
  dplyr::mutate(
    CHIPID = clean_text(
      .data[[sample_id_column]]
    )
  )

# Restrict to baseline CD138-selected NDMM.
n_before_ndmm <- nrow(df)

if (!is.null(sample_group_column)) {
  df <- df |>
    dplyr::filter(
      clean_text(
        .data[[sample_group_column]]
      ) == "RNAS_CD138_NDMM"
    )
} else {
  df <- df |>
    dplyr::filter(
      clean_text(
        .data[[group_column]]
      ) == "NDMM"
    )
}

message(
  "Restricted to baseline CD138-selected NDMM: ",
  n_before_ndmm,
  " -> ",
  nrow(df)
)

if (nrow(df) == 0) {
  stop(
    "No samples remain after the baseline NDMM restriction."
  )
}

# Restrict to TT2_NoThal and TT2_Thal before the molecular filters.
observed_tt_protocols <- sort(
  unique(
    clean_text(
      df[[tt_protocol_column]]
    )
  )
)

n_before_tt <- nrow(df)

df <- df |>
  dplyr::mutate(
    TT_Protocol = clean_text(
      .data[[tt_protocol_column]]
    )
  ) |>
  dplyr::filter(
    TT_Protocol %in% TT_KEEP
  )

message(
  "Restricted to ",
  paste(TT_KEEP, collapse = " or "),
  ": ",
  n_before_tt,
  " -> ",
  nrow(df)
)

if (nrow(df) == 0) {
  stop(
    "No samples remain after TT2 protocol filtering.\nObserved values:\n",
    paste(observed_tt_protocols, collapse = ", ")
  )
}

# Prepare survival and molecular variables.
df_pre_filter <- df |>
  dplyr::mutate(
    time = to_numeric_safe(
      .data[[time_column]]
    ),
    event_raw = to_numeric_safe(
      .data[[event_column]]
    ),
    event = if (
      EVENT_VALUE_ONE_IS_EVENT
    ) {
      event_raw
    } else {
      1 - event_raw
    },
    GEP1q_raw = clean_text(
      .data[[gep1q_column]]
    ),
    GEP1q_class = normalize_gep1q_class(
      GEP1q_raw
    ),
    PolyIG = to_numeric_safe(
      .data[[polyig_column]]
    ),
    GEP70 = to_numeric_safe(
      .data[[gep70_column]]
    ),
    Exclude_Missing_CHIPID = is.na(CHIPID),
    Exclude_Invalid_Time = (
      !is.finite(time) |
        time < 0
    ),
    Exclude_Invalid_Event = (
      !is.finite(event) |
        !event %in% c(0, 1)
    ),
    Exclude_Invalid_GEP1q = (
      !GEP1q_class %in% c(
        "1q2",
        "1q3plus"
      )
    ),
    Exclude_Invalid_PolyIG = !is.finite(PolyIG),
    Exclude_Invalid_GEP70 = !is.finite(GEP70),
    Keep_Complete_Case = (
      !Exclude_Missing_CHIPID &
        !Exclude_Invalid_Time &
        !Exclude_Invalid_Event &
        !Exclude_Invalid_GEP1q &
        !Exclude_Invalid_PolyIG &
        !Exclude_Invalid_GEP70
    )
  )

complete_case_audit_file <- file.path(
  SOURCE_DATA_DIR,
  paste0(
    FIGURE_ID,
    "_complete-case-audit.csv"
  )
)

readr::write_csv(
  df_pre_filter |>
    dplyr::select(
      CHIPID,
      TT_Protocol,
      dplyr::all_of(gep1q_column),
      GEP1q_raw,
      GEP1q_class,
      PolyIG,
      GEP70,
      time,
      event_raw,
      event,
      Exclude_Missing_CHIPID,
      Exclude_Invalid_Time,
      Exclude_Invalid_Event,
      Exclude_Invalid_GEP1q,
      Exclude_Invalid_PolyIG,
      Exclude_Invalid_GEP70,
      Keep_Complete_Case
    ),
  complete_case_audit_file
)

message(
  "\nComplete-case audit:",
  "\n  Rows entering molecular filter: ",
  nrow(df_pre_filter),
  "\n  Missing CHIPID: ",
  sum(df_pre_filter$Exclude_Missing_CHIPID),
  "\n  Invalid survival time: ",
  sum(df_pre_filter$Exclude_Invalid_Time),
  "\n  Invalid event: ",
  sum(df_pre_filter$Exclude_Invalid_Event),
  "\n  Missing/invalid GEP1qcopy_lasso: ",
  sum(df_pre_filter$Exclude_Invalid_GEP1q),
  "\n  Missing/invalid PolyIG: ",
  sum(df_pre_filter$Exclude_Invalid_PolyIG),
  "\n  Missing/invalid GEP70: ",
  sum(df_pre_filter$Exclude_Invalid_GEP70),
  "\n  Complete cases retained: ",
  sum(df_pre_filter$Keep_Complete_Case)
)

# Keep only the favorable molecular group.
df <- df_pre_filter |>
  dplyr::filter(
    Keep_Complete_Case,
    GEP70 <= GEP70_THRESHOLD,
    GEP1q_class == "1q2",
    PolyIG > POLYIG_CUT
  ) |>
  dplyr::mutate(
    TT_Protocol = factor(
      TT_Protocol,
      levels = c(
        TT_REFERENCE,
        TT_COMPARATOR
      )
    )
  )

if (nrow(df) == 0) {
  stop(
    "No patients remain after requiring GEP70 <= ",
    GEP70_THRESHOLD,
    ", GEP1qcopy_lasso == 1q2, and PolyIG > ",
    POLYIG_CUT,
    "."
  )
}

if (anyDuplicated(df$CHIPID)) {
  duplicated_chipids <- df |>
    dplyr::count(
      CHIPID,
      sort = TRUE
    ) |>
    dplyr::filter(
      n > 1
    )
  
  print(duplicated_chipids)
  
  stop(
    "Duplicate CHIPIDs were found in the active comparison cohort."
  )
}

protocol_counts <- table(df$TT_Protocol)

message(
  "\nFinal favorable-group counts:"
)

print(protocol_counts)

if (
  length(protocol_counts) != 2 ||
  any(protocol_counts < 2)
) {
  stop(
    "Both ",
    TT_REFERENCE,
    " and ",
    TT_COMPARATOR,
    " must contain at least two patients.\nCurrent counts: ",
    paste(
      names(protocol_counts),
      as.integer(protocol_counts),
      sep = "=",
      collapse = "; "
    )
  )
}

protocol_event_counts <- df |>
  dplyr::group_by(
    TT_Protocol
  ) |>
  dplyr::summarise(
    n = dplyr::n(),
    events = sum(event == 1),
    censored = sum(event == 0),
    .groups = "drop"
  )

if (sum(protocol_event_counts$events) == 0) {
  stop(
    "No endpoint events occurred in the final comparison cohort."
  )
}

################################################################################
# 8. SURVIVAL AND COX MODELS
################################################################################

survival_object <- survival::Surv(
  time = df$time,
  event = df$event
)

survival_fit <- survival::survfit(
  survival_object ~ TT_Protocol,
  data = df
)

cox_fit <- survival::coxph(
  survival::Surv(
    time,
    event
  ) ~ TT_Protocol,
  data = df,
  ties = "efron"
)

cox_summary <- summary(cox_fit)

if (nrow(cox_summary$coefficients) != 1) {
  stop(
    "The Cox model did not return exactly one TT2 comparison."
  )
}

cox_hr <- unname(
  exp(
    stats::coef(cox_fit)
  )[1]
)

cox_ci <- unname(
  exp(
    stats::confint(cox_fit)
  )[1, ]
)

cox_p <- unname(
  cox_summary$coefficients[
    1,
    "Pr(>|z|)"
  ]
)

cox_concordance <- unname(
  cox_summary$concordance[1]
)

# Two-group log-rank test.
logrank_fit <- survival::survdiff(
  survival::Surv(
    time,
    event
  ) ~ TT_Protocol,
  data = df,
  rho = 0
)

logrank_chisq <- unname(logrank_fit$chisq)
logrank_df <- length(logrank_fit$n) - 1
logrank_p <- stats::pchisq(
  logrank_chisq,
  df = logrank_df,
  lower.tail = FALSE
)

# Optional proportional-hazards diagnostic.
cox_zph <- survival::cox.zph(cox_fit)

ph_test_p <- if (
  "GLOBAL" %in% rownames(cox_zph$table)
) {
  unname(
    cox_zph$table[
      "GLOBAL",
      "p"
    ]
  )
} else {
  unname(
    cox_zph$table[
      1,
      "p"
    ]
  )
}

comparison_label <- paste0(
  "HR (",
  TT_COMPARATOR,
  " vs ",
  TT_REFERENCE,
  ") = ",
  sprintf("%.2f", cox_hr),
  " (95% CI ",
  sprintf("%.2f", cox_ci[1]),
  "–",
  sprintf("%.2f", cox_ci[2]),
  ")",
  "\nCox Wald p = ",
  fmt_p(cox_p),
  "; log-rank p = ",
  fmt_p(logrank_p)
)

################################################################################
# 9. EXPORT SOURCE DATA AND STATISTICS
################################################################################

analysis_cohort_file <- file.path(
  SOURCE_DATA_DIR,
  paste0(
    FIGURE_ID,
    "_analysis-cohort.csv"
  )
)

group_summary_file <- file.path(
  SOURCE_DATA_DIR,
  paste0(
    FIGURE_ID,
    "_group-summary.csv"
  )
)

comparison_statistics_file <- file.path(
  SOURCE_DATA_DIR,
  paste0(
    FIGURE_ID,
    "_comparison-statistics.csv"
  )
)

run_summary_file <- file.path(
  SOURCE_DATA_DIR,
  paste0(
    FIGURE_ID,
    "_run-summary.csv"
  )
)

analysis_cohort_export <- df |>
  dplyr::transmute(
    CHIPID,
    TT_Protocol = as.character(TT_Protocol),
    GEP1qcopy_lasso_Raw = GEP1q_raw,
    GEP1qcopy_lasso_Group = GEP1q_class,
    PolyIG,
    PolyIG_Group = "High",
    PolyIG_Cutpoint = POLYIG_CUT,
    GEP70,
    GEP70_Threshold = GEP70_THRESHOLD,
    Time_Years = time,
    Event = event,
    Endpoint = ENDPOINT,
    Time_Column = time_column,
    Event_Column = event_column,
    TT_Protocol_Column = tt_protocol_column,
    GEP1qcopy_lasso_Column = gep1q_column,
    PolyIG_Column = polyig_column,
    GEP70_Column = gep70_column
  )

comparison_statistics <- tibble::tibble(
  Reference = TT_REFERENCE,
  Comparator = TT_COMPARATOR,
  n_reference = sum(
    df$TT_Protocol == TT_REFERENCE
  ),
  events_reference = sum(
    df$event[
      df$TT_Protocol == TT_REFERENCE
    ] == 1
  ),
  n_comparator = sum(
    df$TT_Protocol == TT_COMPARATOR
  ),
  events_comparator = sum(
    df$event[
      df$TT_Protocol == TT_COMPARATOR
    ] == 1
  ),
  HR_comparator_vs_reference = cox_hr,
  CI_lower_95 = cox_ci[1],
  CI_upper_95 = cox_ci[2],
  Cox_Wald_p_value = cox_p,
  Logrank_chisquare = logrank_chisq,
  Logrank_df = logrank_df,
  Logrank_p_value = logrank_p,
  Cox_concordance = cox_concordance,
  PH_global_p_value = ph_test_p
)

run_summary <- tibble::tibble(
  Step = STEP_ID,
  Figure = FIGURE_ID,
  Endpoint = ENDPOINT,
  n_total = nrow(df),
  events_total = sum(df$event == 1),
  reference_protocol = TT_REFERENCE,
  comparator_protocol = TT_COMPARATOR,
  gep1q_column = gep1q_column,
  required_gep1q_class = "1q2",
  polyig_column = polyig_column,
  polyig_cutpoint = POLYIG_CUT,
  required_polyig_group = "High",
  gep70_column = gep70_column,
  gep70_threshold = GEP70_THRESHOLD,
  event_value_one_is_event = EVENT_VALUE_ONE_IS_EVENT,
  time_column = time_column,
  event_column = event_column
)

readr::write_csv(
  analysis_cohort_export,
  analysis_cohort_file
)

readr::write_csv(
  protocol_event_counts |>
    dplyr::mutate(
      TT_Protocol = as.character(TT_Protocol)
    ),
  group_summary_file
)

readr::write_csv(
  comparison_statistics,
  comparison_statistics_file
)

readr::write_csv(
  run_summary,
  run_summary_file
)

################################################################################
# 10. BUILD THE TWO-CURVE PLOT
################################################################################

x_max <- ceiling(
  max(
    df$time,
    na.rm = TRUE
  ) / BREAK_TIME_BY_YEARS
) * BREAK_TIME_BY_YEARS

if (
  !is.finite(x_max) ||
  x_max <= 0
) {
  stop(
    "Could not calculate a valid x-axis maximum."
  )
}

plot_title <- switch(
  toupper(ENDPOINT),
  "TTP" = "1q2/PolyIG-High Time to Progression in TT2 With vs Without Thalidomide",
  "TTR" = "1q2/PolyIG-High Outcomes in TT2 With vs Without Thalidomide",
  "PFS" = "1q2/PolyIG-High Progression-Free Survival in TT2",
  "OS" = "1q2/PolyIG-High Overall Survival in TT2"
)

plot_subtitle <- paste0(
  "GEP70 standard risk; n = ",
  nrow(df),
  "; PolyIG > ",
  formatC(
    POLYIG_CUT,
    format = "f",
    digits = 2
  )
)

y_axis_label <- switch(
  toupper(ENDPOINT),
  "TTP" = "Cumulative Progression Probability",
  "TTR" = "Cumulative Relapse Probability",
  "PFS" = "Cumulative Progression/Death Probability",
  "OS" = "Cumulative Death Probability"
)

legend_labels <- c(
  paste0(
    "TT2 No Thalidomide (n = ",
    sum(df$TT_Protocol == TT_REFERENCE),
    ")"
  ),
  paste0(
    "TT2 Thalidomide (n = ",
    sum(df$TT_Protocol == TT_COMPARATOR),
    ")"
  )
)

palette_used <- unname(
  PROTOCOL_COLORS[
    c(
      TT_REFERENCE,
      TT_COMPARATOR
    )
  ]
)

linetype_used <- unname(
  PROTOCOL_LINE_TYPES[
    c(
      TT_REFERENCE,
      TT_COMPARATOR
    )
  ]
)

km_plot <- survminer::ggsurvplot(
  fit = survival_fit,
  data = df,
  fun = "event",
  pval = FALSE,
  risk.table = TRUE,
  risk.table.col = "strata",
  risk.table.y.text = TRUE,
  risk.table.y.text.col = TRUE,
  risk.table.height = 0.22,
  risk.table.title = "No. at risk",
  cumevents = FALSE,
  conf.int = SHOW_CONFIDENCE_INTERVAL,
  censor = SHOW_CENSOR_MARKS,
  palette = palette_used,
  linetype = linetype_used,
  size = 1.25,
  legend = "bottom",
  legend.title = NULL,
  legend.labs = legend_labels,
  break.time.by = BREAK_TIME_BY_YEARS,
  xlim = c(
    0,
    x_max
  ),
  ylim = c(
    0,
    1
  ),
  title = plot_title,
  xlab = "Time from First Transplant (Years)",
  ylab = y_axis_label,
  ggtheme = ggplot2::theme_bw(
    base_size = 16
  ),
  tables.theme = survminer::theme_cleantable()
)

km_plot$plot <- km_plot$plot +
  ggplot2::labs(
    subtitle = plot_subtitle
  ) +
  ggplot2::scale_y_continuous(
    limits = c(
      0,
      1
    ),
    breaks = seq(
      0,
      1,
      by = 0.2
    ),
    labels = function(x) {
      sprintf(
        "%.1f",
        x
      )
    }
  ) +
  ggplot2::annotate(
    geom = "label",
    x = 0.97 * x_max,
    y = 0.97,
    label = comparison_label,
    hjust = 1,
    vjust = 1,
    size = 4.0,
    lineheight = 1.0,
    label.size = 0.25,
    fill = "white"
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    ),
    plot.subtitle = ggplot2::element_text(
      hjust = 0.5,
      size = 15
    ),
    axis.title = ggplot2::element_text(
      size = 16
    ),
    axis.text = ggplot2::element_text(
      size = 13
    ),
    legend.text = ggplot2::element_text(
      size = 13
    ),
    legend.position = "bottom"
  )

km_plot$table <- km_plot$table +
  ggplot2::guides(
    color = "none",
    fill = "none",
    linetype = "none"
  ) +
  ggplot2::theme(
    legend.position = "none",
    axis.title.y = ggplot2::element_blank(),
    axis.text = ggplot2::element_text(
      size = 11
    ),
    plot.margin = ggplot2::margin(
      t = 0,
      r = 5.5,
      b = 5.5,
      l = 5.5,
      unit = "pt"
    )
  )

combined_plot <- ggpubr::ggarrange(
  km_plot$plot,
  km_plot$table,
  ncol = 1,
  heights = c(
    2,
    0.42
  ),
  align = "v"
)

print(combined_plot)

################################################################################
# 11. SAVE FIGURE
################################################################################

png_file <- file.path(
  FIGURE_DIR,
  paste0(
    FIGURE_ID,
    "_",
    ENDPOINT,
    "_KM.png"
  )
)

pdf_file <- file.path(
  FIGURE_DIR,
  paste0(
    FIGURE_ID,
    "_",
    ENDPOINT,
    "_KM.pdf"
  )
)

ggplot2::ggsave(
  filename = png_file,
  plot = combined_plot,
  width = PLOT_WIDTH_IN,
  height = PLOT_HEIGHT_IN,
  units = "in",
  dpi = PNG_DPI,
  bg = "white"
)

ggplot2::ggsave(
  filename = pdf_file,
  plot = combined_plot,
  width = PLOT_WIDTH_IN,
  height = PLOT_HEIGHT_IN,
  units = "in",
  device = "pdf",
  bg = "white"
)

################################################################################
# 12. FINAL SUMMARY
################################################################################

message(
  "\n############################################################",
  "\nStep 25b complete — TT2 thalidomide comparison",
  "\n############################################################",
  "\nFinal cohort: ",
  nrow(df),
  "\nTotal events: ",
  sum(df$event == 1),
  "\n",
  TT_REFERENCE,
  ": n = ",
  sum(df$TT_Protocol == TT_REFERENCE),
  ", events = ",
  sum(
    df$event[
      df$TT_Protocol == TT_REFERENCE
    ] == 1
  ),
  "\n",
  TT_COMPARATOR,
  ": n = ",
  sum(df$TT_Protocol == TT_COMPARATOR),
  ", events = ",
  sum(
    df$event[
      df$TT_Protocol == TT_COMPARATOR
    ] == 1
  ),
  "\nHR (",
  TT_COMPARATOR,
  " vs ",
  TT_REFERENCE,
  "): ",
  sprintf("%.3f", cox_hr),
  "\n95% CI: ",
  sprintf("%.3f", cox_ci[1]),
  " to ",
  sprintf("%.3f", cox_ci[2]),
  "\nCox Wald p-value: ",
  format(cox_p, scientific = TRUE, digits = 4),
  "\nLog-rank chi-square: ",
  sprintf("%.3f", logrank_chisq),
  "\nLog-rank p-value: ",
  format(logrank_p, scientific = TRUE, digits = 4),
  "\nPH global p-value: ",
  format(ph_test_p, scientific = TRUE, digits = 4),
  "\nPNG: ",
  png_file,
  "\nPDF: ",
  pdf_file,
  "\nAnalysis cohort: ",
  analysis_cohort_file,
  "\nGroup summary: ",
  group_summary_file,
  "\nComparison statistics: ",
  comparison_statistics_file,
  "\nRun summary: ",
  run_summary_file,
  "\nComplete-case audit: ",
  complete_case_audit_file,
  "\n"
)
