################################################################################
# STEP 33 — MULTIVARIABLE COX REGRESSION
#
# Input:
#   Step 32 TTP analysis cohort
#
# Endpoint:
#   TTP from first transplant
#   event = 1, censored = 0
#
# Model:
#   Manually specified predictors only
#   No automatic or stepwise variable selection
#   One common complete-case cohort
#
# Output:
#   results/table1/Multivariate/
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(survival)
  library(tibble)
})


################################################################################
# 1. FILES / SETTINGS
################################################################################

INPUT_FILE <- file.path(
  "results",
  "table1",
  "Univariate",
  "Step32_UnivariateCox_analysis-cohort_TTP.csv"
)

OUT_DIR <- file.path(
  "results",
  "table1",
  "Multivariate"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

TIES_METHOD <- "efron"

PH_TRANSFORM <- "km"

MIN_N <- 20

MIN_EVENTS <- 5

MIN_EVENTS_PER_PARAMETER_WARNING <- 10


################################################################################
# 2. MANUALLY SPECIFIED MULTIVARIABLE MODEL
#
# Edit this block only if you want to change the model.
################################################################################

VARIABLES_TO_INCLUDE <- c(
  "Age_Great_65",
  "ALBUMIN",
  "HGB116F132M",
  "PLATELETS",
  "Cyto_MM",
  "Tandem_Tx2",
  "HY",
  "t14_16_20",
  "PR",
  "PolyIG_high_bin",
  "GEP1q_3plus_bin"
)


# All variables in the current model are binary 0/1.
BINARY_VARS <- c(
  "Age_Great_65",
  "ALBUMIN",
  "HGB116F132M",
  "PLATELETS",
  "Cyto_MM",
  "Tandem_Tx2",
  "HY",
  "t14_16_20",
  "PR",
  "PolyIG_high_bin",
  "GEP1q_3plus_bin"
)


# No continuous predictors are currently included in this model.
CONTINUOUS_VARS <- character(0)


################################################################################
# 3. READ STEP 32 ANALYSIS COHORT
################################################################################

if (!file.exists(INPUT_FILE)) {
  stop(
    "Cannot find Step 32 TTP analysis cohort:\n",
    INPUT_FILE,
    "\nRun Step 32 first."
  )
}

df <- readr::read_csv(
  INPUT_FILE,
  show_col_types = FALSE
)

required <- c(
  "PATID",
  "time",
  "event",
  VARIABLES_TO_INCLUDE
)

missing_vars <- setdiff(
  required,
  names(df)
)

if (length(missing_vars)) {
  stop(
    "Missing required variable(s):\n",
    paste(
      missing_vars,
      collapse = ", "
    )
  )
}

if (anyDuplicated(df$PATID)) {
  stop(
    "Duplicate PATIDs detected in Step 32 analysis cohort."
  )
}


################################################################################
# 4. PREPARE MODEL VARIABLES
################################################################################

num <- function(x) {
  suppressWarnings(
    as.numeric(
      as.character(x)
    )
  )
}


model_data <- df |>
  select(
    PATID,
    time,
    event,
    all_of(VARIABLES_TO_INCLUDE)
  ) |>
  mutate(
    time = num(time),
    event = num(event)
  )


################################################################################
# 5. FORCE BINARY VARIABLES TO 0/1 FACTORS
#
# 0 = reference
# 1 = comparison
#
# Stop immediately if a supposed binary variable contains anything other than
# 0, 1, or NA.
################################################################################

for (v in intersect(BINARY_VARS, VARIABLES_TO_INCLUDE)) {
  
  x <- num(
    model_data[[v]]
  )
  
  invalid <- unique(
    x[
      !is.na(x) &
        !x %in% c(0, 1)
    ]
  )
  
  if (length(invalid)) {
    stop(
      v,
      " was specified as binary but contains value(s): ",
      paste(
        invalid,
        collapse = ", "
      )
    )
  }
  
  model_data[[v]] <- factor(
    x,
    levels = c(
      0,
      1
    )
  )
}


################################################################################
# 6. CONTINUOUS VARIABLES
################################################################################

for (v in intersect(CONTINUOUS_VARS, VARIABLES_TO_INCLUDE)) {
  
  model_data[[v]] <- num(
    model_data[[v]]
  )
}


################################################################################
# 7. HANDLE ANY OTHER VARIABLES
#
# This section is retained for flexibility if future models include variables
# not explicitly listed in BINARY_VARS or CONTINUOUS_VARS.
#
# numeric -> continuous
# character/factor -> categorical
################################################################################

other_vars <- setdiff(
  VARIABLES_TO_INCLUDE,
  c(
    BINARY_VARS,
    CONTINUOUS_VARS
  )
)

for (v in other_vars) {
  
  x_num <- num(
    model_data[[v]]
  )
  
  x_original <- model_data[[v]]
  
  nonmissing <- sum(
    !is.na(
      x_original
    )
  )
  
  numeric_fraction <- if (nonmissing > 0) {
    
    sum(
      !is.na(
        x_num
      )
    ) / nonmissing
    
  } else {
    
    0
  }
  
  if (numeric_fraction >= 0.90) {
    
    model_data[[v]] <- x_num
    
  } else {
    
    model_data[[v]] <- factor(
      as.character(
        x_original
      )
    )
  }
}


################################################################################
# 8. COMMON COMPLETE-CASE COHORT
################################################################################

valid_endpoint <- (
  is.finite(model_data$time) &
    model_data$time >= 0 &
    model_data$event %in% c(0, 1)
)

complete_predictors <- complete.cases(
  model_data[
    ,
    VARIABLES_TO_INCLUDE,
    drop = FALSE
  ]
)

included <- valid_endpoint &
  complete_predictors


row_audit <- tibble(
  PATID = model_data$PATID,
  Valid_endpoint = valid_endpoint,
  Complete_all_predictors = complete_predictors,
  Included_in_model = included
)


model_complete <- model_data[
  included,
  ,
  drop = FALSE
]


# Remove unused factor levels after complete-case filtering.
for (v in VARIABLES_TO_INCLUDE) {
  
  if (is.factor(model_complete[[v]])) {
    
    model_complete[[v]] <- droplevels(
      model_complete[[v]]
    )
  }
}


N_MODEL <- nrow(
  model_complete
)

N_EVENTS <- sum(
  model_complete$event == 1
)

N_CENSORED <- sum(
  model_complete$event == 0
)


if (N_MODEL < MIN_N) {
  
  stop(
    "Only ",
    N_MODEL,
    " complete cases; minimum is ",
    MIN_N,
    "."
  )
}


if (N_EVENTS < MIN_EVENTS) {
  
  stop(
    "Only ",
    N_EVENTS,
    " progression events; minimum is ",
    MIN_EVENTS,
    "."
  )
}


################################################################################
# 9. COUNT MODEL PARAMETERS
################################################################################

parameters_per_variable <- sapply(
  VARIABLES_TO_INCLUDE,
  function(v) {
    
    x <- model_complete[[v]]
    
    if (is.factor(x)) {
      
      if (nlevels(x) < 2) {
        
        stop(
          v,
          " has fewer than two levels in the complete-case cohort."
        )
      }
      
      nlevels(x) - 1
      
    } else {
      
      if (length(unique(x)) < 2) {
        
        stop(
          v,
          " has no variation in the complete-case cohort."
        )
      }
      
      1
    }
  }
)


N_PARAMETERS <- sum(
  parameters_per_variable
)

EVENTS_PER_PARAMETER <- N_EVENTS /
  N_PARAMETERS


if (N_EVENTS <= N_PARAMETERS) {
  
  stop(
    "Events (",
    N_EVENTS,
    ") are not greater than fitted parameters (",
    N_PARAMETERS,
    ")."
  )
}


if (
  EVENTS_PER_PARAMETER <
  MIN_EVENTS_PER_PARAMETER_WARNING
) {
  
  warning(
    "Events per parameter = ",
    round(
      EVENTS_PER_PARAMETER,
      2
    ),
    ". This is below ",
    MIN_EVENTS_PER_PARAMETER_WARNING,
    "; interpret adjusted estimates cautiously."
  )
}


################################################################################
# 10. FIT MULTIVARIABLE COX MODEL
################################################################################

cox_formula <- reformulate(
  VARIABLES_TO_INCLUDE,
  response = "Surv(time, event)"
)


cox_fit <- coxph(
  cox_formula,
  data = model_complete,
  ties = TIES_METHOD,
  x = TRUE,
  y = TRUE,
  model = TRUE,
  na.action = na.fail
)


cox_summary <- summary(
  cox_fit
)


################################################################################
# 11. ADJUSTED HR TABLE
################################################################################

coef_table <- as.data.frame(
  cox_summary$coefficients
)

ci_table <- as.data.frame(
  cox_summary$conf.int
)


results <- tibble(
  Term = rownames(coef_table),
  HR = ci_table[["exp(coef)"]],
  CI_low = ci_table[["lower .95"]],
  CI_high = ci_table[["upper .95"]],
  p = coef_table[["Pr(>|z|)"]]
) |>
  mutate(
    HR_95CI = sprintf(
      "%.2f (%.2f–%.2f)",
      HR,
      CI_low,
      CI_high
    ),
    N = N_MODEL,
    Events = N_EVENTS
  )


################################################################################
# 12. MAP MODEL TERMS BACK TO ORIGINAL VARIABLES
################################################################################

get_variable <- function(term) {
  
  hits <- VARIABLES_TO_INCLUDE[
    startsWith(
      term,
      VARIABLES_TO_INCLUDE
    )
  ]
  
  if (!length(hits)) {
    
    return(
      term
    )
  }
  
  hits[
    which.max(
      nchar(
        hits
      )
    )
  ]
}


results <- results |>
  mutate(
    Variable = vapply(
      Term,
      get_variable,
      character(1)
    ),
    
    Contrast = case_when(
      Variable %in% BINARY_VARS ~
        "1 vs 0",
      
      Variable %in% CONTINUOUS_VARS ~
        "Per 1-unit increase",
      
      TRUE ~
        "Categorical comparison"
    )
  ) |>
  select(
    Variable,
    Contrast,
    Term,
    N,
    Events,
    HR,
    CI_low,
    CI_high,
    HR_95CI,
    p
  )


################################################################################
# 13. GLOBAL TERM TESTS
#
# Useful especially if future models contain categorical variables with
# more than two levels.
################################################################################

drop1_result <- drop1(
  cox_fit,
  test = "Chisq"
)


term_tests <- as.data.frame(
  drop1_result
) |>
  rownames_to_column(
    "Variable"
  ) |>
  filter(
    Variable != "<none>"
  )


################################################################################
# 14. PROPORTIONAL-HAZARDS TEST
################################################################################

ph_test <- cox.zph(
  cox_fit,
  transform = PH_TRANSFORM
)


ph_table <- as.data.frame(
  ph_test$table
) |>
  rownames_to_column(
    "Term"
  )


################################################################################
# 15. MODEL SUMMARY
################################################################################

model_summary <- tibble(
  Metric = c(
    "Eligible Step 32 TTP cohort",
    "Common complete-case N",
    "Progression events",
    "Censored",
    "Variables",
    "Model parameters",
    "Events per parameter",
    "Concordance",
    "Likelihood-ratio p",
    "Wald p",
    "Score p"
  ),
  
  Value = c(
    nrow(df),
    N_MODEL,
    N_EVENTS,
    N_CENSORED,
    length(VARIABLES_TO_INCLUDE),
    N_PARAMETERS,
    EVENTS_PER_PARAMETER,
    unname(
      cox_summary$concordance[1]
    ),
    unname(
      cox_summary$logtest["pvalue"]
    ),
    unname(
      cox_summary$waldtest["pvalue"]
    ),
    unname(
      cox_summary$sctest["pvalue"]
    )
  )
)


################################################################################
# 16. VARIABLE AUDIT
################################################################################

variable_audit <- tibble(
  Variable = VARIABLES_TO_INCLUDE,
  
  Type = case_when(
    VARIABLES_TO_INCLUDE %in%
      BINARY_VARS ~
      "binary",
    
    VARIABLES_TO_INCLUDE %in%
      CONTINUOUS_VARS ~
      "continuous",
    
    TRUE ~
      "auto"
  ),
  
  Parameters = as.integer(
    parameters_per_variable
  ),
  
  N_nonmissing_before_common_filter = sapply(
    VARIABLES_TO_INCLUDE,
    function(v) {
      
      sum(
        !is.na(
          model_data[[v]]
        )
      )
    }
  ),
  
  N_common_model = N_MODEL
)


################################################################################
# 17. DESIGN-MATRIX CORRELATION
################################################################################

design_matrix <- model.matrix(
  cox_fit
)

design_matrix <- design_matrix[
  ,
  colnames(design_matrix) != "(Intercept)",
  drop = FALSE
]


design_correlation <- if (
  ncol(design_matrix) > 0
) {
  
  suppressWarnings(
    cor(
      design_matrix
    )
  )
  
} else {
  
  matrix(
    numeric(0),
    0,
    0
  )
}


################################################################################
# 18. OUTPUT FILES
################################################################################

HR_FILE <- file.path(
  OUT_DIR,
  "Step33_MultivariableCox_TTP_adjusted_HR.csv"
)

SUMMARY_FILE <- file.path(
  OUT_DIR,
  "Step33_MultivariableCox_TTP_model_summary.csv"
)

TERM_FILE <- file.path(
  OUT_DIR,
  "Step33_MultivariableCox_TTP_term_tests.csv"
)

PH_FILE <- file.path(
  OUT_DIR,
  "Step33_MultivariableCox_TTP_proportional_hazards.csv"
)

VARIABLE_AUDIT_FILE <- file.path(
  OUT_DIR,
  "Step33_MultivariableCox_TTP_variable_audit.csv"
)

ROW_AUDIT_FILE <- file.path(
  OUT_DIR,
  "Step33_MultivariableCox_TTP_row_audit.csv"
)

COMPLETE_CASE_FILE <- file.path(
  OUT_DIR,
  "Step33_MultivariableCox_TTP_complete_case_cohort.csv"
)

CORRELATION_FILE <- file.path(
  OUT_DIR,
  "Step33_MultivariableCox_TTP_design_matrix_correlation.csv"
)

MODEL_FILE <- file.path(
  OUT_DIR,
  "Step33_MultivariableCox_TTP_model.rds"
)

PH_PLOT_FILE <- file.path(
  OUT_DIR,
  "Step33_MultivariableCox_TTP_Schoenfeld_plots.pdf"
)


################################################################################
# 19. SAVE
################################################################################

readr::write_excel_csv(
  results,
  HR_FILE
)

readr::write_excel_csv(
  model_summary,
  SUMMARY_FILE
)

readr::write_excel_csv(
  term_tests,
  TERM_FILE
)

readr::write_excel_csv(
  ph_table,
  PH_FILE
)

readr::write_excel_csv(
  variable_audit,
  VARIABLE_AUDIT_FILE
)

readr::write_excel_csv(
  row_audit,
  ROW_AUDIT_FILE
)

readr::write_csv(
  model_complete,
  COMPLETE_CASE_FILE
)

write.csv(
  design_correlation,
  CORRELATION_FILE,
  row.names = TRUE,
  na = ""
)

saveRDS(
  cox_fit,
  MODEL_FILE
)


################################################################################
# 20. SCHOENFELD RESIDUAL PLOTS
################################################################################

pdf(
  PH_PLOT_FILE,
  width = 8.5,
  height = 7,
  onefile = TRUE
)

plot(
  ph_test
)

dev.off()


################################################################################
# 21. FINAL SUMMARY
################################################################################

cat(
  "\n============================================================\n",
  "STEP 33 — MULTIVARIABLE COX COMPLETE\n",
  "============================================================\n",
  "Endpoint: TTP\n",
  "Eligible Step 32 cohort: ",
  nrow(df),
  "\n",
  "Complete-case model N: ",
  N_MODEL,
  "\n",
  "Progression events: ",
  N_EVENTS,
  "\n",
  "Censored: ",
  N_CENSORED,
  "\n",
  "Variables: ",
  length(VARIABLES_TO_INCLUDE),
  "\n",
  "Parameters: ",
  N_PARAMETERS,
  "\n",
  "Events/parameter: ",
  round(
    EVENTS_PER_PARAMETER,
    2
  ),
  "\n",
  sep = ""
)


cat(
  "\nVariables in model:\n"
)

print(
  VARIABLES_TO_INCLUDE
)


cat(
  "\nAdjusted hazard ratios:\n"
)

print(
  results |>
    select(
      Variable,
      Contrast,
      N,
      Events,
      HR_95CI,
      p
    )
)


cat(
  "\nProportional-hazards test:\n"
)

print(
  ph_test
)


cat(
  "\nFiles saved under:\n",
  OUT_DIR,
  "\n",
  sep = ""
)