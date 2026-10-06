################################################################################
# STEP 27 — UAMS FOLLOW-UP TRUNCATION
#
# Uses the frozen Step 26B Figure 3B analysis cohort.
#
# Produces:
#   4-year, 6-year, 8-year, and full-follow-up TTP figures
#   One CSV containing all pairwise HR / 95% CI / p-values at each follow-up
#
# KM plots:
#   - cumulative progression probability
#   - same group colors/order as Step 26B
#   - overall log-rank p-value
#   - No. at risk table
#   - no HR text printed on KM plots
################################################################################

library(Biobase)
library(survival)
library(survminer)
library(dplyr)
library(readr)
library(ggplot2)
library(ggpubr)

################################################################################
# 1. SETTINGS
################################################################################

FIGURE_ID <- "Fig4_Followup"

STEP26B_DIR <- file.path(
  "results", "figures", "Fig3", "Fig3B_Lasso"
)

INPUT_FILE <- file.path(
  STEP26B_DIR,
  "Fig3B_Lasso_analysis-cohort.csv"
)

ESET_PATH <- file.path(
  "data", "processed", "uams",
  "ESET_uams_all_gep1q_fish_polyig_pcme_binaries_lasso_groups_rnbx-matched.rds"
)

OUTPUT_DIR <- file.path(
  "results", "figures", "Fig4"
)

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

FOLLOWUP_HORIZONS <- c(4, 6, 8)

REFERENCE_GROUP <- "CN2/PolyIG High"

GROUP_COLORS <- c(
  "CN2/PolyIG High" = "steelblue2",
  "CN2/PolyIG Low" = "grey50",
  "CN>=3/PolyIG ALL" = "red3",
  "CN>=3/PolyIG Low" = "red2",
  "CN>=3/PolyIG High" = "purple",
  "Intermediate" = "grey50"
)

PREFERRED_GROUP_ORDER <- c(
  "CN2/PolyIG High",
  "CN2/PolyIG Low",
  "CN>=3/PolyIG ALL",
  "CN>=3/PolyIG Low",
  "CN>=3/PolyIG High",
  "Intermediate"
)

BREAK_TIME_BY_YEARS <- 2

SHOW_LOGRANK_P_VALUE <- TRUE
LOGRANK_P_VALUE_COORD <- c(0.5, 0.85)
LOGRANK_P_VALUE_SIZE <- 7

LEGEND_NROW <- 1L
LEGEND_BYROW <- TRUE
LEGEND_TEXT_SIZE <- 13

PLOT_WIDTH_IN <- 9
PLOT_HEIGHT_IN <- 8
PNG_DPI <- 300

################################################################################
# 2. LOAD FROZEN STEP 26B COHORT
################################################################################

if (file.exists(INPUT_FILE)) {
  
  df <- readr::read_csv(
    INPUT_FILE,
    show_col_types = FALSE
  )
  
} else {
  
  message(
    "Step 26B analysis cohort not found. Rebuilding it from the current ESET."
  )
  
  if (!file.exists(ESET_PATH)) {
    stop(
      "Neither the Step 26B cohort nor current UAMS ESET was found.\n",
      "Expected cohort:\n", INPUT_FILE, "\n",
      "Expected ESET:\n", ESET_PATH
    )
  }
  
  eset <- readRDS(ESET_PATH)
  
  if (!methods::is(eset, "ExpressionSet")) {
    stop("ESET_PATH did not contain an ExpressionSet.")
  }
  
  pdat <- as.data.frame(
    Biobase::pData(eset),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  required_columns <- c(
    "chipid",
    "sample_group",
    "TT_On_for_Tx1",
    "GEP1qcopy_lasso",
    "PolyIG_Score",
    "GEP70",
    "YearsTTP",
    "CensTTP"
  )
  
  missing_columns <- setdiff(
    required_columns,
    colnames(pdat)
  )
  
  if (length(missing_columns)) {
    stop(
      "Missing required pData columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  df <- pdat |>
    transmute(
      CHIPID = chipid,
      sample_group = sample_group,
      TT_On_for_Tx1 = TT_On_for_Tx1,
      time = as.numeric(YearsTTP),
      event = as.numeric(CensTTP),
      GEP1q_class = GEP1qcopy_lasso,
      PolyIG = as.numeric(PolyIG_Score),
      GEP70 = as.numeric(GEP70)
    ) |>
    filter(
      sample_group == "RNAS_CD138_NDMM",
      TT_On_for_Tx1 %in% c(
        "TT3a",
        "TT3b",
        "TT4_S-TT3"
      ),
      !is.na(CHIPID),
      is.finite(time),
      time >= 0,
      event %in% c(0, 1),
      GEP1q_class %in% c(
        "1q2",
        "1q3plus"
      ),
      is.finite(PolyIG),
      is.finite(GEP70),
      GEP70 <= 0.66
    ) |>
    mutate(
      PolyIG_group = if_else(
        PolyIG > 11,
        "High",
        "Low"
      ),
      
      Group = case_when(
        GEP1q_class == "1q2" &
          PolyIG_group == "High" ~
          "CN2/PolyIG High",
        
        GEP1q_class == "1q2" &
          PolyIG_group == "Low" ~
          "CN2/PolyIG Low",
        
        GEP1q_class == "1q3plus" ~
          "CN>=3/PolyIG ALL",
        
        TRUE ~
          NA_character_
      )
    ) |>
    filter(
      !is.na(Group)
    )
}

################################################################################
# 3. STANDARDIZE STEP 26B COHORT COLUMNS
################################################################################

required <- c(
  "CHIPID",
  "Group"
)

missing <- setdiff(
  required,
  names(df)
)

if (length(missing)) {
  stop(
    "Step 26B cohort is missing required columns: ",
    paste(missing, collapse = ", ")
  )
}

# Step 26B uses time/event.
if (!all(c("time", "event") %in% names(df))) {
  stop(
    "Step 26B cohort must contain the exact survival columns 'time' and 'event'."
  )
}

df <- df |>
  mutate(
    CHIPID = as.character(CHIPID),
    Group = as.character(Group),
    Time_Years = as.numeric(time),
    Event = as.integer(event)
  )

if (any(!is.finite(df$Time_Years))) {
  stop("Invalid Time_Years values found.")
}

if (any(df$Time_Years < 0)) {
  stop("Negative survival times found.")
}

if (!all(df$Event %in% c(0L, 1L))) {
  stop("Event must contain only 0 and 1.")
}

################################################################################
# 4. GROUP ORDER
################################################################################

GROUP_LEVELS <- PREFERRED_GROUP_ORDER[
  PREFERRED_GROUP_ORDER %in%
    unique(df$Group)
]

df$Group <- factor(
  df$Group,
  levels = GROUP_LEVELS
)

if (nlevels(droplevels(df$Group)) < 2) {
  stop("At least two groups are required.")
}

if (!REFERENCE_GROUP %in% GROUP_LEVELS) {
  stop(
    "Reference group not found: ",
    REFERENCE_GROUP
  )
}

palette_used <- unname(
  GROUP_COLORS[
    GROUP_LEVELS
  ]
)

if (anyNA(palette_used)) {
  stop("A group is missing a defined color.")
}

################################################################################
# 5. ADMINISTRATIVE CENSORING
################################################################################

truncate_data <- function(data, horizon) {
  
  data |>
    mutate(
      Analysis_Time = pmin(
        Time_Years,
        horizon
      ),
      
      Analysis_Event = if_else(
        Event == 1L &
          Time_Years <= horizon,
        1L,
        0L
      )
    )
}

full_data <- function(data) {
  
  data |>
    mutate(
      Analysis_Time = Time_Years,
      Analysis_Event = Event
    )
}

################################################################################
# 6. ALL PAIRWISE COX COMPARISONS
################################################################################

pairwise_cox <- function(data) {
  
  groups <- levels(
    droplevels(
      data$Group
    )
  )
  
  if (length(groups) < 2) {
    return(tibble())
  }
  
  bind_rows(
    lapply(
      combn(
        groups,
        2,
        simplify = FALSE
      ),
      
      function(pair) {
        
        reference <- pair[1]
        comparison <- pair[2]
        
        d <- data |>
          filter(
            Group %in% pair
          ) |>
          mutate(
            Group_pair = factor(
              as.character(Group),
              levels = pair
            )
          )
        
        n_reference <- sum(
          d$Group_pair == reference
        )
        
        n_comparison <- sum(
          d$Group_pair == comparison
        )
        
        events_reference <- sum(
          d$Analysis_Event[
            d$Group_pair == reference
          ] == 1
        )
        
        events_comparison <- sum(
          d$Analysis_Event[
            d$Group_pair == comparison
          ] == 1
        )
        
        empty <- tibble(
          Reference = reference,
          Comparison = comparison,
          n_reference = n_reference,
          n_comparison = n_comparison,
          events_reference = events_reference,
          events_comparison = events_comparison,
          HR = NA_real_,
          CI_lower = NA_real_,
          CI_upper = NA_real_,
          P_value = NA_real_
        )
        
        if (
          n_reference < 2 ||
          n_comparison < 2 ||
          events_reference + events_comparison == 0
        ) {
          return(empty)
        }
        
        fit <- tryCatch(
          coxph(
            Surv(
              Analysis_Time,
              Analysis_Event
            ) ~ Group_pair,
            data = d
          ),
          error = function(e) NULL
        )
        
        if (is.null(fit)) {
          return(empty)
        }
        
        s <- summary(fit)
        ci <- exp(confint(fit))
        
        tibble(
          Reference = reference,
          Comparison = comparison,
          n_reference = n_reference,
          n_comparison = n_comparison,
          events_reference = events_reference,
          events_comparison = events_comparison,
          HR = unname(
            exp(coef(fit)[1])
          ),
          CI_lower = ci[1, 1],
          CI_upper = ci[1, 2],
          P_value = s$coefficients[
            1,
            "Pr(>|z|)"
          ]
        )
      }
    )
  )
}

################################################################################
# 7. KM PLOT
################################################################################

make_km <- function(
    data,
    label,
    x_max,
    file_label
) {
  
  fit <- survfit(
    Surv(
      Analysis_Time,
      Analysis_Event
    ) ~ Group,
    data = data
  )
  
  logrank <- survdiff(
    Surv(
      Analysis_Time,
      Analysis_Event
    ) ~ Group,
    data = data
  )
  
  logrank_p <- pchisq(
    logrank$chisq,
    df = length(logrank$n) - 1,
    lower.tail = FALSE
  )
  
  title_text <- paste0(
    "Time to Progression by GEP 1q CN × PolyIG (PI/IMiD)\n",
    "Administratively Censored at ",
    label
  )
  
  subtitle_text <- paste0(
    "Standard Risk (GEP70); TT3a, TT3b, TT4-S (n = ",
    nrow(data),
    ")"
  )
  
  km <- survminer::ggsurvplot(
    fit = fit,
    data = data,
    fun = "event",
    
    pval = SHOW_LOGRANK_P_VALUE,
    pval.coord = LOGRANK_P_VALUE_COORD,
    pval.size = LOGRANK_P_VALUE_SIZE,
    
    risk.table = TRUE,
    risk.table.col = "strata",
    risk.table.y.text = TRUE,
    risk.table.y.text.col = TRUE,
    risk.table.height = 0.30,
    risk.table.title = "No. at risk",
    
    cumevents = FALSE,
    conf.int = FALSE,
    censor = TRUE,
    
    palette = palette_used,
    
    legend = "bottom",
    legend.title = "Group",
    legend.labs = GROUP_LEVELS,
    
    break.time.by = BREAK_TIME_BY_YEARS,
    
    xlim = c(0, x_max),
    ylim = c(0, 1),
    
    title = title_text,
    subtitle = subtitle_text,
    
    xlab = "Time from First Transplant (years)",
    ylab = "Cumulative Progression Probability",
    
    ggtheme = theme_bw(
      base_size = 16
    ) +
      theme(
        plot.title = element_text(
          hjust = 0.5,
          face = "bold"
        ),
        plot.subtitle = element_text(
          hjust = 0.5,
          size = 18
        )
      ),
    
    tables.theme =
      survminer::theme_cleantable()
  )
  
  km$plot <- km$plot +
    guides(
      color = guide_legend(
        nrow = LEGEND_NROW,
        byrow = LEGEND_BYROW
      )
    ) +
    theme(
      legend.text = element_text(
        size = LEGEND_TEXT_SIZE
      )
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(
        0,
        1,
        by = 0.2
      ),
      labels = function(x) {
        sprintf("%.1f", x)
      }
    )
  
  # No HR annotation on the KM plot.
  km$table <- km$table +
    guides(
      color = "none",
      fill = "none"
    ) +
    theme(
      legend.position = "none",
      axis.title.y = element_blank(),
      plot.margin = margin(
        t = 0,
        r = 5.5,
        b = 5.5,
        l = 5.5,
        unit = "pt"
      )
    )
  
  combined <- ggarrange(
    km$plot,
    km$table,
    ncol = 1,
    heights = c(
      2,
      0.6
    ),
    align = "v"
  )
  
  print(combined)
  
  png_file <- file.path(
    OUTPUT_DIR,
    paste0(
      FIGURE_ID,
      "_",
      file_label,
      "_TTP_KM.png"
    )
  )
  
  pdf_file <- file.path(
    OUTPUT_DIR,
    paste0(
      FIGURE_ID,
      "_",
      file_label,
      "_TTP_KM.pdf"
    )
  )
  
  ggsave(
    png_file,
    combined,
    width = PLOT_WIDTH_IN,
    height = PLOT_HEIGHT_IN,
    units = "in",
    dpi = PNG_DPI,
    bg = "white"
  )
  
  ggsave(
    pdf_file,
    combined,
    width = PLOT_WIDTH_IN,
    height = PLOT_HEIGHT_IN,
    units = "in",
    device = "pdf",
    bg = "white"
  )
  
  list(
    fit = fit,
    logrank_p = logrank_p,
    png = png_file,
    pdf = pdf_file
  )
}

################################################################################
# 8. RUN 4-, 6-, AND 8-YEAR ANALYSES
################################################################################

results <- lapply(
  FOLLOWUP_HORIZONS,
  function(horizon) {
    
    make_km(
      truncate_data(
        df,
        horizon
      ),
      paste0(
        horizon,
        " Years"
      ),
      horizon,
      paste0(
        horizon,
        "yr"
      )
    )
  }
)

################################################################################
# 9. FULL FOLLOW-UP
################################################################################

full_xmax <- ceiling(
  max(
    df$Time_Years,
    na.rm = TRUE
  ) /
    BREAK_TIME_BY_YEARS
) * BREAK_TIME_BY_YEARS

full_result <- make_km(
  full_data(df),
  "Full Follow-up",
  full_xmax,
  "Full"
)

################################################################################
# 10. ALL PAIRWISE HR / P-VALUE RESULTS
################################################################################

pairwise_results <- bind_rows(
  
  lapply(
    FOLLOWUP_HORIZONS,
    function(horizon) {
      
      x <- pairwise_cox(
        truncate_data(
          df,
          horizon
        )
      )
      
      x$Followup <- paste0(
        horizon,
        "-year"
      )
      
      x
    }
  ),
  
  {
    x <- pairwise_cox(
      full_data(df)
    )
    
    x$Followup <-
      "Full follow-up"
    
    x
  }
)

pairwise_csv <- file.path(
  OUTPUT_DIR,
  "Fig4_Followup_Pairwise_HR_Pvalues_TTP.csv"
)

pairwise_export <- pairwise_results |>
  mutate(
    Comparison_Label = paste0(
      Comparison,
      " vs ",
      Reference
    ),
    
    `HR (95% CI)` = ifelse(
      is.finite(HR) &
        is.finite(CI_lower) &
        is.finite(CI_upper),
      
      sprintf(
        "%.2f (%.2f–%.2f)",
        HR,
        CI_lower,
        CI_upper
      ),
      
      "NA"
    ),
    
    `p-value` = ifelse(
      is.finite(P_value),
      
      ifelse(
        P_value < 0.001,
        "<0.001",
        sprintf(
          "%.3f",
          P_value
        )
      ),
      
      "NA"
    )
  ) |>
  select(
    Followup,
    Comparison = Comparison_Label,
    `HR (95% CI)`,
    `p-value`,
    n_reference,
    n_comparison,
    events_reference,
    events_comparison
  )

readr::write_csv(
  pairwise_export,
  pairwise_csv
)

################################################################################
# 11. FOLLOW-UP SUMMARY
################################################################################

followup_summary <- bind_rows(
  lapply(
    seq_along(FOLLOWUP_HORIZONS),
    function(i) {
      
      horizon <- FOLLOWUP_HORIZONS[i]
      
      d <- truncate_data(
        df,
        horizon
      )
      
      tibble(
        Followup = paste0(horizon, "-year"),
        N = nrow(d),
        Events = sum(d$Analysis_Event == 1),
        At_Risk_At_Horizon = sum(df$Time_Years >= horizon),
        Logrank_P = results[[i]]$logrank_p
      )
    }
  ),
  
  tibble(
    Followup = "Full follow-up",
    N = nrow(df),
    Events = sum(df$Event == 1),
    At_Risk_At_Horizon = NA_integer_,
    Logrank_P = full_result$logrank_p
  )
)

followup_summary_file <- file.path(
  OUTPUT_DIR,
  "Fig4_Followup_Summary_TTP.csv"
)

readr::write_csv(
  followup_summary,
  followup_summary_file
)
################################################################################
# 12. FINAL CHECK
################################################################################

cat(
  "\nStep 27 complete.\n",
  "Endpoint = TTP\n",
  "Frozen Step 26B cohort n = ", nrow(df), "\n",
  "Events = ", sum(df$Event == 1), "\n",
  "Groups:\n"
)

print(
  table(
    df$Group
  )
)

cat(
  "\nPairwise HR/p-value CSV:\n",
  pairwise_csv,
  "\n\nFollow-up summary CSV:\n",
  followup_summary_file,
  "\n\nOutput directory:\n",
  OUTPUT_DIR,
  "\n"
)