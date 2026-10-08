################################################################################
# STEP 29 — FIGURE 5B
# TTP BY FISH 1q COPY NUMBER × unIGH ARD
#
# Cohort:
#   1. Baseline RNAS_CD138_NDMM
#   2. TT3a, TT3b, TT4_S-TT3
#   3. GEP70 <= 0.66
#   4. Valid FISH_1q_num
#   5. Valid IG_ARD
#   6. Valid YearsTTP / CensTTP
#
# Groups:
#   CN2/unIGH High
#   CN2/unIGH Low
#   CN≥3/unIGH ALL
################################################################################

suppressPackageStartupMessages({
  library(Biobase)
  library(survival)
  library(survminer)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggpubr)
  library(tibble)
})

################################################################################
# 0. VERIFY PIPELINE ROOTS
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

################################################################################
# 1. CONFIGURATION
################################################################################

FIGURE_ID <- "Fig5B"
ENDPOINT <- "TTP"

ESET_PATH <- file.path(
  DATA_ROOT, "processed", "ESET_uams_all_gep1q_fish_polyig_pcme_binaries_lasso_groups_rnbx-matched.rds"
)

FIGURE_DIR <- file.path(
  RESULTS_ROOT, "figures", "Fig5", FIGURE_ID
)

dir.create(
  FIGURE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# Exact columns
SAMPLE_ID_COL <- "chipid"
SAMPLE_GROUP_COL <- "sample_group"
TT_COL <- "TT_On_for_Tx1"

FISH1Q_COL <- "FISH_1q_num"
UNIGH_COL <- "IG_ARD"
GEP70_COL <- "GEP70"

TIME_COL <- "YearsTTP"
EVENT_COL <- "CensTTP"

UNIGH_CUT <- 0.29
GEP70_CUT <- 0.66

TT_KEEP <- c(
  "TT3a",
  "TT3b",
  "TT4_S-TT3"
)

################################################################################
# 2. PLOT SETTINGS
################################################################################

GROUP_CN2_HIGH <- "CN2/unIGH High"
GROUP_CN2_LOW <- "CN2/unIGH Low"
GROUP_CN3_ALL <- "CN≥3/unIGH ALL"

GROUP_LEVELS <- c(
  GROUP_CN2_HIGH,
  GROUP_CN2_LOW,
  GROUP_CN3_ALL
)

GROUP_COLORS <- c(
  "CN2/unIGH High" = "steelblue2",
  "CN2/unIGH Low" = "grey50",
  "CN≥3/unIGH ALL" = "red3"
)

SHOW_LOGRANK_P_VALUE <- TRUE
LOGRANK_P_VALUE_COORD <- c(0.5, 0.85)
LOGRANK_P_VALUE_SIZE <- 7

SHOW_HR_LABEL <- TRUE
HR_LABEL_TEXT_SIZE <- 3.6
HR_REFERENCE_GROUP <- GROUP_CN2_HIGH

LEGEND_NROW <- 1L
LEGEND_BYROW <- TRUE
LEGEND_TEXT_SIZE <- 13

BREAK_TIME_BY_YEARS <- 2

PLOT_WIDTH_IN <- 9
PLOT_HEIGHT_IN <- 8
PNG_DPI <- 300

################################################################################
# 3. HELPERS
################################################################################

fmt_p <- function(p) {
  if (!is.finite(p)) return("NA")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

fmt_hr <- function(hr, lo, hi) {
  if (any(!is.finite(c(hr, lo, hi)))) return("NA")
  sprintf("%.2f (%.2f–%.2f)", hr, lo, hi)
}

################################################################################
# 4. LOAD ESET
################################################################################

if (!file.exists(ESET_PATH)) {
  stop(
    "Cannot find UAMS ExpressionSet:\n",
    ESET_PATH
  )
}

eset <- readRDS(
  ESET_PATH
)

if (!methods::is(eset, "ExpressionSet")) {
  stop("ESET_PATH did not contain an ExpressionSet.")
}

pdat <- as.data.frame(
  Biobase::pData(eset),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

################################################################################
# 5. CHECK EXACT REQUIRED COLUMNS
################################################################################

required_columns <- c(
  SAMPLE_ID_COL,
  SAMPLE_GROUP_COL,
  TT_COL,
  FISH1Q_COL,
  UNIGH_COL,
  GEP70_COL,
  TIME_COL,
  EVENT_COL
)

missing_columns <- setdiff(
  required_columns,
  colnames(pdat)
)

if (length(missing_columns)) {
  stop(
    "Missing required pData columns: ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}

################################################################################
# 6. BUILD COHORT
################################################################################

df <- pdat |>
  tibble::as_tibble() |>
  mutate(
    CHIPID = .data[[SAMPLE_ID_COL]],
    sample_group = .data[[SAMPLE_GROUP_COL]],
    TT = .data[[TT_COL]],
    FISH1q = as.numeric(.data[[FISH1Q_COL]]),
    unIGH = as.numeric(.data[[UNIGH_COL]]),
    GEP70 = as.numeric(.data[[GEP70_COL]]),
    time = as.numeric(.data[[TIME_COL]]),
    event = as.numeric(.data[[EVENT_COL]])
  ) |>
  filter(
    sample_group == "RNAS_CD138_NDMM",
    TT %in% TT_KEEP,
    
    is.finite(FISH1q),
    FISH1q == 2 | FISH1q >= 3,
    
    is.finite(unIGH),
    
    is.finite(GEP70),
    GEP70 <= GEP70_CUT,
    
    is.finite(time),
    time >= 0,
    event %in% c(0, 1),
    
    !is.na(CHIPID)
  ) |>
  mutate(
    FISH1q_group = case_when(
      FISH1q == 2 ~ "1q2",
      FISH1q >= 3 ~ "1q>=3",
      TRUE ~ NA_character_
    ),
    
    unIGH_group = if_else(
      unIGH > UNIGH_CUT,
      "High",
      "Low"
    ),
    
    Group = case_when(
      FISH1q_group == "1q2" &
        unIGH_group == "High" ~
        GROUP_CN2_HIGH,
      
      FISH1q_group == "1q2" &
        unIGH_group == "Low" ~
        GROUP_CN2_LOW,
      
      FISH1q_group == "1q>=3" ~
        GROUP_CN3_ALL,
      
      TRUE ~ NA_character_
    ),
    
    Group = factor(
      Group,
      levels = GROUP_LEVELS
    )
  ) |>
  filter(
    !is.na(Group)
  )

df$Group <- droplevels(
  df$Group
)

if (!nrow(df)) {
  stop("No samples remain in the Figure 5B cohort.")
}

if (anyDuplicated(df$CHIPID)) {
  stop("Duplicate CHIPIDs are present in the Figure 5B cohort.")
}

if (nlevels(df$Group) < 2) {
  stop("At least two groups are required.")
}

if (!HR_REFERENCE_GROUP %in% levels(df$Group)) {
  stop(
    "HR reference group is absent: ",
    HR_REFERENCE_GROUP
  )
}

################################################################################
# 7. SURVIVAL MODEL
################################################################################

surv_fit <- survfit(
  Surv(time, event) ~ Group,
  data = df
)

logrank <- survdiff(
  Surv(time, event) ~ Group,
  data = df
)

overall_p <- pchisq(
  logrank$chisq,
  df = length(logrank$n) - 1,
  lower.tail = FALSE
)

################################################################################
# 8. PAIRWISE COX HR / P VALUES
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
  
  pairs <- combn(
    groups,
    2,
    simplify = FALSE
  )
  
  bind_rows(
    lapply(
      pairs,
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
          d$event[
            d$Group_pair == reference
          ] == 1
        )
        
        events_comparison <- sum(
          d$event[
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
          CI_lower_95 = NA_real_,
          CI_upper_95 = NA_real_,
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
            Surv(time, event) ~ Group_pair,
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
          HR = unname(exp(coef(fit)[1])),
          CI_lower_95 = ci[1, 1],
          CI_upper_95 = ci[1, 2],
          P_value = s$coefficients[
            1,
            "Pr(>|z|)"
          ]
        )
      }
    )
  )
}

pairwise <- pairwise_cox(
  df
)

pairwise_export <- pairwise |>
  mutate(
    `HR (95% CI)` = mapply(
      fmt_hr,
      HR,
      CI_lower_95,
      CI_upper_95
    ),
    
    `p-value` = vapply(
      P_value,
      fmt_p,
      character(1)
    ),
    
    `Comparison` = paste0(
      Comparison,
      " vs ",
      Reference
    ),
    
    `Comparator n/events` = paste0(
      n_comparison,
      " / ",
      events_comparison
    ),
    
    `Reference n/events` = paste0(
      n_reference,
      " / ",
      events_reference
    )
  ) |>
  select(
    Comparison,
    `HR (95% CI)`,
    `p-value`,
    `Comparator n/events`,
    `Reference n/events`
  )

pairwise_file <- file.path(
  FIGURE_DIR,
  "Fig5B_pairwise_HR_pvalues.csv"
)

readr::write_excel_csv(
  pairwise_export,
  pairwise_file
)

################################################################################
# 9. HR LABEL VS REFERENCE
################################################################################

ordered_levels <- c(
  HR_REFERENCE_GROUP,
  setdiff(
    levels(df$Group),
    HR_REFERENCE_GROUP
  )
)

cox_reference_data <- df |>
  mutate(
    Group_hr = factor(
      as.character(Group),
      levels = ordered_levels
    )
  )

cox_reference <- coxph(
  Surv(time, event) ~ Group_hr,
  data = cox_reference_data
)

cox_summary <- summary(
  cox_reference
)

ci <- exp(
  confint(
    cox_reference
  )
)

comparison_groups <- sub(
  "^Group_hr",
  "",
  rownames(
    cox_summary$coefficients
  )
)

hr_lines <- paste0(
  comparison_groups,
  " vs ",
  HR_REFERENCE_GROUP,
  ": HR = ",
  sprintf(
    "%.2f",
    exp(
      coef(
        cox_reference
      )
    )
  ),
  " (",
  sprintf(
    "%.2f",
    ci[, 1]
  ),
  "–",
  sprintf(
    "%.2f",
    ci[, 2]
  ),
  "), p = ",
  vapply(
    cox_summary$coefficients[
      ,
      "Pr(>|z|)"
    ],
    fmt_p,
    character(1)
  )
)

hr_label <- paste(
  c(
    paste0(
      "Reference: ",
      HR_REFERENCE_GROUP
    ),
    hr_lines
  ),
  collapse = "\n"
)

################################################################################
# 10. FIGURE 5B — KAPLAN–MEIER PLOT
################################################################################

group_levels_present <- levels(
  df$Group
)

palette_used <- unname(
  GROUP_COLORS[
    group_levels_present
  ]
)

if (anyNA(palette_used)) {
  stop(
    "A Figure 5B group is missing a defined color."
  )
}

x_max <- ceiling(
  max(
    df$time,
    na.rm = TRUE
  ) /
    BREAK_TIME_BY_YEARS
) * BREAK_TIME_BY_YEARS

km_plot <- survminer::ggsurvplot(
  fit = surv_fit,
  data = df,
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
  legend.labs = group_levels_present,
  
  break.time.by = BREAK_TIME_BY_YEARS,
  
  xlim = c(
    0,
    x_max
  ),
  
  ylim = c(
    0,
    1
  ),
  
  title =
    "Time to Progression by FISH 1q CN × unIGH (PI/IMiD)",
  
  subtitle =
    paste0(
      "Standard Risk (GEP70); TT3a, TT3b, TT4-S (n = ",
      nrow(df),
      ")"
    ),
  
  xlab =
    "Time from First Transplant (years)",
  
  ylab =
    "Cumulative Progression Probability",
  
  ggtheme =
    ggplot2::theme_bw(
      base_size = 16
    ) +
    ggplot2::theme(
      plot.title =
        ggplot2::element_text(
          hjust = 0.5,
          face = "bold"
        ),
      
      plot.subtitle =
        ggplot2::element_text(
          hjust = 0.5,
          size = 18
        )
    ),
  
  tables.theme =
    survminer::theme_cleantable()
)

km_plot$plot <- km_plot$plot +
  ggplot2::guides(
    color = ggplot2::guide_legend(
      nrow = LEGEND_NROW,
      byrow = LEGEND_BYROW
    )
  ) +
  ggplot2::theme(
    legend.text = ggplot2::element_text(
      size = LEGEND_TEXT_SIZE
    )
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 1),
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
  )

################################################################################
# 11. HR LABEL
################################################################################

if (SHOW_HR_LABEL) {
  
  plot_build <- ggplot2::ggplot_build(
    km_plot$plot
  )
  
  x_range <-
    plot_build$layout$panel_params[[1]]$x.range
  
  y_range <-
    plot_build$layout$panel_params[[1]]$y.range
  
  km_plot$plot <- km_plot$plot +
    ggplot2::annotate(
      geom = "text",
      
      x = x_range[2],
      
      y =
        y_range[2] -
        0.04 *
        (
          y_range[2] -
            y_range[1]
        ),
      
      label = hr_label,
      
      hjust = 1.08,
      vjust = 1,
      
      size = HR_LABEL_TEXT_SIZE,
      lineheight = 0.95
    )
}

################################################################################
# 12. RISK TABLE
################################################################################

km_plot$table <- km_plot$table +
  ggplot2::guides(
    color = "none",
    fill = "none"
  ) +
  ggplot2::theme(
    legend.position = "none",
    
    axis.title.y =
      ggplot2::element_blank(),
    
    plot.margin =
      ggplot2::margin(
        t = 0,
        r = 5.5,
        b = 5.5,
        l = 5.5,
        unit = "pt"
      )
  )

################################################################################
# 13. COMBINE
################################################################################

combined_plot <- ggpubr::ggarrange(
  km_plot$plot,
  km_plot$table,
  ncol = 1,
  heights = c(
    2,
    0.6
  ),
  align = "v"
)

print(
  combined_plot
)

################################################################################
# 14. SAVE FIGURE 5B
################################################################################

png_file <- file.path(
  FIGURE_DIR,
  "Fig5B_FISH1q_x_unIGH_ARD_TTP_KM.png"
)

pdf_file <- file.path(
  FIGURE_DIR,
  "Fig5B_FISH1q_x_unIGH_ARD_TTP_KM.pdf"
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
# 15. SAVE FIGURE 5B ANALYSIS COHORT
#
# Exact FISH1q × unIGH cohort used for Fig5B.
# Step 30 can read this file directly.
################################################################################

analysis_cohort_csv <- file.path(
  FIGURE_DIR,
  "Fig5B_FISH1q_x_unIGH_ARD_TTP_analysis-cohort.csv"
)

analysis_cohort_export <- df |>
  mutate(
    Group = as.character(
      Group
    ),
    Time_Years = time,
    Event = event
  )

readr::write_csv(
  analysis_cohort_export,
  analysis_cohort_csv
)

################################################################################
# 16. FINAL
################################################################################

cat(
  "\nStep 29 — Figure 5B complete.\n",
  
  "Endpoint: TTP\n",
  
  "Final cohort: ",
  nrow(df),
  
  "\nEvents: ",
  sum(df$event == 1),
  
  "\nOverall log-rank p: ",
  fmt_p(overall_p),
  
  "\nFISH1q column: ",
  FISH1Q_COL,
  
  "\nunIGH column: ",
  UNIGH_COL,
  
  "\nGEP70 column: ",
  GEP70_COL,
  
  "\nTTP time column: ",
  TIME_COL,
  
  "\nTTP event column: ",
  EVENT_COL,
  
  "\nAnalysis cohort CSV: ",
  analysis_cohort_csv,
  
  "\nPairwise HR/P CSV: ",
  pairwise_file,
  
  "\nPNG: ",
  png_file,
  
  "\nPDF: ",
  pdf_file,
  
  "\n"
)