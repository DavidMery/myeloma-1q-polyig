################################################################################
# STEP 28 — FIGURE 5A
# TTP BY GEP1qcopy_lasso × PolyIG
#
# Cohort:
#   1. Baseline RNAS_CD138_NDMM
#   2. TT3a, TT3b, TT4_S-TT3
#   3. GEP70 <= 0.66
#   4. Valid FISH_1q_num and IG_ARD
#   5. Valid GEP1qcopy_lasso and PolyIG_Score
#   6. Valid YearsTTP and CensTTP
#
# The analysis-cohort CSV retains all original pData columns plus derived
# analysis columns.
################################################################################

library(Biobase)
library(survival)
library(survminer)
library(dplyr)
library(tibble)
library(readr)
library(ggplot2)
library(ggpubr)


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

FIGURE_ID <- "Fig5A"
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

# Exact pData columns
SAMPLE_ID_COL <- "chipid"
SAMPLE_GROUP_COL <- "sample_group"
TT_COL <- "TT_On_for_Tx1"

GEP1Q_COL <- "GEP1qcopy_lasso"
POLYIG_COL <- "PolyIG_Score"
GEP70_COL <- "GEP70"

FISH1Q_COL <- "FISH_1q_num"
UNIGH_COL <- "IG_ARD"

TIME_COL <- "YearsTTP"
EVENT_COL <- "CensTTP"

GEP70_CUTOFF <- 0.66
POLYIG_CUT <- 11

TT_KEEP <- c(
  "TT3a",
  "TT3b",
  "TT4_S-TT3"
)

################################################################################
# 2. PLOT SETTINGS
################################################################################

GROUP_CN2_LOW <- "CN2/PolyIG Low"
GROUP_CN2_HIGH <- "CN2/PolyIG High"
GROUP_CN3_ALL <- "CN≥3/PolyIG ALL"

GROUP_LEVELS <- c(
  GROUP_CN2_HIGH,
  GROUP_CN2_LOW,
  GROUP_CN3_ALL
)

GROUP_COLORS <- c(
  "CN2/PolyIG High" = "steelblue2",
  "CN2/PolyIG Low" = "grey50",
  "CN≥3/PolyIG ALL" = "red3"
)

BREAK_TIME_BY_YEARS <- 2

HR_REFERENCE_GROUP <- GROUP_CN2_HIGH

SHOW_LOGRANK_P_VALUE <- TRUE
LOGRANK_P_VALUE_COORD <- c(0.5, 0.85)
LOGRANK_P_VALUE_SIZE <- 7

SHOW_HR_LABEL <- TRUE
HR_LABEL_TEXT_SIZE <- 3.6

LEGEND_NROW <- 1L
LEGEND_BYROW <- TRUE
LEGEND_TEXT_SIZE <- 13

PLOT_WIDTH_IN <- 9
PLOT_HEIGHT_IN <- 8
PNG_DPI <- 300

################################################################################
# 3. HELPERS
################################################################################

fmt_p <- function(p) {
  if (!is.finite(p)) return("NA")
  
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
# 4. ALL PAIRWISE COX COMPARISONS
################################################################################

pairwise_cox <- function(d) {
  
  lev <- levels(
    droplevels(
      d$Group
    )
  )
  
  if (length(lev) < 2) {
    return(tibble())
  }
  
  pairs <- combn(
    lev,
    2,
    simplify = FALSE
  )
  
  bind_rows(
    lapply(
      pairs,
      function(pair) {
        
        ref <- pair[1]
        cmp <- pair[2]
        
        x <- d |>
          filter(
            Group %in% c(ref, cmp)
          ) |>
          mutate(
            Group_pair = factor(
              as.character(Group),
              levels = c(ref, cmp)
            )
          )
        
        n_ref <- sum(
          x$Group_pair == ref
        )
        
        n_cmp <- sum(
          x$Group_pair == cmp
        )
        
        e_ref <- sum(
          x$event[
            x$Group_pair == ref
          ] == 1
        )
        
        e_cmp <- sum(
          x$event[
            x$Group_pair == cmp
          ] == 1
        )
        
        if (
          n_ref < 2 ||
          n_cmp < 2 ||
          e_ref + e_cmp == 0
        ) {
          return(
            tibble(
              reference = ref,
              comparison = cmp,
              n_reference = n_ref,
              n_comparison = n_cmp,
              events_reference = e_ref,
              events_comparison = e_cmp,
              HR = NA_real_,
              CI_lower_95 = NA_real_,
              CI_upper_95 = NA_real_,
              p_value = NA_real_,
              status = "insufficient_n_or_events"
            )
          )
        }
        
        fit <- tryCatch(
          coxph(
            Surv(time, event) ~ Group_pair,
            data = x
          ),
          error = function(e) NULL
        )
        
        if (is.null(fit)) {
          return(
            tibble(
              reference = ref,
              comparison = cmp,
              n_reference = n_ref,
              n_comparison = n_cmp,
              events_reference = e_ref,
              events_comparison = e_cmp,
              HR = NA_real_,
              CI_lower_95 = NA_real_,
              CI_upper_95 = NA_real_,
              p_value = NA_real_,
              status = "coxph_error"
            )
          )
        }
        
        s <- summary(fit)
        ci <- exp(confint(fit))
        
        tibble(
          reference = ref,
          comparison = cmp,
          n_reference = n_ref,
          n_comparison = n_cmp,
          events_reference = e_ref,
          events_comparison = e_cmp,
          HR = exp(coef(fit))[1],
          CI_lower_95 = ci[1, 1],
          CI_upper_95 = ci[1, 2],
          p_value = s$coefficients[
            1,
            "Pr(>|z|)"
          ],
          status = "ok"
        )
      }
    )
  )
}

################################################################################
# 5. LOAD EXPRESSIONSET
################################################################################

if (!file.exists(ESET_PATH)) {
  stop(
    "Cannot find ExpressionSet:\n",
    ESET_PATH
  )
}

eset <- readRDS(
  ESET_PATH
)

if (!methods::is(eset, "ExpressionSet")) {
  stop(
    "ESET_PATH did not contain an ExpressionSet."
  )
}

# Keep complete pData so the final cohort CSV retains all original columns.
pdat <- as.data.frame(
  Biobase::pData(eset),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

################################################################################
# 6. CHECK EXACT REQUIRED COLUMNS
################################################################################

required_columns <- c(
  SAMPLE_ID_COL,
  SAMPLE_GROUP_COL,
  TT_COL,
  GEP1Q_COL,
  POLYIG_COL,
  GEP70_COL,
  FISH1Q_COL,
  UNIGH_COL,
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
# 7. BUILD MATCHED FISH1q × unIGH COHORT
################################################################################

df <- pdat |>
  tibble::as_tibble() |>
  mutate(
    CHIPID = .data[[SAMPLE_ID_COL]],
    TT_protocol = .data[[TT_COL]],
    
    FISH1q = as.numeric(
      .data[[FISH1Q_COL]]
    ),
    
    unIGH_ARD = as.numeric(
      .data[[UNIGH_COL]]
    ),
    
    GEP1qcopy_lasso = .data[[GEP1Q_COL]],
    
    PolyIG = as.numeric(
      .data[[POLYIG_COL]]
    ),
    
    GEP70 = as.numeric(
      .data[[GEP70_COL]]
    ),
    
    time = as.numeric(
      .data[[TIME_COL]]
    ),
    
    event = as.numeric(
      .data[[EVENT_COL]]
    )
  ) |>
  filter(
    .data[[SAMPLE_GROUP_COL]] == "RNAS_CD138_NDMM",
    
    TT_protocol %in% TT_KEEP,
    
    # Matched FISH1q × unIGH clinical cohort.
    is.finite(FISH1q),
    is.finite(unIGH_ARD),
    
    # GEP70 standard risk.
    is.finite(GEP70),
    GEP70 <= GEP70_CUTOFF,
    
    # TTP.
    is.finite(time),
    time >= 0,
    event %in% c(0, 1),
    
    # GEP1q / PolyIG.
    GEP1qcopy_lasso %in% c(
      "1q2",
      "1q3plus"
    ),
    is.finite(PolyIG),
    
    !is.na(CHIPID)
  ) |>
  mutate(
    PolyIG_group = if_else(
      PolyIG > POLYIG_CUT,
      "High",
      "Low"
    )
  )

if (!nrow(df)) {
  stop(
    "No samples remain after the matched FISH1q × unIGH cohort filter."
  )
}

if (anyDuplicated(df$CHIPID)) {
  stop(
    "Duplicate CHIPIDs are present in the Fig5A cohort."
  )
}

################################################################################
# 8. FIGURE 5A GROUPS
################################################################################

df <- df |>
  mutate(
    Group = case_when(
      
      GEP1qcopy_lasso == "1q2" &
        PolyIG_group == "High" ~
        GROUP_CN2_HIGH,
      
      GEP1qcopy_lasso == "1q2" &
        PolyIG_group == "Low" ~
        GROUP_CN2_LOW,
      
      GEP1qcopy_lasso == "1q3plus" ~
        GROUP_CN3_ALL,
      
      TRUE ~
        NA_character_
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

if (nlevels(df$Group) < 2) {
  stop(
    "Fewer than two Fig5A groups remain."
  )
}

if (!HR_REFERENCE_GROUP %in% levels(df$Group)) {
  stop(
    "HR reference group is absent: ",
    HR_REFERENCE_GROUP
  )
}

################################################################################
# 9. SURVIVAL
################################################################################

fit <- survfit(
  Surv(time, event) ~ Group,
  data = df
)

groups <- levels(
  df$Group
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
# 10. ALL PAIRWISE COX COMPARISONS
################################################################################

pairwise <- pairwise_cox(
  df
)

pairwise_csv <- file.path(
  FIGURE_DIR,
  "Fig5A_GEP1q_lasso_x_PolyIG_TTP_pairwise-HR.csv"
)

readr::write_excel_csv(
  pairwise,
  pairwise_csv
)

################################################################################
# 11. COX HR LABEL
################################################################################

hr_fit <- coxph(
  Surv(time, event) ~
    relevel(
      Group,
      ref = HR_REFERENCE_GROUP
    ),
  data = df
)

hr_summary <- summary(
  hr_fit
)

hr_ci <- exp(
  confint(
    hr_fit
  )
)

hr_names <- rownames(
  hr_summary$coefficients
)

hr_names <- sub(
  "^relevel\\(Group, ref = HR_REFERENCE_GROUP\\)",
  "",
  hr_names
)

hr_label_lines <- paste0(
  hr_names,
  " vs ",
  HR_REFERENCE_GROUP,
  ": HR = ",
  sprintf(
    "%.2f",
    exp(coef(hr_fit))
  ),
  " (",
  sprintf(
    "%.2f",
    hr_ci[, 1]
  ),
  "–",
  sprintf(
    "%.2f",
    hr_ci[, 2]
  ),
  "), p = ",
  vapply(
    hr_summary$coefficients[
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
    hr_label_lines
  ),
  collapse = "\n"
)

################################################################################
# 12. KM PLOT
################################################################################

palette_used <- unname(
  GROUP_COLORS[
    groups
  ]
)

if (anyNA(palette_used)) {
  stop(
    "A Fig5A group is missing a defined color."
  )
}

x_max <- ceiling(
  max(
    df$time,
    na.rm = TRUE
  ) /
    BREAK_TIME_BY_YEARS
) *
  BREAK_TIME_BY_YEARS

km <- survminer::ggsurvplot(
  fit = fit,
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
  legend.labs = groups,
  
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
    "Time to Progression by GEP 1q CN × PolyIG (PI/IMiD)",
  
  xlab =
    "Time from First Transplant (years)",
  
  ylab =
    "Cumulative Progression Probability",
  
  ggtheme =
    theme_bw(
      base_size = 16
    ) +
    theme(
      plot.title =
        element_text(
          hjust = 0.5,
          face = "bold"
        ),
      
      plot.subtitle =
        element_text(
          hjust = 0.5,
          size = 18
        )
    ),
  
  tables.theme =
    survminer::theme_cleantable()
)

km$plot <- km$plot +
  labs(
    subtitle =
      paste0(
        "Standard Risk (GEP70); TT3a, TT3b, TT4-S",
        " (n = ",
        nrow(df),
        ")"
      )
  ) +
  guides(
    color =
      guide_legend(
        nrow = LEGEND_NROW,
        byrow = LEGEND_BYROW
      )
  ) +
  theme(
    legend.text =
      element_text(
        size = LEGEND_TEXT_SIZE
      )
  ) +
  scale_y_continuous(
    limits = c(
      0,
      1
    ),
    breaks =
      seq(
        0,
        1,
        0.2
      ),
    labels =
      function(x) {
        sprintf(
          "%.1f",
          x
        )
      }
  )

################################################################################
# 13. HR LABEL ON PLOT
################################################################################

if (SHOW_HR_LABEL) {
  
  b <- ggplot2::ggplot_build(
    km$plot
  )
  
  xr <- b$layout$panel_params[[1]]$x.range
  yr <- b$layout$panel_params[[1]]$y.range
  
  km$plot <- km$plot +
    annotate(
      "text",
      x = xr[2],
      y =
        yr[2] -
        0.04 *
        diff(yr),
      label = hr_label,
      hjust = 1.08,
      vjust = 1,
      size = HR_LABEL_TEXT_SIZE,
      lineheight = 0.95
    )
}

################################################################################
# 14. RISK TABLE
################################################################################

km$table <- km$table +
  guides(
    color = "none",
    fill = "none"
  ) +
  theme(
    legend.position = "none",
    
    axis.title.y =
      element_blank(),
    
    plot.margin =
      margin(
        0,
        5.5,
        5.5,
        5.5,
        unit = "pt"
      )
  )

################################################################################
# 15. COMBINE
################################################################################

combined <- ggpubr::ggarrange(
  km$plot,
  km$table,
  ncol = 1,
  heights = c(
    2,
    0.6
  ),
  align = "v"
)

print(
  combined
)

################################################################################
# 16. SAVE FIGURE
################################################################################

png_file <- file.path(
  FIGURE_DIR,
  "Fig5A_GEP1q_lasso_x_PolyIG_TTP_KM.png"
)

pdf_file <- file.path(
  FIGURE_DIR,
  "Fig5A_GEP1q_lasso_x_PolyIG_TTP_KM.pdf"
)

ggplot2::ggsave(
  png_file,
  combined,
  width = PLOT_WIDTH_IN,
  height = PLOT_HEIGHT_IN,
  units = "in",
  dpi = PNG_DPI,
  bg = "white"
)

ggplot2::ggsave(
  pdf_file,
  combined,
  width = PLOT_WIDTH_IN,
  height = PLOT_HEIGHT_IN,
  units = "in",
  device = "pdf",
  bg = "white"
)

################################################################################
# 17. SAVE COMPLETE ANALYSIS COHORT
#
# All original pData columns remain.
# Derived analysis columns are appended.
################################################################################

analysis_cohort_csv <- file.path(
  FIGURE_DIR,
  "Fig5A_GEP1q_lasso_x_PolyIG_TTP_analysis-cohort.csv"
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
# 18. FINAL SUMMARY
################################################################################

message(
  "\n############################################################",
  "\nStep 28 / Fig5A complete",
  "\n############################################################",
  
  "\nEndpoint: TTP",
  
  "\nAnalysis: GEP1q LASSO × PolyIG",
  
  "\nMatched FISH1q × unIGH cohort: yes",
  
  "\nFISH1q column: ",
  FISH1Q_COL,
  
  "\nunIGH ARD column: ",
  UNIGH_COL,
  
  "\nGEP1q column: ",
  GEP1Q_COL,
  
  "\nPolyIG column: ",
  POLYIG_COL,
  
  "\nGEP70 column: ",
  GEP70_COL,
  
  "\nTTP time column: ",
  TIME_COL,
  
  "\nTTP event column: ",
  EVENT_COL,
  
  "\nFinal n: ",
  nrow(df),
  
  "\nEvents: ",
  sum(df$event == 1),
  
  "\nOverall log-rank p: ",
  fmt_p(overall_p),
  
  "\nPairwise HR/P CSV: ",
  pairwise_csv,
  
  "\nAnalysis cohort CSV: ",
  analysis_cohort_csv,
  
  "\nPNG: ",
  png_file,
  
  "\nPDF: ",
  pdf_file,
  
  "\n"
)