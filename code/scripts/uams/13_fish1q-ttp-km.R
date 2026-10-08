################################################################################
# STEP 13 — FIGURE S2A
# FISH chromosome 1q copy number and Time to Progression
#
# Cohort: baseline CD138-selected NDMM, TT2–TT4, GEP70 standard risk,
#         available FISH 1q copy number.
#
# Required first:
#   Run through code/run_pipeline.R so DATA_ROOT and RESULTS_ROOT are defined.
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
# 1. SETTINGS
################################################################################

FIGURE_ID <- "FISH1q"

ESET_PATH <- file.path(
  DATA_ROOT,
  "processed",
  "ESET_uams_all_gep1q_fish_polyig_pcme.rds"
)

FIGURE_DIR <- file.path(
  RESULTS_ROOT,
  "figures",
  "Extended",
  FIGURE_ID
)

SOURCE_DATA_DIR <- file.path(
  FIGURE_DIR,
  "source-data"
)

TT_KEEP <- c(
  "TT2_NoThal",
  "TT2_Thal",
  "TT3a",
  "TT3b",
  "TT4_S-TT3"
)

GEP70_CUTOFF <- 0.66

BREAK_TIME_BY_YEARS <- 2

PLOT_WIDTH_IN <- 8
PLOT_HEIGHT_IN <- 6
PNG_DPI <- 300

COLOR_1Q2 <- "steelblue3"
COLOR_1Q_GAIN <- "red3"

if (!file.exists(ESET_PATH)) {
  
  stop(
    "Processed UAMS ExpressionSet not found: ",
    ESET_PATH
  )
}

dir.create(
  SOURCE_DATA_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


################################################################################
# 2. LOAD DATA
################################################################################

eset <- readRDS(
  ESET_PATH
)

df <- Biobase::pData(
  eset
)


################################################################################
# 3. RESTRICT COHORT
################################################################################

df <- df |>
  filter(
    sample_group == "RNAS_CD138_NDMM",
    tt_on_for_tx1 %in% TT_KEEP
  )


################################################################################
# 4. PREPARE VARIABLES
################################################################################

df <- df |>
  mutate(
    CHIPID = as.character(
      chipid
    ),
    
    GEP70 = as.numeric(
      GEP70
    ),
    
    FISH1q = as.numeric(
      FISH_1q_num
    ),
    
    time = as.numeric(
      YearsTTP
    ),
    
    event = as.numeric(
      CensTTP
    )
  ) |>
  filter(
    !is.na(
      CHIPID
    ),
    
    is.finite(
      GEP70
    ),
    
    GEP70 <= GEP70_CUTOFF,
    
    is.finite(
      FISH1q
    ),
    
    FISH1q == 2 |
      FISH1q >= 3,
    
    is.finite(
      time
    ),
    
    time >= 0,
    
    event %in% c(
      0,
      1
    )
  ) |>
  mutate(
    FISH1q_group = ifelse(
      FISH1q == 2,
      "1q2",
      "1q≥3"
    ),
    
    FISH1q_group = factor(
      FISH1q_group,
      levels = c(
        "1q2",
        "1q≥3"
      )
    )
  )


################################################################################
# 5. CHECK FINAL COHORT
################################################################################

if (anyDuplicated(
  df$CHIPID
)) {
  
  stop(
    "Duplicate CHIPIDs found in final cohort."
  )
}

if (any(
  df$GEP70 > GEP70_CUTOFF,
  na.rm = TRUE
)) {
  
  stop(
    "A GEP70-high sample entered the final cohort."
  )
}


################################################################################
# 6. SURVIVAL MODELS
################################################################################

survival_object <- Surv(
  time = df$time,
  event = df$event
)

survival_fit <- survfit(
  survival_object ~ FISH1q_group,
  data = df
)

cox_fit <- coxph(
  survival_object ~ FISH1q_group,
  data = df
)

cox_summary <- summary(
  cox_fit
)

hazard_ratio <- unname(
  exp(
    coef(
      cox_fit
    )
  )[1]
)

confidence_interval <- exp(
  confint(
    cox_fit
  )
)[1, ]

cox_p_value <- unname(
  cox_summary$coef[
    1,
    "Pr(>|z|)"
  ]
)

hr_label <- paste0(
  "FISH 1q≥3 vs 1q2: HR = ",
  sprintf(
    "%.2f",
    hazard_ratio
  ),
  " (",
  sprintf(
    "%.2f",
    confidence_interval[1]
  ),
  "–",
  sprintf(
    "%.2f",
    confidence_interval[2]
  ),
  ")",
  "\np = ",
  format.pval(
    cox_p_value,
    digits = 3
  )
)


################################################################################
# 7. SOURCE DATA
################################################################################

analysis_cohort_export <- df |>
  transmute(
    CHIPID,
    GEP70,
    TT_Protocol = tt_on_for_tx1,
    FISH1qCopies = FISH1q,
    FISH1q_Group = as.character(
      FISH1q_group
    ),
    Time_Years = time,
    Event = event,
    Endpoint = "TTP"
  )

cox_statistics_export <- tibble(
  Figure = FIGURE_ID,
  Endpoint = "TTP",
  Comparison = "FISH 1q≥3 vs 1q2",
  n = nrow(
    df
  ),
  events = sum(
    df$event == 1
  ),
  HR = hazard_ratio,
  CI_lower_95 = confidence_interval[1],
  CI_upper_95 = confidence_interval[2],
  p_value = cox_p_value
)

group_counts_export <- data.frame(
  FISH1q_group = names(
    table(
      df$FISH1q_group
    )
  ),
  
  n = as.integer(
    table(
      df$FISH1q_group
    )
  )
)

write_csv(
  analysis_cohort_export,
  file.path(
    SOURCE_DATA_DIR,
    paste0(
      FIGURE_ID,
      "_FISH1q_TTP_analysis-cohort.csv"
    )
  )
)

write_csv(
  cox_statistics_export,
  file.path(
    SOURCE_DATA_DIR,
    paste0(
      FIGURE_ID,
      "_FISH1q_TTP_cox-statistics.csv"
    )
  )
)

write_csv(
  group_counts_export,
  file.path(
    SOURCE_DATA_DIR,
    paste0(
      FIGURE_ID,
      "_FISH1q_TTP_group-counts.csv"
    )
  )
)


################################################################################
# 8. PLOT
################################################################################

max_time <- max(
  df$time,
  na.rm = TRUE
)

x_max <- ceiling(
  max_time / BREAK_TIME_BY_YEARS
) * BREAK_TIME_BY_YEARS

plot_title <- paste0(
  "Time to Progression",
  " by Chromosome 1q (FISH) Copy\n",
  "Standard Risk (GEP70); TT2–TT4",
  " (n = ",
  nrow(
    df
  ),
  ")"
)

km_plot <- survminer::ggsurvplot(
  fit = survival_fit,
  data = df,
  
  fun = "event",
  
  pval = TRUE,
  
  pval.coord = c(
    BREAK_TIME_BY_YEARS / 2,
    0.85
  ),
  
  pval.size = 5,
  
  risk.table = TRUE,
  
  risk.table.col = "strata",
  
  risk.table.y.text = TRUE,
  
  risk.table.y.text.col = TRUE,
  
  risk.table.height = 0.18,
  
  risk.table.title = "No. at risk",
  
  cumevents = FALSE,
  
  palette = c(
    COLOR_1Q2,
    COLOR_1Q_GAIN
  ),
  
  legend.title = "FISH 1q copy class",
  
  legend.labs = c(
    "1q2",
    "1q≥3"
  ),
  
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
  
  xlab = "Time from First Transplant (years)",
  
  ylab = "Cumulative Progression Probability",
  
  ggtheme = ggplot2::theme_bw(
    base_size = 16
  ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = 18
      )
    ),
  
  tables.theme = survminer::theme_cleantable(
    base_size = 13
  )
)

km_plot$plot <- km_plot$plot +
  ggplot2::annotate(
    geom = "text",
    x = x_max * 0.98,
    y = 0.96,
    label = hr_label,
    hjust = 1,
    vjust = 1,
    size = 4.5
  )

km_plot$table <- km_plot$table +
  ggplot2::guides(
    color = "none",
    fill = "none"
  ) +
  ggplot2::theme(
    legend.position = "none",
    axis.title.y = ggplot2::element_blank(),
    
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
    0.4
  ),
  
  align = "v"
)

print(
  combined_plot
)


################################################################################
# 9. SAVE
################################################################################

png_file <- file.path(
  FIGURE_DIR,
  paste0(
    FIGURE_ID,
    "_FISH1q_TTP_KM.png"
  )
)

pdf_file <- file.path(
  FIGURE_DIR,
  paste0(
    FIGURE_ID,
    "_FISH1q_TTP_KM.pdf"
  )
)

ggsave(
  png_file,
  combined_plot,
  width = PLOT_WIDTH_IN,
  height = PLOT_HEIGHT_IN,
  units = "in",
  dpi = PNG_DPI,
  bg = "white"
)

ggsave(
  pdf_file,
  combined_plot,
  width = PLOT_WIDTH_IN,
  height = PLOT_HEIGHT_IN,
  units = "in",
  device = grDevices::cairo_pdf,
  bg = "white"
)


################################################################################
# 10. SUMMARY
################################################################################

message(
  "\nStep 13 complete — Figure S2A",
  "\nEndpoint: TTP",
  "\nFinal cohort: ",
  nrow(
    df
  ),
  "\nEvents: ",
  sum(
    df$event == 1
  ),
  "\nFISH 1q2: ",
  sum(
    df$FISH1q_group == "1q2"
  ),
  "\nFISH 1q≥3: ",
  sum(
    df$FISH1q_group == "1q≥3"
  ),
  "\nHR: ",
  sprintf(
    "%.3f",
    hazard_ratio
  ),
  "\n95% CI: ",
  sprintf(
    "%.3f",
    confidence_interval[1]
  ),
  "–",
  sprintf(
    "%.3f",
    confidence_interval[2]
  ),
  "\nCox p-value: ",
  format(
    cox_p_value,
    scientific = TRUE
  ),
  "\nPNG: ",
  png_file,
  "\nPDF: ",
  pdf_file,
  "\n"
)