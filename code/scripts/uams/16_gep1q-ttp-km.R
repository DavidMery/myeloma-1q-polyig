################################################################################
# STEP 16B — FIGURE S3A
# Kaplan-Meier analysis of TTP by GEP1qcopy_lasso class
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

FIGURE_ID <- "GEP1q"

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

COLOR_CN2 <- "steelblue3"
COLOR_CN3PLUS <- "red3"


################################################################################
# 2. VERIFY REQUIRED INPUT
################################################################################

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
# 3. LOAD DATA
################################################################################

eset <- readRDS(
  ESET_PATH
)

df <- Biobase::pData(
  eset
)


################################################################################
# 4. RESTRICT COHORT
################################################################################

df <- df |>
  filter(
    sample_group == "RNAS_CD138_NDMM",
    tt_on_for_tx1 %in% TT_KEEP
  )


################################################################################
# 5. PREPARE VARIABLES
################################################################################

df <- df |>
  mutate(
    CHIPID = as.character(
      chipid
    ),
    
    GEP70 = as.numeric(
      GEP70
    ),
    
    GEP1q_raw = trimws(
      as.character(
        GEP1qcopy_lasso
      )
    ),
    
    time = as.numeric(
      YearsTTP
    ),
    
    event = as.numeric(
      CensTTP
    )
  )


################################################################################
# 6. RESTRICT TO GEP70 STANDARD RISK
################################################################################

df <- df |>
  filter(
    is.finite(
      GEP70
    ),
    
    GEP70 <= GEP70_CUTOFF
  )


################################################################################
# 7. DEFINE GEP1q CLASS
################################################################################

df <- df |>
  mutate(
    GEP1q_class = case_when(
      
      GEP1q_raw %in% c(
        "1q2",
        "CN2",
        "2"
      ) ~ "CN2",
      
      GEP1q_raw %in% c(
        "1q3plus",
        "CN≥3",
        "CN>=3",
        "3+",
        "3",
        "4"
      ) ~ "CN≥3",
      
      TRUE ~ NA_character_
    ),
    
    GEP1q_group = factor(
      GEP1q_class,
      levels = c(
        "CN2",
        "CN≥3"
      )
    )
  )


################################################################################
# 8. FINAL ANALYSIS COHORT
################################################################################

df <- df |>
  filter(
    !is.na(
      CHIPID
    ),
    
    is.finite(
      time
    ),
    
    time >= 0,
    
    event %in% c(
      0,
      1
    ),
    
    !is.na(
      GEP1q_group
    )
  )


if (anyDuplicated(
  df$CHIPID
)) {
  
  stop(
    "Duplicate CHIPIDs found in final cohort."
  )
}


if (nlevels(
  droplevels(
    df$GEP1q_group
  )
) != 2) {
  
  stop(
    "Both CN2 and CN≥3 groups must contain samples."
  )
}


################################################################################
# 9. SURVIVAL ANALYSIS
################################################################################

survival_object <- Surv(
  time = df$time,
  event = df$event
)

survival_fit <- survfit(
  survival_object ~ GEP1q_group,
  data = df
)

cox_fit <- coxph(
  survival_object ~ GEP1q_group,
  data = df
)

cox_summary <- summary(
  cox_fit
)

hazard_ratio <- exp(
  coef(
    cox_fit
  )
)[1]

confidence_interval <- exp(
  confint(
    cox_fit
  )
)[1, ]

cox_p_value <- cox_summary$coef[
  1,
  "Pr(>|z|)"
]


################################################################################
# 10. HR LABEL
################################################################################

hr_label <- paste0(
  "GEP1q CN≥3 vs CN2: HR = ",
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
# 11. SOURCE DATA
################################################################################

analysis_cohort_export <- df |>
  transmute(
    CHIPID,
    
    GEP70,
    
    TT_Protocol = tt_on_for_tx1,
    
    GEP1qcopy_lasso_Raw = GEP1q_raw,
    
    GEP1qcopy_lasso_Group = as.character(
      GEP1q_group
    ),
    
    Time_Years = time,
    
    Event = event,
    
    Endpoint = "TTP"
  )


cox_statistics_export <- tibble(
  Figure = FIGURE_ID,
  
  Step = "16b",
  
  Endpoint = "TTP",
  
  Comparison = "GEP1qcopy_lasso CN≥3 vs CN2",
  
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


group_counts_export <- tibble(
  GEP1q_group = names(
    table(
      df$GEP1q_group
    )
  ),
  
  n = as.integer(
    table(
      df$GEP1q_group
    )
  )
)


write_csv(
  analysis_cohort_export,
  file.path(
    SOURCE_DATA_DIR,
    paste0(
      FIGURE_ID,
      "_GEP1q_lasso_TTP_analysis-cohort.csv"
    )
  )
)


write_csv(
  cox_statistics_export,
  file.path(
    SOURCE_DATA_DIR,
    paste0(
      FIGURE_ID,
      "_GEP1q_lasso_TTP_cox-statistics.csv"
    )
  )
)


write_csv(
  group_counts_export,
  file.path(
    SOURCE_DATA_DIR,
    paste0(
      FIGURE_ID,
      "_GEP1q_lasso_TTP_group-counts.csv"
    )
  )
)


################################################################################
# 12. PLOT
################################################################################

max_time <- max(
  df$time,
  na.rm = TRUE
)

x_max <- ceiling(
  max_time /
    BREAK_TIME_BY_YEARS
) *
  BREAK_TIME_BY_YEARS


plot_title <- paste0(
  "Time to Progression by Penalized GEP1q Copy Class\n",
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
  
  conf.int = FALSE,
  
  censor = TRUE,
  
  palette = c(
    COLOR_CN2,
    COLOR_CN3PLUS
  ),
  
  legend.title = "GEP1q",
  
  legend.labs = c(
    "CN2",
    "CN≥3"
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
  
  ggtheme = theme_bw(
    base_size = 16
  ) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 18
      )
    ),
  
  tables.theme = theme_cleantable(
    base_size = 13
  )
)


km_plot$plot <- km_plot$plot +
  annotate(
    geom = "text",
    
    x = x_max * 0.98,
    
    y = 0.96,
    
    label = hr_label,
    
    hjust = 1,
    
    vjust = 1,
    
    size = 4.5
  )


km_plot$table <- km_plot$table +
  guides(
    color = "none",
    colour = "none",
    fill = "none",
    shape = "none"
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


combined_plot <- ggarrange(
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
# 13. SAVE FIGURE
################################################################################

png_file <- file.path(
  FIGURE_DIR,
  paste0(
    FIGURE_ID,
    "_GEP1q_lasso_TTP_KM.png"
  )
)

pdf_file <- file.path(
  FIGURE_DIR,
  paste0(
    FIGURE_ID,
    "_GEP1q_lasso_TTP_KM.pdf"
  )
)


ggsave(
  filename = png_file,
  plot = combined_plot,
  width = PLOT_WIDTH_IN,
  height = PLOT_HEIGHT_IN,
  units = "in",
  dpi = PNG_DPI,
  bg = "white"
)


ggsave(
  filename = pdf_file,
  plot = combined_plot,
  width = PLOT_WIDTH_IN,
  height = PLOT_HEIGHT_IN,
  units = "in",
  device = "pdf",
  bg = "white"
)


################################################################################
# 14. SUMMARY
################################################################################

message(
  "\nStep 16b complete",
  
  "\nFinal cohort: ",
  nrow(
    df
  ),
  
  "\nEvents: ",
  sum(
    df$event == 1
  ),
  
  "\nCN2: ",
  sum(
    df$GEP1q_group == "CN2"
  ),
  
  "\nCN≥3: ",
  sum(
    df$GEP1q_group == "CN≥3"
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