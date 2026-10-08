################################################################################
# STEP 09 — FIGURE 1A
#
# Purpose:
#   1. Restrict to RNAS CD138+ NDMM treated on TT2-TT5.
#   2. Divide patients by GEP70 ≤ 0.66 versus > 0.66.
#   3. Fit Kaplan-Meier curves for TTP.
#   4. Calculate hazard ratio and Cox P value.
#   5. Plot cumulative event probability with numbers at risk.
#   6. Save the figure and essential source data.
#
# Required first:
#   Run through code/run_pipeline.R so DATA_ROOT and RESULTS_ROOT are defined.
################################################################################

library(Biobase)
library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(ggpubr)


################################################################################
# 1. VERIFY PIPELINE ROOTS
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
# 2. SETTINGS
################################################################################

GEP70_CUTOFF <- 0.66

TT_KEEP <- c(
  "TT2_NoThal",
  "TT2_Thal",
  "TT3a",
  "TT3b",
  "TT4_S-TT3",
  "TT5"
)


################################################################################
# 3. FILES
################################################################################

ESET_FILE <- file.path(
  DATA_ROOT,
  "processed",
  "ESET_uams_all_gep1q_fish_polyig_pcme.rds"
)

OUT_DIR <- file.path(
  RESULTS_ROOT,
  "figures",
  "Fig1",
  "Fig1A"
)

if (!file.exists(ESET_FILE)) {
  
  stop(
    "Processed UAMS ExpressionSet not found: ",
    ESET_FILE
  )
}

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


################################################################################
# 4. LOAD PHENOTYPE DATA
################################################################################

eset <- readRDS(
  ESET_FILE
)

df <- pData(
  eset
)


################################################################################
# 5. SELECT TTP COLUMNS
################################################################################

TIME_COL <- "YearsTTP"
EVENT_COL <- "CensTTP"


################################################################################
# 6. RESTRICT TO RNAS CD138+ NDMM AND TT2-TT5
################################################################################

df <- df |>
  filter(
    sample_group == "RNAS_CD138_NDMM",
    tt_on_for_tx1 %in% TT_KEEP
  )


################################################################################
# 7. PREPARE SURVIVAL DATA
################################################################################

df <- df |>
  mutate(
    GEP70 = as.numeric(GEP70),
    time = as.numeric(.data[[TIME_COL]]),
    event = as.numeric(.data[[EVENT_COL]])
  ) |>
  filter(
    is.finite(GEP70),
    is.finite(time),
    time >= 0,
    event %in% c(0, 1)
  ) |>
  mutate(
    GEP70_group = ifelse(
      GEP70 <= GEP70_CUTOFF,
      "GEP70 ≤ 0.66",
      "GEP70 > 0.66"
    ),
    GEP70_group = factor(
      GEP70_group,
      levels = c(
        "GEP70 ≤ 0.66",
        "GEP70 > 0.66"
      )
    )
  )


################################################################################
# 8. FIT KAPLAN-MEIER MODEL
################################################################################

surv_object <- Surv(
  time = df$time,
  event = df$event
)

fit <- survfit(
  surv_object ~ GEP70_group,
  data = df
)


################################################################################
# 9. FIT COX MODEL
################################################################################

cox <- coxph(
  surv_object ~ GEP70_group,
  data = df
)

cox_summary <- summary(
  cox
)

HR <- exp(
  coef(
    cox
  )
)[1]

CI <- exp(
  confint(
    cox
  )
)[1, ]

P_VALUE <- cox_summary$coef[
  1,
  "Pr(>|z|)"
]

HR_LABEL <- paste0(
  "GEP70 > 0.66 vs GEP70 ≤ 0.66: HR = ",
  sprintf("%.2f", HR),
  " (",
  sprintf("%.2f", CI[1]),
  "–",
  sprintf("%.2f", CI[2]),
  ")\nP = ",
  format.pval(
    P_VALUE,
    digits = 3
  )
)


################################################################################
# 10. CREATE CUMULATIVE-EVENT KAPLAN-MEIER PLOT
################################################################################

x_max <- ceiling(
  max(
    df$time,
    na.rm = TRUE
  )
) + 1

plot_title <- paste0(
  "Time to Progression",
  " by High and Standard Risk (GEP70)\n",
  "Total Therapy: TT2-TT5",
  " (n = ",
  nrow(df),
  ")"
)

p <- ggsurvplot(
  fit,
  data = df,
  
  fun = "event",
  
  pval = TRUE,
  pval.coord = c(
    1,
    0.85
  ),
  pval.size = 5,
  
  risk.table = TRUE,
  risk.table.col = "strata",
  risk.table.y.text.col = TRUE,
  risk.table.height = 0.18,
  risk.table.title = "No. at risk",
  
  cumevents = FALSE,
  
  palette = c(
    "steelblue3",
    "black"
  ),
  
  legend.title = "GEP70 group",
  legend.labs = c(
    "GEP70 ≤ 0.66",
    "GEP70 > 0.66"
  ),
  
  break.time.by = 2,
  
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
        face = "bold"
      )
    ),
  
  tables.theme = theme_cleantable()
)


################################################################################
# 11. ADD HAZARD-RATIO ANNOTATION
################################################################################

plot_build <- ggplot_build(
  p$plot
)

x_range <- plot_build$layout$panel_params[[1]]$x.range
y_range <- plot_build$layout$panel_params[[1]]$y.range

p$plot <- p$plot +
  annotate(
    "text",
    x = x_range[2],
    y = y_range[2] - 0.05 * diff(y_range),
    label = HR_LABEL,
    hjust = 1.05,
    vjust = 1,
    size = 5
  )


################################################################################
# 12. COMBINE PLOT AND RISK TABLE
################################################################################

combined <- ggarrange(
  p$plot,
  
  p$table +
    guides(
      color = "none",
      fill = "none"
    ) +
    theme(
      legend.position = "none",
      
      plot.margin = margin(
        t = 0,
        r = 5.5,
        b = 5.5,
        l = 5.5,
        unit = "pt"
      )
    ),
  
  ncol = 1,
  
  heights = c(
    2,
    0.4
  )
)


################################################################################
# 13. SAVE FIGURE
################################################################################

FIGURE_STUB <- "Fig1A_GEP70_TTP"

ggsave(
  file.path(
    OUT_DIR,
    paste0(
      FIGURE_STUB,
      ".png"
    )
  ),
  combined,
  width = 8.5,
  height = 6,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(
    OUT_DIR,
    paste0(
      FIGURE_STUB,
      ".pdf"
    )
  ),
  combined,
  width = 8.5,
  height = 6,
  device = grDevices::pdf,
  useDingbats = FALSE,
  bg = "white"
)


################################################################################
# 14. SAVE SOURCE DATA
################################################################################

source_data <- data.frame(
  sample_id = rownames(df),
  GEP70 = df$GEP70,
  GEP70_group = df$GEP70_group,
  time = df$time,
  event = df$event
)

write.csv(
  source_data,
  file.path(
    OUT_DIR,
    paste0(
      FIGURE_STUB,
      "_source-data.csv"
    )
  ),
  row.names = FALSE
)


################################################################################
# 15. SAVE STATISTICS
################################################################################

statistics <- data.frame(
  Endpoint = "TTP",
  N = nrow(df),
  
  GEP70_standard_n = sum(
    df$GEP70_group == "GEP70 ≤ 0.66"
  ),
  
  GEP70_high_n = sum(
    df$GEP70_group == "GEP70 > 0.66"
  ),
  
  Hazard_Ratio = HR,
  CI_Lower = CI[1],
  CI_Upper = CI[2],
  Cox_P = P_VALUE
)

write.csv(
  statistics,
  file.path(
    OUT_DIR,
    paste0(
      FIGURE_STUB,
      "_statistics.csv"
    )
  ),
  row.names = FALSE
)


################################################################################
# 16. COMPLETE
################################################################################

cat(
  "\n",
  "============================================================\n",
  "Step 09 complete — Figure 1A created successfully\n",
  "============================================================\n",
  sep = ""
)

cat(
  "Endpoint: TTP\n"
)

cat(
  "Samples analyzed: ",
  nrow(df),
  "\n",
  sep = ""
)

cat(
  "GEP70 ≤ 0.66: ",
  sum(
    df$GEP70_group == "GEP70 ≤ 0.66"
  ),
  "\n",
  sep = ""
)

cat(
  "GEP70 > 0.66: ",
  sum(
    df$GEP70_group == "GEP70 > 0.66"
  ),
  "\n",
  sep = ""
)

cat(
  "HR: ",
  round(
    HR,
    2
  ),
  "\n",
  sep = ""
)

cat(
  "95% CI: ",
  sprintf(
    "%.2f–%.2f",
    CI[1],
    CI[2]
  ),
  "\n",
  sep = ""
)

cat(
  "P: ",
  format.pval(
    P_VALUE,
    digits = 3
  ),
  "\n",
  sep = ""
)

cat(
  "Figure/output directory: ",
  OUT_DIR,
  "\n",
  sep = ""
)