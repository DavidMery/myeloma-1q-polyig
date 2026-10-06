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
################################################################################

library(Biobase)
library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(ggpubr)


################################################################################
# 1. Settings
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
# 2. Files
################################################################################

ESET_FILE <- file.path(
  "data",
  "processed",
  "uams",
  "ESET_uams_all_gep1q_fish_polyig_pcme.rds"
)

OUT_DIR <- file.path(
  "results",
  "figures",
  "Fig1",
  "Fig1A"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


################################################################################
# 3. Load phenotype data
################################################################################

eset <- readRDS(
  ESET_FILE
)

df <- pData(
  eset
)


################################################################################
# 4. Select TTP columns
################################################################################

TIME_COL <- "YearsTTP"
EVENT_COL <- "CensTTP"


################################################################################
# 5. Restrict to RNAS CD138+ NDMM and TT2-TT5
################################################################################

df <- df |>
  filter(
    sample_group == "RNAS_CD138_NDMM",
    tt_on_for_tx1 %in% TT_KEEP
  )


################################################################################
# 6. Prepare survival data
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
# 7. Fit Kaplan-Meier model
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
# 8. Fit Cox model
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
# 9. Create cumulative-event Kaplan-Meier plot
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
# 10. Add hazard-ratio annotation
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
# 11. Combine plot and risk table
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
# 12. Save figure
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
# 13. Save source data
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
# 14. Save statistics
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
# 15. Complete
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