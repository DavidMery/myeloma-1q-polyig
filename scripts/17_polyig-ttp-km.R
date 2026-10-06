################################################################################
# STEP 17 — FIGURE S3B
# PolyIG Kaplan-Meier analysis of Time to Progression
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
# 1. SETTINGS
################################################################################

FIGURE_ID <- "PolyIG"

ESET_PATH <- file.path(
  "data",
  "processed",
  "uams",
  "ESET_uams_all_gep1q_fish_polyig_pcme.rds"
)

FIGURE_DIR <- file.path(
  "results",
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


# PolyIG cutoff
# Set FALSE to use the prespecified cutoff of 11.
USE_OPTIMAL_CUTPOINT <- TRUE

POLYIG_CUT <- 11

CUTPOINT_MINPROP <- 0.15


# Plot settings
BREAK_TIME_BY_YEARS <- 2

PLOT_WIDTH_IN <- 8
PLOT_HEIGHT_IN <- 6
PNG_DPI <- 300

COLOR_POLYIG_LOW <- "red3"
COLOR_POLYIG_HIGH <- "steelblue3"


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
    
    PolyIG = as.numeric(
      PolyIG_Score
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
      PolyIG
    ),
    
    is.finite(
      time
    ),
    
    time >= 0,
    
    event %in% c(
      0,
      1
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
# 6. FIND POLYIG CUTPOINT
################################################################################

cutpoint_source <- "prespecified"

if (USE_OPTIMAL_CUTPOINT) {
  
  score_values <- sort(
    unique(
      df$PolyIG
    )
  )
  
  candidate_cutpoints <- (
    score_values[-1] +
      score_values[-length(score_values)]
  ) / 2
  
  minimum_n <- ceiling(
    CUTPOINT_MINPROP *
      nrow(
        df
      )
  )
  
  scan <- lapply(
    candidate_cutpoints,
    function(cutpoint) {
      
      high <- df$PolyIG > cutpoint
      
      n_low <- sum(
        !high
      )
      
      n_high <- sum(
        high
      )
      
      if (
        n_low < minimum_n ||
        n_high < minimum_n
      ) {
        return(
          NULL
        )
      }
      
      fit <- tryCatch(
        coxph(
          Surv(
            time,
            event
          ) ~ high,
          data = df
        ),
        error = function(e) {
          NULL
        }
      )
      
      if (is.null(
        fit
      )) {
        return(
          NULL
        )
      }
      
      fit_summary <- summary(
        fit
      )
      
      hr <- exp(
        coef(
          fit
        )
      )[1]
      
      ci <- exp(
        confint(
          fit
        )
      )[1, ]
      
      p <- fit_summary$coef[
        1,
        "Pr(>|z|)"
      ]
      
      tibble(
        cutpoint = cutpoint,
        n_low = n_low,
        n_high = n_high,
        HR = hr,
        CI_lower_95 = ci[1],
        CI_upper_95 = ci[2],
        p_value = p,
        score_metric = abs(
          log(
            hr
          )
        )
      )
    }
  )
  
  scan <- bind_rows(
    scan
  )
  
  if (!nrow(
    scan
  )) {
    stop(
      "No valid PolyIG cutpoints were available."
    )
  }
  
  scan <- scan |>
    arrange(
      desc(
        score_metric
      ),
      p_value,
      cutpoint
    )
  
  POLYIG_CUT <- scan$cutpoint[1]
  
  cutpoint_source <- "optimal_hr"
  
  write_csv(
    scan,
    file.path(
      SOURCE_DATA_DIR,
      paste0(
        FIGURE_ID,
        "_PolyIG_TTP_cutpoint-scan.csv"
      )
    )
  )
}


################################################################################
# 7. DEFINE POLYIG GROUPS
################################################################################

df <- df |>
  mutate(
    PolyIG_group = ifelse(
      PolyIG > POLYIG_CUT,
      "PolyIG High",
      "PolyIG Low"
    ),
    
    PolyIG_group = factor(
      PolyIG_group,
      levels = c(
        "PolyIG Low",
        "PolyIG High"
      )
    )
  )


group_counts <- table(
  df$PolyIG_group
)

if (
  length(
    group_counts
  ) != 2 ||
  any(
    group_counts == 0
  )
) {
  stop(
    "Both PolyIG groups must contain samples."
  )
}


################################################################################
# 8. SURVIVAL ANALYSIS
################################################################################

survival_object <- Surv(
  time = df$time,
  event = df$event
)

survival_fit <- survfit(
  survival_object ~ PolyIG_group,
  data = df
)

cox_fit <- coxph(
  survival_object ~ PolyIG_group,
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
# 9. HR LABEL
################################################################################

hr_label <- paste0(
  "PolyIG High vs PolyIG Low: HR = ",
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
# 10. SOURCE DATA
################################################################################

analysis_cohort_export <- df |>
  transmute(
    CHIPID,
    
    GEP70,
    
    TT_Protocol = tt_on_for_tx1,
    
    PolyIG,
    
    PolyIG_Group = as.character(
      PolyIG_group
    ),
    
    Time_Years = time,
    
    Event = event,
    
    Endpoint = "TTP"
  )


cox_statistics_export <- tibble(
  Figure = FIGURE_ID,
  
  Step = 17L,
  
  Endpoint = "TTP",
  
  Comparison = "PolyIG High vs PolyIG Low",
  
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
  PolyIG_group = names(
    table(
      df$PolyIG_group
    )
  ),
  
  n = as.integer(
    table(
      df$PolyIG_group
    )
  )
)


cutpoint_summary_export <- tibble(
  Figure = FIGURE_ID,
  
  score_column = "PolyIG_Score",
  
  cutpoint = POLYIG_CUT,
  
  cutpoint_source = cutpoint_source,
  
  minimum_group_proportion = if (
    USE_OPTIMAL_CUTPOINT
  ) {
    CUTPOINT_MINPROP
  } else {
    NA_real_
  },
  
  n_total = nrow(
    df
  ),
  
  n_low = sum(
    df$PolyIG_group == "PolyIG Low"
  ),
  
  n_high = sum(
    df$PolyIG_group == "PolyIG High"
  ),
  
  HR_high_vs_low = hazard_ratio,
  
  CI_lower_95 = confidence_interval[1],
  
  CI_upper_95 = confidence_interval[2],
  
  p_value = cox_p_value
)


write_csv(
  analysis_cohort_export,
  file.path(
    SOURCE_DATA_DIR,
    paste0(
      FIGURE_ID,
      "_PolyIG_TTP_analysis-cohort.csv"
    )
  )
)


write_csv(
  cox_statistics_export,
  file.path(
    SOURCE_DATA_DIR,
    paste0(
      FIGURE_ID,
      "_PolyIG_TTP_cox-statistics.csv"
    )
  )
)


write_csv(
  group_counts_export,
  file.path(
    SOURCE_DATA_DIR,
    paste0(
      FIGURE_ID,
      "_PolyIG_TTP_group-counts.csv"
    )
  )
)


write_csv(
  cutpoint_summary_export,
  file.path(
    SOURCE_DATA_DIR,
    paste0(
      FIGURE_ID,
      "_PolyIG_TTP_cutpoint-summary.csv"
    )
  )
)


################################################################################
# 11. PLOT
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
  "Time to Progression by PolyIG (GEP)\n",
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
    COLOR_POLYIG_LOW,
    COLOR_POLYIG_HIGH
  ),
  
  legend.title = "PolyIG group",
  
  legend.labs = c(
    "PolyIG Low",
    "PolyIG High"
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
# 12. SAVE FIGURE
################################################################################

png_file <- file.path(
  FIGURE_DIR,
  paste0(
    FIGURE_ID,
    "_PolyIG_TTP_KM.png"
  )
)

pdf_file <- file.path(
  FIGURE_DIR,
  paste0(
    FIGURE_ID,
    "_PolyIG_TTP_KM.pdf"
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
# 13. SUMMARY
################################################################################

message(
  "\nStep 17 complete",
  
  "\nFinal cohort: ",
  nrow(
    df
  ),
  
  "\nEvents: ",
  sum(
    df$event == 1
  ),
  
  "\nPolyIG cutoff: ",
  sprintf(
    "%.4f",
    POLYIG_CUT
  ),
  
  "\nCutpoint source: ",
  cutpoint_source,
  
  "\nPolyIG Low: ",
  sum(
    df$PolyIG_group == "PolyIG Low"
  ),
  
  "\nPolyIG High: ",
  sum(
    df$PolyIG_group == "PolyIG High"
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