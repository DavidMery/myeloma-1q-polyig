################################################################################
# STEP 14 — FIGURE S2B
# Kaplan-Meier Time to Progression by clinical unIGH ARD
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

FIGURE_ID <- "unIGH"

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

UNIGH_ARD_CUT <- 0.15

USE_OPTIMAL_CUTPOINT <- TRUE
CUTPOINT_MINPROP <- 0.15
SAVE_CUTPOINT_SCAN <- TRUE

BREAK_TIME_BY_YEARS <- 2

PLOT_WIDTH_IN <- 8.5
PLOT_HEIGHT_IN <- 6
PNG_DPI <- 300

COLOR_LOWER_ARD <- "red3"
COLOR_HIGHER_ARD <- "steelblue3"

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
    
    unIGH_ARD = as.numeric(
      IG_ARD
    ),
    
    FISH1q_2copies = as.numeric(
      FISH_1q_2copies
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
      unIGH_ARD
    ),
    
    is.finite(
      FISH1q_2copies
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
# 6. FIND unIGH ARD CUTPOINT
################################################################################

cutpoint_source <- "manual"

if (USE_OPTIMAL_CUTPOINT) {
  
  score <- df$unIGH_ARD
  
  x <- sort(
    unique(
      score
    )
  )
  
  candidate_cutpoints <- (
    x[-1] +
      x[-length(x)]
  ) / 2
  
  minimum_n <- ceiling(
    CUTPOINT_MINPROP *
      length(
        score
      )
  )
  
  scan <- lapply(
    candidate_cutpoints,
    function(cutpoint) {
      
      higher_group <- score > cutpoint
      
      n_lower <- sum(
        !higher_group
      )
      
      n_higher <- sum(
        higher_group
      )
      
      if (
        n_lower < minimum_n ||
        n_higher < minimum_n
      ) {
        return(
          NULL
        )
      }
      
      fit <- tryCatch(
        coxph(
          Surv(
            df$time,
            df$event
          ) ~ factor(
            higher_group
          )
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
      
      data.frame(
        cutpoint = cutpoint,
        n_lower = n_lower,
        n_higher = n_higher,
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
      "No valid unIGH ARD cutpoint was found."
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
  
  UNIGH_ARD_CUT <- scan$cutpoint[1]
  
  cutpoint_source <- "optimal_hr"
  
  if (SAVE_CUTPOINT_SCAN) {
    
    write_csv(
      scan,
      file.path(
        SOURCE_DATA_DIR,
        paste0(
          FIGURE_ID,
          "_unIGH_ARD_TTP_cutpoint-scan.csv"
        )
      )
    )
  }
}


################################################################################
# 7. BUILD GROUPS
################################################################################

cutpoint_display <- formatC(
  UNIGH_ARD_CUT,
  format = "f",
  digits = 2
)

lower_ard_label <- paste0(
  "unIGH ≤ ",
  cutpoint_display
)

higher_ard_label <- paste0(
  "unIGH > ",
  cutpoint_display
)

df$unIGH_ARD_group <- factor(
  ifelse(
    df$unIGH_ARD > UNIGH_ARD_CUT,
    higher_ard_label,
    lower_ard_label
  ),
  levels = c(
    lower_ard_label,
    higher_ard_label
  )
)

group_counts <- table(
  df$unIGH_ARD_group
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
    "Both unIGH ARD groups must contain samples."
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
  survival_object ~ unIGH_ARD_group,
  data = df
)

cox_fit <- coxph(
  survival_object ~ unIGH_ARD_group,
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

hr_label <- paste0(
  "Higher vs lower unIGH: HR = ",
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
# 9. SOURCE DATA
################################################################################

analysis_cohort_export <- df |>
  transmute(
    CHIPID,
    GEP70,
    TT_Protocol = tt_on_for_tx1,
    FISH_1q_2copies,
    unIGH_ARD,
    unIGH_ARD_Group = as.character(
      unIGH_ARD_group
    ),
    Time_Years = time,
    Event = event,
    Endpoint = "TTP"
  )

cox_statistics_export <- tibble(
  Figure = FIGURE_ID,
  Endpoint = "TTP",
  Comparison = "Higher vs lower unIGH ARD",
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
  unIGH_ARD_group = names(
    table(
      df$unIGH_ARD_group
    )
  ),
  
  n = as.integer(
    table(
      df$unIGH_ARD_group
    )
  )
)

cutpoint_summary_export <- tibble(
  Figure = FIGURE_ID,
  Endpoint = "TTP",
  score_column = "IG_ARD",
  cutpoint = UNIGH_ARD_CUT,
  cutpoint_source = cutpoint_source,
  minimum_group_proportion = CUTPOINT_MINPROP,
  n_total = nrow(
    df
  ),
  lower_group = lower_ard_label,
  higher_group = higher_ard_label,
  n_lower = sum(
    df$unIGH_ARD_group ==
      lower_ard_label
  ),
  n_higher = sum(
    df$unIGH_ARD_group ==
      higher_ard_label
  ),
  HR_higher_vs_lower = hazard_ratio,
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
      "_unIGH_ARD_TTP_analysis-cohort.csv"
    )
  )
)

write_csv(
  cox_statistics_export,
  file.path(
    SOURCE_DATA_DIR,
    paste0(
      FIGURE_ID,
      "_unIGH_ARD_TTP_cox-statistics.csv"
    )
  )
)

write_csv(
  group_counts_export,
  file.path(
    SOURCE_DATA_DIR,
    paste0(
      FIGURE_ID,
      "_unIGH_ARD_TTP_group-counts.csv"
    )
  )
)

write_csv(
  cutpoint_summary_export,
  file.path(
    SOURCE_DATA_DIR,
    paste0(
      FIGURE_ID,
      "_unIGH_ARD_TTP_cutpoint-summary.csv"
    )
  )
)


################################################################################
# 10. FIGURE
################################################################################

max_time <- max(
  df$time,
  na.rm = TRUE
)

x_max <- ceiling(
  max_time /
    BREAK_TIME_BY_YEARS
) * BREAK_TIME_BY_YEARS

plot_title <- paste0(
  "Time to Progress by Uninvolved Immunoglobulin (unIGH)\n",
  "Standard Risk (GEP70); TT2–TT4 (n = ",
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
    COLOR_LOWER_ARD,
    COLOR_HIGHER_ARD
  ),
  
  legend.title = "Clinical unIGH",
  
  legend.labs = c(
    lower_ard_label,
    higher_ard_label
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
  
  ylab = "Cumulative Progress Probability",
  
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
  
  tables.theme = survminer::theme_cleantable()
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
# 11. SAVE
################################################################################

png_file <- file.path(
  FIGURE_DIR,
  paste0(
    FIGURE_ID,
    "_unIGH_ARD_TTP_KM.png"
  )
)

pdf_file <- file.path(
  FIGURE_DIR,
  paste0(
    FIGURE_ID,
    "_unIGH_ARD_TTP_KM.pdf"
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
  device = "pdf",
  bg = "white"
)


################################################################################
# 12. SUMMARY
################################################################################

message(
  "\nStep 14 complete — Figure S2B",
  "\nFinal cohort: ",
  nrow(
    df
  ),
  "\nEvents: ",
  sum(
    df$event == 1
  ),
  "\nCutpoint: ",
  sprintf(
    "%.4f",
    UNIGH_ARD_CUT
  ),
  "\nCutpoint source: ",
  cutpoint_source,
  "\nLower unIGH: ",
  unname(
    group_counts[
      lower_ard_label
    ]
  ),
  "\nHigher unIGH: ",
  unname(
    group_counts[
      higher_ard_label
    ]
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