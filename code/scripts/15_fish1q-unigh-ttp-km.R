################################################################################
# STEP 15 — FIGURE S2C
# FISH chromosome 1q × clinical unIGH ARD
# Time to Progression
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

FIGURE_ID <- "Fish1q_x_unIGH"

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


# Cohort
TT_KEEP <- c(
  "TT3a",
  "TT3b",
  "TT4-S",
  "TT4_S-TT3"
)

GEP70_CUTOFF <- 0.66


# FISH 1q
FISH_DIPLOID_VALUE <- 2
FISH_GAIN_MINIMUM <- 3


# Clinical unIGH ARD
UNIGH_ARD_CUT <- 0.29


# TRUE  = combine all FISH 1q >=3-copy cases into one group.
# FALSE = retain separate 1q≥3/unIGH High and Low groups.
COMBINE_1Q3PLUS <- FALSE


# Cox reference group
HR_REFERENCE_GROUP <- "1q2/unIGH High"


# Plot settings
SHOW_LOGRANK_P_VALUE <- TRUE

LOGRANK_P_VALUE_COORD <- c(
  0.5,
  0.85
)

LOGRANK_P_VALUE_SIZE <- 5

SHOW_HR_LABEL <- TRUE
HR_LABEL_TEXT_SIZE <- 3.6

LEGEND_NROW <- 1L
LEGEND_BYROW <- TRUE
LEGEND_TEXT_SIZE <- 13

BREAK_TIME_BY_YEARS <- 2

PLOT_WIDTH_IN <- 9
PLOT_HEIGHT_IN <- 8
PNG_DPI <- 300

SAVE_PAIRWISE_CSV <- TRUE


# Group colors
GROUP_COLORS <- if (COMBINE_1Q3PLUS) {
  
  c(
    "1q2/unIGH Low" = "grey50",
    "1q2/unIGH High" = "steelblue2",
    "1q≥3" = "purple"
  )
  
} else {
  
  c(
    "1q2/unIGH Low" = "grey50",
    "1q2/unIGH High" = "steelblue2",
    "1q≥3/unIGH Low" = "red2",
    "1q≥3/unIGH High" = "purple"
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
    
    unIGH_ARD = as.numeric(
      IG_ARD
    ),
    
    time = as.numeric(
      YearsTTP
    ),
    
    event = as.numeric(
      CensTTP
    )
  )


################################################################################
# 5. RESTRICT TO GEP70 STANDARD RISK
################################################################################

n_gep70_high <- sum(
  is.finite(
    df$GEP70
  ) &
    df$GEP70 > GEP70_CUTOFF
)

df <- df |>
  filter(
    is.finite(
      GEP70
    ),
    GEP70 <= GEP70_CUTOFF
  )


################################################################################
# 6. FINAL ANALYSIS COHORT
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
    
    is.finite(
      FISH1q
    ),
    
    FISH1q == FISH_DIPLOID_VALUE |
      FISH1q >= FISH_GAIN_MINIMUM,
    
    is.finite(
      unIGH_ARD
    )
  )


if (anyDuplicated(
  df$CHIPID
)) {
  stop(
    "Duplicate CHIPIDs found in final cohort."
  )
}


################################################################################
# 7. DEFINE FISH 1q AND unIGH GROUPS
################################################################################

df <- df |>
  mutate(
    FISH1q_group = case_when(
      FISH1q == FISH_DIPLOID_VALUE ~
        "1q2",
      
      FISH1q >= FISH_GAIN_MINIMUM ~
        "1q≥3",
      
      TRUE ~ NA_character_
    ),
    
    unIGH_group = ifelse(
      unIGH_ARD > UNIGH_ARD_CUT,
      "High",
      "Low"
    ),
    
    Group = if (COMBINE_1Q3PLUS) {
      
      case_when(
        FISH1q == 2 &
          unIGH_ARD > UNIGH_ARD_CUT ~
          "1q2/unIGH High",
        
        FISH1q == 2 &
          unIGH_ARD <= UNIGH_ARD_CUT ~
          "1q2/unIGH Low",
        
        FISH1q >= 3 ~
          "1q≥3",
        
        TRUE ~ NA_character_
      )
      
    } else {
      
      case_when(
        FISH1q == 2 &
          unIGH_ARD > UNIGH_ARD_CUT ~
          "1q2/unIGH High",
        
        FISH1q == 2 &
          unIGH_ARD <= UNIGH_ARD_CUT ~
          "1q2/unIGH Low",
        
        FISH1q >= 3 &
          unIGH_ARD > UNIGH_ARD_CUT ~
          "1q≥3/unIGH High",
        
        FISH1q >= 3 &
          unIGH_ARD <= UNIGH_ARD_CUT ~
          "1q≥3/unIGH Low",
        
        TRUE ~ NA_character_
      )
    }
  )


GROUP_LEVELS <- if (COMBINE_1Q3PLUS) {
  
  c(
    "1q2/unIGH High",
    "1q2/unIGH Low",
    "1q≥3"
  )
  
} else {
  
  c(
    "1q2/unIGH High",
    "1q2/unIGH Low",
    "1q≥3/unIGH High",
    "1q≥3/unIGH Low"
  )
}


df$Group <- factor(
  df$Group,
  levels = GROUP_LEVELS
)

df$Group <- droplevels(
  df$Group
)

group_levels_present <- levels(
  df$Group
)

group_counts <- table(
  df$Group
)


if (length(
  group_levels_present
) < 2) {
  stop(
    "At least two nonempty groups are required."
  )
}


if (!HR_REFERENCE_GROUP %in%
    group_levels_present) {
  stop(
    "Reference group is absent: ",
    HR_REFERENCE_GROUP
  )
}


################################################################################
# 8. CHECK PLOT COLORS
################################################################################

palette_used <- unname(
  GROUP_COLORS[
    group_levels_present
  ]
)

if (anyNA(
  palette_used
)) {
  stop(
    "One or more active groups have no color."
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
  survival_object ~ Group,
  data = df
)


################################################################################
# 10. COX MODEL VS REFERENCE GROUP
################################################################################

ordered_levels <- c(
  HR_REFERENCE_GROUP,
  setdiff(
    group_levels_present,
    HR_REFERENCE_GROUP
  )
)

df$Group_hr <- factor(
  as.character(
    df$Group
  ),
  levels = ordered_levels
)

cox_fit <- coxph(
  survival_object ~ Group_hr,
  data = df
)

cox_summary <- summary(
  cox_fit
)

cox_ci <- exp(
  confint(
    cox_fit
  )
)

cox_table <- tibble(
  reference = HR_REFERENCE_GROUP,
  
  comparison = sub(
    "^Group_hr",
    "",
    rownames(
      cox_summary$coefficients
    )
  ),
  
  HR = exp(
    coef(
      cox_fit
    )
  ),
  
  CI_lower_95 = cox_ci[
    ,
    1
  ],
  
  CI_upper_95 = cox_ci[
    ,
    2
  ],
  
  p_value = cox_summary$coefficients[
    ,
    "Pr(>|z|)"
  ]
)


################################################################################
# 11. HR LABEL
################################################################################

format_p <- function(p) {
  
  if (p < 0.001) {
    
    format(
      p,
      scientific = TRUE,
      digits = 2
    )
    
  } else {
    
    formatC(
      p,
      format = "f",
      digits = 3
    )
  }
}


hr_label <- paste(
  c(
    paste0(
      "Reference: ",
      HR_REFERENCE_GROUP
    ),
    
    paste0(
      cox_table$comparison,
      " vs ",
      HR_REFERENCE_GROUP,
      ": HR = ",
      sprintf(
        "%.2f",
        cox_table$HR
      ),
      " (",
      sprintf(
        "%.2f",
        cox_table$CI_lower_95
      ),
      "–",
      sprintf(
        "%.2f",
        cox_table$CI_upper_95
      ),
      ")",
      "\np = ",
      vapply(
        cox_table$p_value,
        format_p,
        character(1)
      )
    )
  ),
  collapse = "\n"
)


################################################################################
# 12. PAIRWISE COX ANALYSES
################################################################################

group_pairs <- combn(
  group_levels_present,
  2,
  simplify = FALSE
)

pairwise_table <- bind_rows(
  lapply(
    group_pairs,
    function(pair) {
      
      reference <- pair[1]
      comparison <- pair[2]
      
      pair_data <- df |>
        filter(
          Group %in% c(
            reference,
            comparison
          )
        ) |>
        mutate(
          Group_pair = factor(
            as.character(
              Group
            ),
            levels = c(
              reference,
              comparison
            )
          )
        )
      
      fit <- coxph(
        Surv(
          time,
          event
        ) ~ Group_pair,
        data = pair_data
      )
      
      fit_summary <- summary(
        fit
      )
      
      ci <- exp(
        confint(
          fit
        )
      )[1, ]
      
      tibble(
        reference = reference,
        comparison = comparison,
        
        n_reference = sum(
          pair_data$Group_pair ==
            reference
        ),
        
        n_comparison = sum(
          pair_data$Group_pair ==
            comparison
        ),
        
        events_reference = sum(
          pair_data$event[
            pair_data$Group_pair ==
              reference
          ] == 1
        ),
        
        events_comparison = sum(
          pair_data$event[
            pair_data$Group_pair ==
              comparison
          ] == 1
        ),
        
        HR = exp(
          coef(
            fit
          )
        )[1],
        
        CI_lower_95 = ci[1],
        
        CI_upper_95 = ci[2],
        
        p_value = fit_summary$coefficients[
          1,
          "Pr(>|z|)"
        ]
      )
    }
  )
)


################################################################################
# 13. SOURCE DATA
################################################################################

analysis_cohort_export <- df |>
  transmute(
    CHIPID,
    
    TT_Protocol = tt_on_for_tx1,
    
    FISH1qCopies = FISH1q,
    
    FISH1q_Group = FISH1q_group,
    
    unIGH_ARD,
    
    unIGH_ARD_Group = unIGH_group,
    
    GEP70,
    
    Group = as.character(
      Group
    ),
    
    Time_Years = time,
    
    Event = event,
    
    Endpoint = "TTP"
  )


group_summary <- df |>
  group_by(
    Group
  ) |>
  summarise(
    n = n(),
    events = sum(
      event == 1
    ),
    censored = sum(
      event == 0
    ),
    .groups = "drop"
  ) |>
  mutate(
    Group = as.character(
      Group
    )
  )


write_csv(
  analysis_cohort_export,
  file.path(
    SOURCE_DATA_DIR,
    paste0(
      FIGURE_ID,
      "_FISH1q_x_unIGH_ARD_TTP_analysis-cohort.csv"
    )
  )
)


write_csv(
  group_summary,
  file.path(
    SOURCE_DATA_DIR,
    paste0(
      FIGURE_ID,
      "_FISH1q_x_unIGH_ARD_TTP_group-summary.csv"
    )
  )
)


write_csv(
  cox_table,
  file.path(
    SOURCE_DATA_DIR,
    paste0(
      FIGURE_ID,
      "_FISH1q_x_unIGH_ARD_TTP_cox-vs-reference.csv"
    )
  )
)


if (SAVE_PAIRWISE_CSV) {
  
  write_csv(
    pairwise_table,
    file.path(
      SOURCE_DATA_DIR,
      paste0(
        FIGURE_ID,
        "_FISH1q_x_unIGH_ARD_TTP_pairwise-cox.csv"
      )
    )
  )
}


################################################################################
# 14. PLOT
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


plot_title <-
  "Time to Progression by FISH 1q Copy × unIGH"


plot_subtitle <- paste0(
  "Total Therapy: TT3a, TT3b, TT4-S",
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
  
  pval = SHOW_LOGRANK_P_VALUE,
  
  pval.coord = LOGRANK_P_VALUE_COORD,
  
  pval.size = LOGRANK_P_VALUE_SIZE,
  
  risk.table = TRUE,
  
  risk.table.col = "strata",
  
  risk.table.y.text = TRUE,
  
  risk.table.y.text.col = TRUE,
  
  risk.table.height = 0.22,
  
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
      ),
      
      plot.subtitle = element_text(
        hjust = 0.5,
        size = 18
      )
    ),
  
  tables.theme = theme_cleantable()
)


km_plot$plot <- km_plot$plot +
  labs(
    subtitle = plot_subtitle
  ) +
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
    limits = c(
      0,
      1
    ),
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


if (SHOW_HR_LABEL) {
  
  plot_build <- ggplot_build(
    km_plot$plot
  )
  
  x_range <-
    plot_build$layout$panel_params[[1]]$x.range
  
  y_range <-
    plot_build$layout$panel_params[[1]]$y.range
  
  km_plot$plot <- km_plot$plot +
    annotate(
      geom = "text",
      
      x = x_range[2],
      
      y = y_range[2] -
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


km_plot$table <- km_plot$table +
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
# 15. SAVE FIGURE
################################################################################

png_file <- file.path(
  FIGURE_DIR,
  paste0(
    FIGURE_ID,
    "_FISH1q_x_unIGH_ARD_TTP_KM.png"
  )
)

pdf_file <- file.path(
  FIGURE_DIR,
  paste0(
    FIGURE_ID,
    "_FISH1q_x_unIGH_ARD_TTP_KM.pdf"
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
# 16. SUMMARY
################################################################################

message(
  "\nStep 15 complete — Figure S2C",
  
  "\nFinal cohort: ",
  nrow(
    df
  ),
  
  "\nEvents: ",
  sum(
    df$event == 1
  ),
  
  "\nCombine 1q≥3 groups: ",
  COMBINE_1Q3PLUS,
  
  "\nunIGH ARD cutoff: ",
  sprintf(
    "%.2f",
    UNIGH_ARD_CUT
  ),
  
  "\nGEP70-high excluded: ",
  n_gep70_high,
  
  "\nReference group: ",
  HR_REFERENCE_GROUP,
  
  "\nGroups: ",
  paste(
    names(
      group_counts
    ),
    unname(
      group_counts
    ),
    sep = "=",
    collapse = "; "
  ),
  
  "\nPNG: ",
  png_file,
  
  "\nPDF: ",
  pdf_file,
  
  "\n"
)