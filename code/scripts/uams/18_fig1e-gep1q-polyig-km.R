################################################################################
# STEP 18 — FIGURE 1E LASSO VERSION
#
# Purpose:
#   1. Load the UAMS ExpressionSet containing GEP1q and PolyIG scores.
#   2. Select baseline RNAS CD138+ NDMM samples treated on TT2–TT5.
#   3. Restrict to patients with complete TTP, GEP1q, PolyIG, and GEP70 data.
#   4. Define GEP70-standard-risk molecular groups.
#   5. Generate cumulative progression curves.
#   6. Calculate Cox hazard ratios and all pairwise Cox comparisons.
#
# Required first:
#   Run through code/run_pipeline.R so DATA_ROOT and RESULTS_ROOT are defined.
################################################################################

suppressPackageStartupMessages({
  library(Biobase)
  library(survival)
  library(survminer)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
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
# 1. SETTINGS
################################################################################

ESET_PATH <- file.path(
  DATA_ROOT,
  "processed",
  "ESET_uams_all_gep1q_fish_polyig_pcme.rds"
)

FIGURE_DIR <- file.path(
  RESULTS_ROOT,
  "figures",
  "Fig1",
  "Fig1E"
)

POLYIG_CUT <- 11
GEP70_THRESHOLD <- 0.66

COMBINE_1Q3PLUS <- TRUE
INCLUDE_GEP70_COMPARATOR <- TRUE
FILTER_TT_PROTOCOLS <- TRUE

TT_KEEP <- c(
  "TT2_NoThal",
  "TT2_Thal",
  "TT3a",
  "TT3b",
  "TT4_S-TT3",
  "TT5"
)

BREAK_TIME_BY_YEARS <- 2

PLOT_WIDTH_IN <- 10.5
PLOT_HEIGHT_IN <- 8
PNG_DPI <- 300


################################################################################
# 2. HELPER
################################################################################

fmt_p <- function(p) {
  
  if (!is.finite(p)) {
    
    return(
      "NA"
    )
  }
  
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
# 3. LOAD DATA
################################################################################

if (!file.exists(ESET_PATH)) {
  
  stop(
    "Cannot find ExpressionSet:\n",
    ESET_PATH
  )
}

dir.create(
  FIGURE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

eset <- readRDS(
  ESET_PATH
)

if (!methods::is(
  eset,
  "ExpressionSet"
)) {
  
  stop(
    "ESET_PATH did not contain an ExpressionSet."
  )
}

df <- tibble::as_tibble(
  as.data.frame(
    Biobase::pData(
      eset
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
)


################################################################################
# 4. REQUIRED COLUMNS
################################################################################

REQUIRED_COLUMNS <- c(
  "sample_group",
  "TT_On_for_Tx1",
  "GEP1qcopy_lasso",
  "PolyIG_Score",
  "GEP70",
  "YearsTTP",
  "CensTTP"
)

missing_columns <- setdiff(
  REQUIRED_COLUMNS,
  colnames(df)
)

if (length(
  missing_columns
) > 0) {
  
  stop(
    "Missing required column(s): ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}


################################################################################
# 5. COHORT
################################################################################

# Baseline CD138-selected NDMM
df <- df |>
  filter(
    sample_group == "RNAS_CD138_NDMM"
  )

# Total Therapy 2–5
if (FILTER_TT_PROTOCOLS) {
  
  df <- df |>
    filter(
      TT_On_for_Tx1 %in% TT_KEEP
    )
}


################################################################################
# 6. SURVIVAL + MOLECULAR VARIABLES
################################################################################

df <- df |>
  mutate(
    time = YearsTTP,
    event = CensTTP,
    GEP1q_class = GEP1qcopy_lasso,
    PolyIG = PolyIG_Score
  ) |>
  filter(
    is.finite(
      time
    ),
    
    time >= 0,
    
    event %in% c(
      0,
      1
    ),
    
    GEP1q_class %in% c(
      "1q2",
      "1q3plus"
    ),
    
    is.finite(
      PolyIG
    ),
    
    is.finite(
      GEP70
    )
  ) |>
  mutate(
    GEP70_High = GEP70 > GEP70_THRESHOLD,
    
    GEP70_Low = GEP70 <= GEP70_THRESHOLD,
    
    PolyIG_group = if_else(
      PolyIG > POLYIG_CUT,
      "High",
      "Low"
    )
  )

if (!INCLUDE_GEP70_COMPARATOR) {
  
  df <- df |>
    filter(
      GEP70_Low
    )
}


################################################################################
# 7. COMPOSITE GROUPS
################################################################################

LABEL_1Q2_LOW <- "CN2/PolyIG Low"
LABEL_1Q2_HIGH <- "CN2/PolyIG High"

LABEL_1Q3PLUS_LOW <- "CN>=3/PolyIG Low"
LABEL_1Q3PLUS_HIGH <- "CN>=3/PolyIG High"

GEP70_HIGH_LABEL <- "GEP70 > 0.66"


df <- df |>
  mutate(
    
    Group = if (COMBINE_1Q3PLUS) {
      
      case_when(
        
        INCLUDE_GEP70_COMPARATOR &
          GEP70_High ~
          GEP70_HIGH_LABEL,
        
        GEP70_Low &
          GEP1q_class == "1q2" &
          PolyIG_group == "Low" ~
          LABEL_1Q2_LOW,
        
        GEP70_Low &
          GEP1q_class == "1q2" &
          PolyIG_group == "High" ~
          LABEL_1Q2_HIGH,
        
        GEP70_Low &
          GEP1q_class == "1q3plus" ~
          "CN>=3",
        
        TRUE ~
          NA_character_
      )
      
    } else {
      
      case_when(
        
        INCLUDE_GEP70_COMPARATOR &
          GEP70_High ~
          GEP70_HIGH_LABEL,
        
        GEP70_Low &
          GEP1q_class == "1q2" &
          PolyIG_group == "Low" ~
          LABEL_1Q2_LOW,
        
        GEP70_Low &
          GEP1q_class == "1q2" &
          PolyIG_group == "High" ~
          LABEL_1Q2_HIGH,
        
        GEP70_Low &
          GEP1q_class == "1q3plus" &
          PolyIG_group == "Low" ~
          LABEL_1Q3PLUS_LOW,
        
        GEP70_Low &
          GEP1q_class == "1q3plus" &
          PolyIG_group == "High" ~
          LABEL_1Q3PLUS_HIGH,
        
        TRUE ~
          NA_character_
      )
    }
  ) |>
  filter(
    !is.na(
      Group
    )
  )


################################################################################
# 8. GROUP ORDER
################################################################################

GROUP_LEVELS <- if (COMBINE_1Q3PLUS) {
  
  c(
    LABEL_1Q2_HIGH,
    LABEL_1Q2_LOW,
    "CN>=3",
    GEP70_HIGH_LABEL
  )
  
} else {
  
  c(
    LABEL_1Q2_HIGH,
    LABEL_1Q2_LOW,
    LABEL_1Q3PLUS_HIGH,
    LABEL_1Q3PLUS_LOW,
    GEP70_HIGH_LABEL
  )
}

df$Group <- factor(
  df$Group,
  levels = GROUP_LEVELS
)

df$Group <- droplevels(
  df$Group
)

if (nlevels(
  df$Group
) < 2) {
  
  stop(
    "At least two groups are required."
  )
}


################################################################################
# 9. SURVIVAL FIT
################################################################################

survival_object <- survival::Surv(
  df$time,
  df$event
)

survival_fit <- survival::survfit(
  survival_object ~ Group,
  data = df
)


################################################################################
# 10. OVERALL LOG-RANK TEST
################################################################################

logrank_test <- survival::survdiff(
  survival_object ~ Group,
  data = df
)

overall_p_value <- 1 - pchisq(
  logrank_test$chisq,
  df = length(
    logrank_test$n
  ) - 1
)


################################################################################
# 11. COX MODEL
# Favorable CN2/PolyIG-high is the reference group.
################################################################################

HR_REFERENCE_GROUP <- LABEL_1Q2_HIGH

hr_data <- df |>
  mutate(
    Group_hr = factor(
      as.character(
        Group
      ),
      levels = c(
        HR_REFERENCE_GROUP,
        setdiff(
          levels(
            Group
          ),
          HR_REFERENCE_GROUP
        )
      )
    )
  )

cox_fit <- survival::coxph(
  survival::Surv(
    time,
    event
  ) ~ Group_hr,
  data = hr_data
)

cox_summary <- summary(
  cox_fit
)

hr <- exp(
  coef(
    cox_fit
  )
)

ci <- exp(
  confint(
    cox_fit
  )
)

cox_p <- cox_summary$coefficients[
  ,
  "Pr(>|z|)"
]

hr_table <- tibble::tibble(
  Reference = HR_REFERENCE_GROUP,
  
  Comparison = sub(
    "^Group_hr",
    "",
    names(
      hr
    )
  ),
  
  HR = as.numeric(
    hr
  ),
  
  CI_lower_95 = as.numeric(
    ci[
      ,
      1
    ]
  ),
  
  CI_upper_95 = as.numeric(
    ci[
      ,
      2
    ]
  ),
  
  P_value = as.numeric(
    cox_p
  )
)


################################################################################
# 12. ALL PAIRWISE COX HAZARD RATIOS
################################################################################

group_levels_present <- levels(
  df$Group
)

pairwise_groups <- utils::combn(
  group_levels_present,
  2,
  simplify = FALSE
)

pairwise_hr_table <- dplyr::bind_rows(
  
  lapply(
    pairwise_groups,
    
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
      
      pair_fit <- survival::coxph(
        survival::Surv(
          time,
          event
        ) ~ Group_pair,
        data = pair_data
      )
      
      pair_summary <- summary(
        pair_fit
      )
      
      pair_hr <- exp(
        coef(
          pair_fit
        )
      )[1]
      
      pair_ci <- exp(
        confint(
          pair_fit
        )
      )[1, ]
      
      pair_p <- pair_summary$coefficients[
        1,
        "Pr(>|z|)"
      ]
      
      tibble::tibble(
        
        Reference = reference,
        
        Comparison = comparison,
        
        n_Reference = sum(
          pair_data$Group_pair ==
            reference
        ),
        
        n_Comparison = sum(
          pair_data$Group_pair ==
            comparison
        ),
        
        Events_Reference = sum(
          pair_data$event[
            pair_data$Group_pair ==
              reference
          ] == 1
        ),
        
        Events_Comparison = sum(
          pair_data$event[
            pair_data$Group_pair ==
              comparison
          ] == 1
        ),
        
        HR = pair_hr,
        
        CI_lower_95 = pair_ci[1],
        
        CI_upper_95 = pair_ci[2],
        
        P_value = pair_p
      )
    }
  )
)


################################################################################
# 13. SAVE COX TABLES
################################################################################

write.csv(
  hr_table,
  file.path(
    FIGURE_DIR,
    "Fig1E_LASSO_HR_table.csv"
  ),
  row.names = FALSE
)

write.csv(
  pairwise_hr_table,
  file.path(
    FIGURE_DIR,
    "Fig1E_LASSO_all_pairwise_HR_table.csv"
  ),
  row.names = FALSE
)


################################################################################
# 14. FIGURE SETTINGS
################################################################################

x_max <- ceiling(
  max(
    df$time,
    na.rm = TRUE
  ) / BREAK_TIME_BY_YEARS
) * BREAK_TIME_BY_YEARS

if (!is.finite(
  x_max
) ||
x_max <= 0) {
  
  stop(
    "Could not calculate a valid x-axis maximum."
  )
}

plot_title <- paste0(
  "Time to Progression by Penalized GEP1q Copy × PolyIG\n",
  "Standard Risk (GEP70); TT2–TT5",
  " (n = ",
  nrow(
    df
  ),
  ")"
)

y_axis_label <- "Cumulative Progression Probability"


################################################################################
# 15. PLOT COLORS
################################################################################

GROUP_COLORS <- if (COMBINE_1Q3PLUS) {
  
  c(
    "grey50",
    "steelblue2",
    "red",
    "black"
  )
  
} else {
  
  c(
    "grey50",
    "steelblue2",
    "red2",
    "purple",
    "black"
  )
}

GROUP_COLOR_LABELS <- if (COMBINE_1Q3PLUS) {
  
  c(
    LABEL_1Q2_LOW,
    LABEL_1Q2_HIGH,
    "CN>=3",
    GEP70_HIGH_LABEL
  )
  
} else {
  
  c(
    LABEL_1Q2_LOW,
    LABEL_1Q2_HIGH,
    LABEL_1Q3PLUS_LOW,
    LABEL_1Q3PLUS_HIGH,
    GEP70_HIGH_LABEL
  )
}

palette_used <- GROUP_COLORS[
  match(
    levels(
      df$Group
    ),
    GROUP_COLOR_LABELS
  )
]

if (any(
  is.na(
    palette_used
  )
)) {
  
  stop(
    "Missing color for group(s): ",
    paste(
      levels(
        df$Group
      )[
        is.na(
          palette_used
        )
      ],
      collapse = ", "
    )
  )
}


################################################################################
# 16. KAPLAN–MEIER / CUMULATIVE PROGRESSION PLOT
################################################################################

km_plot <- survminer::ggsurvplot(
  
  fit = survival_fit,
  
  data = df,
  
  fun = "event",
  
  pval = FALSE,
  
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
  
  legend.title = "GEP1q × PolyIG group",
  
  legend.labs = levels(
    df$Group
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
  
  ylab = y_axis_label,
  
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


################################################################################
# 17. ADD OVERALL LOG-RANK P VALUE
################################################################################

km_plot$plot <- km_plot$plot +
  
  ggplot2::annotate(
    geom = "text",
    
    x = x_max * 0.98,
    
    y = 0.96,
    
    label = paste0(
      "Overall log-rank P = ",
      fmt_p(
        overall_p_value
      )
    ),
    
    hjust = 1,
    
    vjust = 1,
    
    size = 4.5
  )


################################################################################
# 18. RISK TABLE FORMATTING
################################################################################

km_plot$table <- km_plot$table +
  
  ggplot2::guides(
    color = "none",
    colour = "none",
    fill = "none",
    shape = "none"
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


################################################################################
# 19. COMBINE PLOT + RISK TABLE
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
# 20. SAVE FIGURE
################################################################################

png_file <- file.path(
  FIGURE_DIR,
  "Fig1E_LASSO_GEP1q_lasso_x_PolyIG_TTP_KM.png"
)

pdf_file <- file.path(
  FIGURE_DIR,
  "Fig1E_LASSO_GEP1q_lasso_x_PolyIG_TTP_KM.pdf"
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
# 21. COMPLETE
################################################################################

message(
  "\nStep 18 complete.",
  "\nEndpoint = TTP",
  "\nPolyIG column = PolyIG_Score",
  "\nGEP1q column = GEP1qcopy_lasso",
  "\nCombine >=3-copy 1q = ",
  COMBINE_1Q3PLUS,
  "\nFinal n = ",
  nrow(
    df
  ),
  "\nEvents = ",
  sum(
    df$event == 1
  ),
  "\nPNG = ",
  png_file,
  "\nPDF = ",
  pdf_file
)