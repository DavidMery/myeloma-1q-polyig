################################################################################
# STEP 26B — FIGURE 3B
# TIME TO PROGRESSION BY GEP1q COPY NUMBER × POLYIG
################################################################################

library(Biobase)
library(survival)
library(survminer)
library(dplyr)
library(readr)
library(ggplot2)
library(ggpubr)

################################################################################
# 1. CONFIGURATION
################################################################################

FIGURE_ID <- "Fig3B_Lasso"
ENDPOINT <- "TTP"

ESET_PATH <- file.path(
  "data", "processed", "uams",
  "ESET_uams_all_gep1q_fish_polyig_pcme_binaries_lasso_groups_rnbx-matched.rds"
)

FIGURE_DIR <- file.path("results", "figures", "Fig3", FIGURE_ID)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

# Exact columns
SAMPLE_ID_COL <- "chipid"
SAMPLE_GROUP_COL <- "sample_group"
TT_COL <- "TT_On_for_Tx1"

GEP1Q_COL <- "GEP1qcopy_lasso"
POLYIG_COL <- "PolyIG_Score"
GEP70_COL <- "GEP70"

TIME_COL <- "YearsTTP"
EVENT_COL <- "CensTTP"

POLYIG_CUT <- 11
USE_OPTIMAL_CUTPOINT <- FALSE
OPTIMAL_CUT_METHOD <- "hr"
CUTPOINT_MINPROP <- 0.15

REQUIRE_BASELINE_RNAS_NDMM <- TRUE
FILTER_TT_PROTOCOLS <- TRUE

TT_KEEP <- c(
  "TT3a",
  "TT3b",
  "TT4_S-TT3"
)

GEP70_THRESHOLD <- 0.66
INCLUDE_GEP70_COMPARATOR <- FALSE
GEP70_HIGH_LABEL <- "GEP70 > 0.66"

# Combine both >=3-copy 1q groups.
COMBINE_1Q3PLUS_GROUPS <- TRUE

# Alternative: combine CN2/PolyIG-low with CN>=3/PolyIG-high.
COMBINE_1Q2LOW_1Q3PLUSHIGH <- FALSE

if (COMBINE_1Q3PLUS_GROUPS && COMBINE_1Q2LOW_1Q3PLUSHIGH) {
  stop("Only one group-combination toggle can be TRUE.")
}

SHOW_LOGRANK_P_VALUE <- TRUE
LOGRANK_P_VALUE_COORD <- c(0.5, 0.85)
LOGRANK_P_VALUE_SIZE <- 7

SHOW_HR_LABEL <- TRUE
HR_LABEL_TEXT_SIZE <- 3.6
HR_REFERENCE_GROUP <- "CN2/PolyIG High"

LEGEND_NROW <- 1L
LEGEND_BYROW <- TRUE
LEGEND_TEXT_SIZE <- 13

BREAK_TIME_BY_YEARS <- 2

PLOT_WIDTH_IN <- 9
PLOT_HEIGHT_IN <- 8
PNG_DPI <- 300

################################################################################
# 2. HELPERS
################################################################################

fmt_p <- function(p) {
  if (!is.finite(p)) return("NA")
  
  if (p < 0.001) {
    format(p, scientific = TRUE, digits = 2)
  } else {
    formatC(p, format = "f", digits = 3)
  }
}

clean_group <- function(x) {
  x <- gsub("\n", " ", as.character(x), fixed = TRUE)
  trimws(gsub("[[:space:]]+", " ", x))
}

################################################################################
# 3. GROUP LABELS AND COLORS
################################################################################

LABEL_1Q2_LOW <- "CN2/PolyIG Low"
LABEL_1Q2_HIGH <- "CN2/PolyIG High"

LABEL_1Q3PLUS_LOW <- "CN>=3/PolyIG Low"
LABEL_1Q3PLUS_HIGH <- "CN>=3/PolyIG High"
LABEL_1Q3PLUS_ALL <- "CN>=3/PolyIG ALL"

LABEL_INTERMEDIATE <- "Intermediate"

if (COMBINE_1Q2LOW_1Q3PLUSHIGH) {
  BASE_GROUP_LEVELS <- c(
    LABEL_1Q2_HIGH,
    LABEL_INTERMEDIATE,
    LABEL_1Q3PLUS_LOW
  )
} else if (COMBINE_1Q3PLUS_GROUPS) {
  BASE_GROUP_LEVELS <- c(
    LABEL_1Q2_HIGH,
    LABEL_1Q2_LOW,
    LABEL_1Q3PLUS_ALL
  )
} else {
  BASE_GROUP_LEVELS <- c(
    LABEL_1Q2_HIGH,
    LABEL_1Q2_LOW,
    LABEL_1Q3PLUS_HIGH,
    LABEL_1Q3PLUS_LOW
  )
}

GROUP_LEVELS <- if (INCLUDE_GEP70_COMPARATOR) {
  c(
    BASE_GROUP_LEVELS,
    GEP70_HIGH_LABEL
  )
} else {
  BASE_GROUP_LEVELS
}

GROUP_COLORS <- c(
  "CN2/PolyIG Low" = "grey50",
  "CN2/PolyIG High" = "steelblue2",
  "CN>=3/PolyIG Low" = "red2",
  "CN>=3/PolyIG High" = "purple",
  "CN>=3/PolyIG ALL" = "red3",
  "Intermediate" = "grey50",
  "GEP70 > 0.66" = "black"
)

################################################################################
# 4. LOAD ESET
################################################################################

if (!file.exists(ESET_PATH)) {
  stop("ESET not found:\n", ESET_PATH)
}

eset <- readRDS(ESET_PATH)

if (!methods::is(eset, "ExpressionSet")) {
  stop("ESET_PATH did not contain an ExpressionSet.")
}

df <- as.data.frame(
  Biobase::pData(eset),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_columns <- c(
  SAMPLE_ID_COL,
  SAMPLE_GROUP_COL,
  TT_COL,
  GEP1Q_COL,
  POLYIG_COL,
  GEP70_COL,
  TIME_COL,
  EVENT_COL
)

missing_columns <- setdiff(
  required_columns,
  colnames(df)
)

if (length(missing_columns)) {
  stop(
    "Missing required pData columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

################################################################################
# 5. BUILD COHORT
################################################################################

df <- df |>
  transmute(
    CHIPID = .data[[SAMPLE_ID_COL]],
    sample_group = .data[[SAMPLE_GROUP_COL]],
    TT_On_for_Tx1 = .data[[TT_COL]],
    time = as.numeric(.data[[TIME_COL]]),
    event = as.numeric(.data[[EVENT_COL]]),
    GEP1q_class = .data[[GEP1Q_COL]],
    PolyIG = as.numeric(.data[[POLYIG_COL]]),
    GEP70 = as.numeric(.data[[GEP70_COL]])
  )

if (REQUIRE_BASELINE_RNAS_NDMM) {
  df <- df |>
    filter(
      sample_group == "RNAS_CD138_NDMM"
    )
}

if (FILTER_TT_PROTOCOLS) {
  df <- df |>
    filter(
      TT_On_for_Tx1 %in% TT_KEEP
    )
}

df <- df |>
  filter(
    !is.na(CHIPID),
    is.finite(time),
    time >= 0,
    event %in% c(0, 1),
    GEP1q_class %in% c("1q2", "1q3plus"),
    is.finite(PolyIG),
    is.finite(GEP70)
  ) |>
  mutate(
    GEP70_High = GEP70 > GEP70_THRESHOLD,
    GEP70_Low = GEP70 <= GEP70_THRESHOLD
  )

if (!INCLUDE_GEP70_COMPARATOR) {
  df <- df |>
    filter(
      GEP70_Low
    )
}

if (!nrow(df)) {
  stop("No samples remain after cohort filtering.")
}

################################################################################
# 6. OPTIONAL POLYIG CUTPOINT
################################################################################

if (USE_OPTIMAL_CUTPOINT) {
  
  d <- df |>
    filter(
      GEP70_Low
    )
  
  cuts <- sort(
    unique(
      d$PolyIG
    )
  )
  
  if (length(cuts) < 2) {
    stop("PolyIG has fewer than two unique values.")
  }
  
  cuts <- (
    cuts[-1] +
      cuts[-length(cuts)]
  ) / 2
  
  min_n <- ceiling(
    CUTPOINT_MINPROP *
      nrow(d)
  )
  
  scan <- bind_rows(
    lapply(
      cuts,
      function(cut) {
        
        high <- d$PolyIG > cut
        
        if (
          sum(!high) < min_n ||
          sum(high) < min_n
        ) {
          return(NULL)
        }
        
        fit_cut <- tryCatch(
          coxph(
            Surv(time, event) ~ factor(high),
            data = d
          ),
          error = function(e) NULL
        )
        
        if (is.null(fit_cut)) {
          return(NULL)
        }
        
        s <- summary(fit_cut)
        hr <- exp(coef(fit_cut)[1])
        ci <- exp(confint(fit_cut)[1, ])
        
        data.frame(
          cutpoint = cut,
          HR = hr,
          CI_lower_95 = ci[1],
          CI_upper_95 = ci[2],
          p_value = s$coef[1, "Pr(>|z|)"]
        )
      }
    )
  )
  
  if (!nrow(scan)) {
    stop("No valid PolyIG cutpoint remained.")
  }
  
  if (OPTIMAL_CUT_METHOD == "hr") {
    
    POLYIG_CUT <- scan |>
      mutate(
        metric = abs(log(HR))
      ) |>
      arrange(
        desc(metric),
        p_value
      ) |>
      slice(1) |>
      pull(cutpoint)
    
  } else {
    
    POLYIG_CUT <- scan |>
      arrange(
        p_value
      ) |>
      slice(1) |>
      pull(cutpoint)
  }
}

################################################################################
# 7. DEFINE GROUPS
################################################################################

df <- df |>
  mutate(
    PolyIG_group = if_else(
      PolyIG > POLYIG_CUT,
      "High",
      "Low"
    ),
    
    Group = case_when(
      INCLUDE_GEP70_COMPARATOR &
        GEP70_High ~
        GEP70_HIGH_LABEL,
      
      GEP70_Low &
        COMBINE_1Q2LOW_1Q3PLUSHIGH &
        GEP1q_class == "1q2" &
        PolyIG_group == "Low" ~
        LABEL_INTERMEDIATE,
      
      GEP70_Low &
        COMBINE_1Q2LOW_1Q3PLUSHIGH &
        GEP1q_class == "1q3plus" &
        PolyIG_group == "High" ~
        LABEL_INTERMEDIATE,
      
      GEP70_Low &
        COMBINE_1Q3PLUS_GROUPS &
        GEP1q_class == "1q3plus" ~
        LABEL_1Q3PLUS_ALL,
      
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
  stop("At least two groups are required.")
}

if (!HR_REFERENCE_GROUP %in% levels(df$Group)) {
  stop(
    "HR_REFERENCE_GROUP is absent: ",
    HR_REFERENCE_GROUP
  )
}

################################################################################
# 8. SURVIVAL MODEL
################################################################################

survival_fit <- survfit(
  Surv(time, event) ~ Group,
  data = df
)

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
    "A group is missing a defined color."
  )
}

logrank <- survdiff(
  Surv(time, event) ~ Group,
  data = df
)

overall_p <- 1 - pchisq(
  logrank$chisq,
  length(logrank$n) - 1
)

################################################################################
# 9. HRs VS REFERENCE
################################################################################

ordered_levels <- c(
  HR_REFERENCE_GROUP,
  setdiff(
    group_levels_present,
    HR_REFERENCE_GROUP
  )
)

reference_data <- df |>
  mutate(
    Group_hr = factor(
      as.character(Group),
      levels = ordered_levels
    )
  )

reference_fit <- coxph(
  Surv(time, event) ~ Group_hr,
  data = reference_data
)

reference_summary <- summary(
  reference_fit
)

reference_ci <- exp(
  confint(
    reference_fit
  )
)

reference_groups <- sub(
  "^Group_hr",
  "",
  rownames(
    reference_summary$coefficients
  )
)

reference_hr <- exp(
  coef(
    reference_fit
  )
)

reference_p <- reference_summary$coefficients[
  ,
  "Pr(>|z|)"
]

hr_label <- paste(
  c(
    paste0(
      "Reference: ",
      clean_group(
        HR_REFERENCE_GROUP
      )
    ),
    
    paste0(
      reference_groups,
      " vs ",
      clean_group(
        HR_REFERENCE_GROUP
      ),
      ": HR = ",
      sprintf(
        "%.2f (%.2f–%.2f)",
        reference_hr,
        reference_ci[, 1],
        reference_ci[, 2]
      ),
      "\np = ",
      vapply(
        reference_p,
        fmt_p,
        character(1)
      )
    )
  ),
  collapse = "\n"
)

################################################################################
# 10. ALL PAIRWISE HRs
################################################################################

pairwise_table <- bind_rows(
  lapply(
    combn(
      group_levels_present,
      2,
      simplify = FALSE
    ),
    
    function(pair) {
      
      reference_group <- pair[1]
      comparison_group <- pair[2]
      
      pair_data <- df |>
        filter(
          Group %in% pair
        ) |>
        mutate(
          Group_pair = factor(
            as.character(Group),
            levels = pair
          )
        )
      
      nr <- sum(
        pair_data$Group_pair ==
          reference_group
      )
      
      nc <- sum(
        pair_data$Group_pair ==
          comparison_group
      )
      
      er <- sum(
        pair_data$event[
          pair_data$Group_pair ==
            reference_group
        ] == 1
      )
      
      ec <- sum(
        pair_data$event[
          pair_data$Group_pair ==
            comparison_group
        ] == 1
      )
      
      empty <- data.frame(
        reference = reference_group,
        comparison = comparison_group,
        n_reference = nr,
        n_comparison = nc,
        events_reference = er,
        events_comparison = ec,
        HR = NA_real_,
        CI_lower_95 = NA_real_,
        CI_upper_95 = NA_real_,
        p_value = NA_real_
      )
      
      if (
        nr < 2 ||
        nc < 2 ||
        er + ec == 0
      ) {
        return(empty)
      }
      
      fit_pair <- tryCatch(
        coxph(
          Surv(time, event) ~ Group_pair,
          data = pair_data
        ),
        error = function(e) NULL
      )
      
      if (is.null(fit_pair)) {
        return(empty)
      }
      
      s <- summary(
        fit_pair
      )
      
      ci <- exp(
        confint(
          fit_pair
        )
      )
      
      data.frame(
        reference = reference_group,
        comparison = comparison_group,
        n_reference = nr,
        n_comparison = nc,
        events_reference = er,
        events_comparison = ec,
        HR = exp(coef(fit_pair)[1]),
        CI_lower_95 = ci[1, 1],
        CI_upper_95 = ci[1, 2],
        p_value = s$coef[1, "Pr(>|z|)"]
      )
    }
  )
)

################################################################################
# 11. FIGURE LABELS
################################################################################

x_max <- ceiling(
  max(
    df$time,
    na.rm = TRUE
  ) /
    BREAK_TIME_BY_YEARS
) * BREAK_TIME_BY_YEARS

plot_title <-
  "Time to Progression by GEP 1q CN × PolyIG (PI/IMiD)"

protocol_subtitle <- if (FILTER_TT_PROTOCOLS) {
  "Standard Risk (GEP70); TT3a, TT3b, TT4-S"
} else {
  "Standard Risk (GEP70); all treatment protocols"
}

plot_subtitle <- paste0(
  protocol_subtitle,
  " (n = ",
  nrow(df),
  ")"
)

y_axis_label <-
  "Cumulative Progression Probability"

################################################################################
# 12. KAPLAN–MEIER PLOT
################################################################################

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
  xlim = c(0, x_max),
  ylim = c(0, 1),
  
  title = plot_title,
  xlab = "Time from First Transplant (years)",
  ylab = y_axis_label,
  
  ggtheme = ggplot2::theme_bw(
    base_size = 16
  ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold"
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5,
        size = 18
      )
    ),
  
  tables.theme =
    survminer::theme_cleantable()
)

km_plot$plot <- km_plot$plot +
  ggplot2::labs(
    subtitle = plot_subtitle
  ) +
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
    labels = function(x) sprintf("%.1f", x)
  )

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
      y = y_range[2] -
        0.04 * (
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
  heights = c(2, 0.6),
  align = "v"
)

print(
  combined_plot
)

################################################################################
# 13. SAVE FIGURE 3B
################################################################################

png_file <- file.path(
  FIGURE_DIR,
  paste0(
    FIGURE_ID,
    "_GEP1q_lasso_x_PolyIG_",
    ENDPOINT,
    "_KM.png"
  )
)

pdf_file <- file.path(
  FIGURE_DIR,
  paste0(
    FIGURE_ID,
    "_GEP1q_lasso_x_PolyIG_",
    ENDPOINT,
    "_KM.pdf"
  )
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
# 14. PAIRWISE HR / P-VALUE CSV
################################################################################

pairwise_file <- file.path(
  FIGURE_DIR,
  paste0(
    FIGURE_ID,
    "_pairwise-HR-pvalues.csv"
  )
)

pairwise_export <- pairwise_table |>
  mutate(
    Reference = clean_group(
      reference
    ),
    
    Comparison_Group = clean_group(
      comparison
    ),
    
    Comparison = paste0(
      Comparison_Group,
      " vs ",
      Reference
    ),
    
    `HR (95% CI)` = ifelse(
      is.finite(HR) &
        is.finite(CI_lower_95) &
        is.finite(CI_upper_95),
      
      sprintf(
        "%.2f (%.2f–%.2f)",
        HR,
        CI_lower_95,
        CI_upper_95
      ),
      
      "NA"
    ),
    
    `p-value` = vapply(
      p_value,
      fmt_p,
      character(1)
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

readr::write_csv(
  pairwise_export,
  pairwise_file
)

################################################################################
# 15. SAVE ANALYSIS COHORT
################################################################################

cohort_file <- file.path(
  FIGURE_DIR,
  paste0(
    FIGURE_ID,
    "_analysis-cohort.csv"
  )
)

readr::write_csv(
  df,
  cohort_file
)

################################################################################
# 16. FINAL SUMMARY
################################################################################

message(
  "\nStep 26B complete.",
  "\nEndpoint: TTP",
  "\nFinal cohort: ", nrow(df),
  "\nEvents: ", sum(df$event == 1),
  "\nPolyIG cutoff: ", formatC(POLYIG_CUT, format = "f", digits = 4),
  "\nOverall log-rank p: ", fmt_p(overall_p),
  "\nActive groups: ", paste(group_levels_present, collapse = "; "),
  "\nPNG: ", png_file,
  "\nPDF: ", pdf_file,
  "\nPairwise HR/p-value CSV: ", pairwise_file,
  "\nAnalysis cohort CSV: ", cohort_file,
  "\n"
)