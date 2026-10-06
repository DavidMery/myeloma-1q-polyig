################################################################################
# STEP 25B — FIGURE 3A
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
# 1. SETTINGS
################################################################################

FIGURE_ID <- "Fig3A_Lasso"
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
CUTPOINT_MINPROP <- 0.15

TT_KEEP <- c(
  "TT2_NoThal",
  "TT2_Thal"
)

FILTER_TT_PROTOCOLS <- TRUE
REQUIRE_BASELINE_RNAS_NDMM <- TRUE

GEP70_THRESHOLD <- 0.66
INCLUDE_GEP70_COMPARATOR <- FALSE

# Combine CN>=3/PolyIG-low and CN>=3/PolyIG-high.
COMBINE_1Q3PLUS_GROUPS <- TRUE

# Alternative grouping option.
COMBINE_1Q2LOW_1Q3PLUSHIGH <- FALSE

if (COMBINE_1Q3PLUS_GROUPS && COMBINE_1Q2LOW_1Q3PLUSHIGH) {
  stop("Only one group-combination toggle can be TRUE.")
}

BREAK_TIME_BY_YEARS <- 2

SHOW_LOGRANK_P_VALUE <- TRUE
LOGRANK_P_VALUE_COORD <- c(0.5, 0.85)
LOGRANK_P_VALUE_SIZE <- 7

SHOW_HR_LABEL <- TRUE
HR_LABEL_TEXT_SIZE <- 3.6
HR_REFERENCE_GROUP <- "CN2/PolyIG High"

PLOT_WIDTH_IN <- 9
PLOT_HEIGHT_IN <- 8
PNG_DPI <- 300

################################################################################
# 2. HELPERS
################################################################################

p_format <- function(x) {
  if (!is.finite(x)) return("NA")
  if (x < 0.001) {
    format(x, scientific = TRUE, digits = 2)
  } else {
    formatC(x, format = "f", digits = 3)
  }
}

################################################################################
# 3. LOAD DATA
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

missing_columns <- setdiff(required_columns, colnames(df))

if (length(missing_columns)) {
  stop(
    "Missing required pData columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

################################################################################
# 4. BUILD ANALYSIS COHORT
################################################################################

df <- df |>
  transmute(
    CHIPID = .data[[SAMPLE_ID_COL]],
    sample_group = .data[[SAMPLE_GROUP_COL]],
    TT_On_for_Tx1 = .data[[TT_COL]],
    time = as.numeric(.data[[TIME_COL]]),
    event = as.numeric(.data[[EVENT_COL]]),
    GEP1q = .data[[GEP1Q_COL]],
    PolyIG = as.numeric(.data[[POLYIG_COL]]),
    GEP70 = as.numeric(.data[[GEP70_COL]])
  )

if (REQUIRE_BASELINE_RNAS_NDMM) {
  df <- df |>
    filter(sample_group == "RNAS_CD138_NDMM")
}

if (FILTER_TT_PROTOCOLS) {
  df <- df |>
    filter(TT_On_for_Tx1 %in% TT_KEEP)
}

df <- df |>
  filter(
    !is.na(CHIPID),
    is.finite(time),
    time >= 0,
    event %in% c(0, 1),
    GEP1q %in% c("1q2", "1q3plus"),
    is.finite(PolyIG),
    is.finite(GEP70)
  )

################################################################################
# 5. OPTIONAL POLYIG CUTPOINT SEARCH
################################################################################

if (USE_OPTIMAL_CUTPOINT) {
  
  d <- df |>
    filter(GEP70 <= GEP70_THRESHOLD)
  
  cuts <- sort(unique(d$PolyIG))
  cuts <- (cuts[-1] + cuts[-length(cuts)]) / 2
  
  min_n <- ceiling(
    CUTPOINT_MINPROP * nrow(d)
  )
  
  scan <- lapply(cuts, function(cut) {
    
    high <- d$PolyIG > cut
    
    if (sum(high) < min_n || sum(!high) < min_n) {
      return(NULL)
    }
    
    fit_cut <- tryCatch(
      coxph(
        Surv(time, event) ~ high,
        data = d
      ),
      error = function(e) NULL
    )
    
    if (is.null(fit_cut)) return(NULL)
    
    s <- summary(fit_cut)
    
    data.frame(
      cutpoint = cut,
      HR = exp(coef(fit_cut)[1]),
      p = s$coef[1, "Pr(>|z|)"]
    )
  })
  
  scan <- bind_rows(scan)
  
  if (!nrow(scan)) {
    stop("No valid PolyIG cutpoint was found.")
  }
  
  POLYIG_CUT <- scan |>
    arrange(p, desc(abs(log(HR)))) |>
    slice(1) |>
    pull(cutpoint)
}

################################################################################
# 6. DEFINE GROUPS
################################################################################

CN2_LOW <- "CN2/PolyIG Low"
CN2_HIGH <- "CN2/PolyIG High"
CN3_LOW <- "CN>=3/PolyIG Low"
CN3_HIGH <- "CN>=3/PolyIG High"
CN3_ALL <- "CN>=3/PolyIG ALL"
INTERMEDIATE <- "Intermediate"
GEP70_HIGH <- "GEP70 > 0.66"

if (COMBINE_1Q2LOW_1Q3PLUSHIGH) {
  group_levels <- c(
    CN2_HIGH,
    INTERMEDIATE,
    CN3_LOW
  )
} else if (COMBINE_1Q3PLUS_GROUPS) {
  group_levels <- c(
    CN2_HIGH,
    CN2_LOW,
    CN3_ALL
  )
} else {
  group_levels <- c(
    CN2_HIGH,
    CN2_LOW,
    CN3_HIGH,
    CN3_LOW
  )
}

if (INCLUDE_GEP70_COMPARATOR) {
  group_levels <- c(
    group_levels,
    GEP70_HIGH
  )
}

group_colors <- c(
  "CN2/PolyIG Low" = "grey50",
  "CN2/PolyIG High" = "steelblue2",
  "CN>=3/PolyIG Low" = "red2",
  "CN>=3/PolyIG High" = "purple",
  "CN>=3/PolyIG ALL" = "red3",
  "Intermediate" = "grey50",
  "GEP70 > 0.66" = "black"
)

df <- df |>
  mutate(
    PolyIG_group = if_else(
      PolyIG > POLYIG_CUT,
      "High",
      "Low"
    ),
    
    Group = case_when(
      INCLUDE_GEP70_COMPARATOR &
        GEP70 > GEP70_THRESHOLD ~
        GEP70_HIGH,
      
      GEP70 <= GEP70_THRESHOLD &
        COMBINE_1Q2LOW_1Q3PLUSHIGH &
        GEP1q == "1q2" &
        PolyIG_group == "Low" ~
        INTERMEDIATE,
      
      GEP70 <= GEP70_THRESHOLD &
        COMBINE_1Q2LOW_1Q3PLUSHIGH &
        GEP1q == "1q3plus" &
        PolyIG_group == "High" ~
        INTERMEDIATE,
      
      GEP70 <= GEP70_THRESHOLD &
        COMBINE_1Q3PLUS_GROUPS &
        GEP1q == "1q3plus" ~
        CN3_ALL,
      
      GEP70 <= GEP70_THRESHOLD &
        GEP1q == "1q2" &
        PolyIG_group == "Low" ~
        CN2_LOW,
      
      GEP70 <= GEP70_THRESHOLD &
        GEP1q == "1q2" &
        PolyIG_group == "High" ~
        CN2_HIGH,
      
      GEP70 <= GEP70_THRESHOLD &
        GEP1q == "1q3plus" &
        PolyIG_group == "Low" ~
        CN3_LOW,
      
      GEP70 <= GEP70_THRESHOLD &
        GEP1q == "1q3plus" &
        PolyIG_group == "High" ~
        CN3_HIGH,
      
      TRUE ~
        NA_character_
    ),
    
    Group = factor(
      Group,
      levels = group_levels
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

################################################################################
# 7. SURVIVAL ANALYSIS
################################################################################

fit <- survfit(
  Surv(time, event) ~ Group,
  data = df
)

lr <- survdiff(
  Surv(time, event) ~ Group,
  data = df
)

overall_p <- 1 - pchisq(
  lr$chisq,
  length(lr$n) - 1
)

################################################################################
# 8. ALL PAIRWISE HRs
################################################################################

pairwise <- bind_rows(
  lapply(
    combn(
      levels(df$Group),
      2,
      simplify = FALSE
    ),
    function(pair) {
      
      d <- df |>
        filter(
          Group %in% pair
        ) |>
        mutate(
          PairGroup = factor(
            as.character(Group),
            levels = pair
          )
        )
      
      m <- coxph(
        Surv(time, event) ~ PairGroup,
        data = d
      )
      
      s <- summary(m)
      ci <- exp(confint(m))
      
      data.frame(
        Reference = pair[1],
        Comparison = pair[2],
        HR = exp(coef(m)[1]),
        CI_lower_95 = ci[1, 1],
        CI_upper_95 = ci[1, 2],
        p_value = s$coef[1, "Pr(>|z|)"],
        stringsAsFactors = FALSE
      )
    }
  )
)

################################################################################
# 9. HRs VS REFERENCE
################################################################################

if (!HR_REFERENCE_GROUP %in% levels(df$Group)) {
  stop(
    "HR reference group is absent: ",
    HR_REFERENCE_GROUP
  )
}

ref_levels <- c(
  HR_REFERENCE_GROUP,
  setdiff(
    levels(df$Group),
    HR_REFERENCE_GROUP
  )
)

ref_fit <- coxph(
  Surv(time, event) ~ factor(Group, levels = ref_levels),
  data = df
)

ref_s <- summary(ref_fit)
ref_ci <- exp(confint(ref_fit))

ref_names <- sub(
  "^factor\\(Group, levels = ref_levels\\)",
  "",
  rownames(ref_s$coefficients)
)

hr_label <- paste(
  c(
    paste0(
      "Reference: ",
      HR_REFERENCE_GROUP
    ),
    paste0(
      ref_names,
      " vs ",
      HR_REFERENCE_GROUP,
      ": HR = ",
      sprintf(
        "%.2f (%.2f–%.2f)",
        exp(coef(ref_fit)),
        ref_ci[, 1],
        ref_ci[, 2]
      ),
      "\np = ",
      vapply(
        ref_s$coefficients[, "Pr(>|z|)"],
        p_format,
        character(1)
      )
    )
  ),
  collapse = "\n"
)

################################################################################
# 10. PLOT
################################################################################

palette_used <- unname(
  group_colors[
    levels(df$Group)
  ]
)

if (anyNA(palette_used)) {
  stop("Missing group color.")
}

x_max <- ceiling(
  max(df$time, na.rm = TRUE) /
    BREAK_TIME_BY_YEARS
) * BREAK_TIME_BY_YEARS

plot_title <-
  "Time to Progression by GEP 1q CN × PolyIG (No PI)"

plot_subtitle <- paste0(
  "Standard Risk (GEP70); ",
  if (FILTER_TT_PROTOCOLS) {
    "TT2"
  } else {
    "all treatment protocols"
  },
  " (n = ",
  nrow(df),
  ")"
)

km <- ggsurvplot(
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
  legend.labs = levels(df$Group),
  
  break.time.by = BREAK_TIME_BY_YEARS,
  xlim = c(0, x_max),
  ylim = c(0, 1),
  
  title = plot_title,
  xlab = "Time from First Transplant (years)",
  ylab = "Cumulative Progression Probability",
  
  ggtheme = theme_bw(base_size = 16) +
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

km$plot <- km$plot +
  labs(
    subtitle = plot_subtitle
  ) +
  guides(
    color = guide_legend(
      nrow = 1,
      byrow = TRUE
    )
  ) +
  theme(
    legend.text = element_text(
      size = 13
    )
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(
      0,
      1,
      by = 0.2
    ),
    labels = function(x) sprintf("%.1f", x)
  )

if (SHOW_HR_LABEL) {
  
  b <- ggplot_build(
    km$plot
  )
  
  xr <- b$layout$panel_params[[1]]$x.range
  yr <- b$layout$panel_params[[1]]$y.range
  
  km$plot <- km$plot +
    annotate(
      "text",
      x = xr[2],
      y = yr[2] - 0.04 * (yr[2] - yr[1]),
      label = hr_label,
      hjust = 1.08,
      vjust = 1,
      size = HR_LABEL_TEXT_SIZE,
      lineheight = 0.95
    )
}

km$table <- km$table +
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
  km$plot,
  km$table,
  ncol = 1,
  heights = c(2, 0.6),
  align = "v"
)

print(
  combined_plot
)

################################################################################
# 11. OUTPUTS
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

pairwise_csv <- file.path(
  FIGURE_DIR,
  paste0(
    FIGURE_ID,
    "_pairwise-HR-pvalues.csv"
  )
)

cohort_csv <- file.path(
  FIGURE_DIR,
  paste0(
    FIGURE_ID,
    "_analysis-cohort.csv"
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

write_csv(
  pairwise,
  pairwise_csv
)

write_csv(
  df,
  cohort_csv
)

message(
  "\nStep 25B complete.",
  "\nEndpoint = TTP",
  "\nN = ", nrow(df),
  "\nEvents = ", sum(df$event == 1),
  "\nPolyIG cutoff = ", POLYIG_CUT,
  "\nOverall log-rank p = ", p_format(overall_p),
  "\nPNG: ", png_file,
  "\nPDF: ", pdf_file,
  "\nPairwise HR/p-value CSV: ", pairwise_csv,
  "\nAnalysis cohort CSV: ", cohort_csv,
  "\n"
)