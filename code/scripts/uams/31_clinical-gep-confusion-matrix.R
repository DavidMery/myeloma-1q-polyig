################################################################################
# STEP 31 — FIGURE 5D
# Clinical FISH1q × unIGH vs GEP1q × PolyIG confusion matrix
#
# Rows    = Clinical FISH1q + unIGH
# Columns = GEP1q + PolyIG
#
# Cell labels:
#   Patient count
#   Percentage within clinical row
#
# Black outline = exact agreement
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
})

################################################################################
# 0. VERIFY PIPELINE ROOTS
################################################################################

if (!exists("RESULTS_ROOT")) {
  stop(
    "RESULTS_ROOT is not defined. ",
    "Run this script through code/run_pipeline.R."
  )
}

dir.create(
  RESULTS_ROOT,
  recursive = TRUE,
  showWarnings = FALSE
)

################################################################################
# 1. FILES / SETTINGS
################################################################################

MATCHED_FILE <- file.path(
  RESULTS_ROOT, "figures", "Fig5", "Fig5C",
  "Fig5C_matched-analysis-cohort.csv"
)

OUT_DIR <- file.path(RESULTS_ROOT, "figures", "Fig5", "Fig5D")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

RISK_LEVELS <- c("Favorable", "Intermediate", "Adverse")

PLOT_WIDTH_IN <- 8.5
PLOT_HEIGHT_IN <- 7.5
PNG_DPI <- 300

################################################################################
# 2. READ MATCHED CLASSIFICATIONS FROM STEP 30
################################################################################

if (!file.exists(MATCHED_FILE)) {
  stop(
    "Step 30 matched cohort file not found:\n",
    MATCHED_FILE
  )
}

matched <- read_csv(MATCHED_FILE, show_col_types = FALSE) |>
  transmute(
    CHIPID = as.character(CHIPID),
    Clinical_Risk = factor(as.character(Clinical_Risk),
                           levels = RISK_LEVELS, ordered = TRUE),
    GEP_Risk = factor(as.character(GEP_Risk),
                      levels = RISK_LEVELS, ordered = TRUE)
  ) |>
  filter(!is.na(Clinical_Risk), !is.na(GEP_Risk))

if (anyDuplicated(matched$CHIPID)) {
  stop("Duplicate CHIPIDs are present in the Step 30 matched cohort.")
}

################################################################################
# 3. CONFUSION MATRIX
################################################################################

confusion_matrix <- table(
  Clinical = matched$Clinical_Risk,
  GEP = matched$GEP_Risk
)

n_total <- sum(confusion_matrix)
exact_agreement_n <- sum(diag(confusion_matrix))
exact_agreement_rate <- exact_agreement_n / n_total

################################################################################
# 4. COHEN'S KAPPA
################################################################################

row_prop <- rowSums(confusion_matrix) / n_total
col_prop <- colSums(confusion_matrix) / n_total
expected_agreement <- sum(row_prop * col_prop)

cohen_kappa <- (
  exact_agreement_rate - expected_agreement
) / (
  1 - expected_agreement
)

################################################################################
# 5. LONG-FORM MATRIX DATA
################################################################################

plot_data <- as.data.frame(confusion_matrix) |>
  as_tibble() |>
  rename(
    Clinical_Risk = Clinical,
    GEP_Risk = GEP,
    n = Freq
  ) |>
  mutate(
    Clinical_Risk = factor(Clinical_Risk,
                           levels = RISK_LEVELS,
                           ordered = TRUE),
    GEP_Risk = factor(GEP_Risk,
                      levels = RISK_LEVELS,
                      ordered = TRUE)
  ) |>
  group_by(Clinical_Risk) |>
  mutate(
    clinical_total = sum(n),
    row_percent = 100 * n / clinical_total
  ) |>
  ungroup() |>
  mutate(
    exact_agreement = as.character(Clinical_Risk) == as.character(GEP_Risk),
    cell_label = paste0(n, "\n", sprintf("%.1f%%", row_percent)),
    
    # Reverse clinical order so Favorable appears at top
    Clinical_Risk_Plot = factor(
      as.character(Clinical_Risk),
      levels = rev(RISK_LEVELS)
    ),
    
    GEP_Risk_Plot = factor(
      as.character(GEP_Risk),
      levels = RISK_LEVELS
    )
  )

################################################################################
# 6. SUBTITLE
################################################################################

plot_subtitle <- paste0(
  "Matched patients: n = ", n_total,
  " | Exact agreement = ", sprintf("%.1f%%", 100 * exact_agreement_rate),
  " | Cohen's \u03BA = ", sprintf("%.3f", cohen_kappa)
)

################################################################################
# 7. CONFUSION-MATRIX PLOT
################################################################################

confusion_plot <- ggplot(
  plot_data,
  aes(
    x = GEP_Risk_Plot,
    y = Clinical_Risk_Plot,
    fill = row_percent
  )
) +
  
  # Main heatmap
  geom_tile(
    color = "white",
    linewidth = 1.2
  ) +
  
  # Black border around exact-agreement cells
  geom_tile(
    data = filter(plot_data, exact_agreement),
    fill = NA,
    color = "black",
    linewidth = 1.4
  ) +
  
  # Count + within-row percentage
  geom_text(
    aes(label = cell_label),
    size = 5.2,
    fontface = "bold",
    lineheight = 1.05
  ) +
  
  scale_fill_gradient(
    low = "white",
    high = "steelblue4",
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    labels = function(x) paste0(x, "%"),
    name = "Within-clinical\nrow percentage"
  ) +
  
  coord_fixed() +
  
  labs(
    title = "Clinical and GEP-Derived Risk-Group Agreement",
    subtitle = plot_subtitle,
    x = "GEP1q + PolyIG Risk Group",
    y = "Clinical FISH1q + unIGH Risk Group",
    caption = paste0(
      "Cell labels show patient count and percentage within each clinical row. ",
      "Black outlines mark exact agreement."
    )
  ) +
  
  theme_bw(base_size = 15) +
  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 19),
    plot.subtitle = element_text(hjust = 0.5, size = 12.5),
    plot.caption = element_text(hjust = 0, size = 10),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 13),
    panel.grid = element_blank(),
    legend.title = element_text(face = "bold"),
    plot.margin = margin(12, 18, 12, 18)
  )

print(confusion_plot)

################################################################################
# 8. SAVE FIGURE
################################################################################

PNG_FILE <- file.path(
  OUT_DIR,
  "Fig5D_Clinical_vs_GEP_ConfusionMatrix.png"
)

PDF_FILE <- file.path(
  OUT_DIR,
  "Fig5D_Clinical_vs_GEP_ConfusionMatrix.pdf"
)

ggsave(
  PNG_FILE,
  confusion_plot,
  width = PLOT_WIDTH_IN,
  height = PLOT_HEIGHT_IN,
  units = "in",
  dpi = PNG_DPI,
  bg = "white"
)

ggsave(
  PDF_FILE,
  confusion_plot,
  width = PLOT_WIDTH_IN,
  height = PLOT_HEIGHT_IN,
  units = "in",
  device = "pdf",
  bg = "white"
)

################################################################################
# 9. SAVE MATRIX DATA
################################################################################

COUNTS_FILE <- file.path(
  OUT_DIR,
  "Fig5D_Clinical_vs_GEP_ConfusionMatrix_data.csv"
)

write_csv(
  plot_data |>
    transmute(
      Clinical_Risk = as.character(Clinical_Risk),
      GEP_Risk = as.character(GEP_Risk),
      n,
      clinical_total,
      row_percent,
      exact_agreement
    ),
  COUNTS_FILE
)

################################################################################
# 10. FINAL
################################################################################

cat(
  "\nStep 31 — Figure 5D complete\n",
  "Matched n = ", n_total, "\n",
  "Exact agreement = ", exact_agreement_n, "/", n_total,
  " (", sprintf("%.1f%%", 100 * exact_agreement_rate), ")\n",
  "Cohen's kappa = ", sprintf("%.3f", cohen_kappa), "\n",
  "PNG: ", PNG_FILE, "\n",
  "PDF: ", PDF_FILE, "\n",
  "Data: ", COUNTS_FILE, "\n",
  sep = ""
)

cat("\nConfusion matrix:\n")
print(confusion_matrix)