################################################################################
# STEP 30 — FIGURE 5C
# Clinical FISH1q × unIGH vs GEP1q × PolyIG
# Overlay in the same matched patients
################################################################################

suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(grid)
  library(gridExtra)
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

GEP_FILE <- file.path(
  RESULTS_ROOT, "figures", "Fig5", "Fig5A",
  "Fig5A_GEP1q_lasso_x_PolyIG_TTP_analysis-cohort.csv"
)

CLINICAL_FILE <- file.path(
  RESULTS_ROOT, "figures", "Fig5", "Fig5B",
  "Fig5B_FISH1q_x_unIGH_ARD_TTP_analysis-cohort.csv"
)

OUT_DIR <- file.path(RESULTS_ROOT, "figures", "Fig5", "Fig5C")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

RISK_LEVELS <- c("Favorable", "Intermediate", "Adverse")
RISK_COLORS <- c(
  "Favorable" = "steelblue2",
  "Intermediate" = "grey50",
  "Adverse" = "red3"
)

BREAK_TIME_BY <- 2

################################################################################
# 2. READ STEP 28 AND STEP 29 COHORTS
################################################################################

if (!file.exists(GEP_FILE)) {
  stop(
    "Step 28 analysis cohort not found: ",
    GEP_FILE
  )
}

if (!file.exists(CLINICAL_FILE)) {
  stop(
    "Step 29 analysis cohort not found: ",
    CLINICAL_FILE
  )
}

gep <- read_csv(GEP_FILE, show_col_types = FALSE) |>
  transmute(
    CHIPID = as.character(CHIPID),
    time_gep = as.numeric(time),
    event_gep = as.numeric(event),
    GEP_Group = as.character(Group)
  )

clinical <- read_csv(CLINICAL_FILE, show_col_types = FALSE) |>
  transmute(
    CHIPID = as.character(CHIPID),
    time_clinical = as.numeric(time),
    event_clinical = as.numeric(event),
    Clinical_Group = as.character(Group)
  )

################################################################################
# 3. MATCH PATIENTS AND DEFINE COMMON RISK GROUPS
################################################################################

matched <- inner_join(gep, clinical, by = "CHIPID")

stopifnot(
  all(abs(matched$time_gep - matched$time_clinical) < 1e-8),
  all(matched$event_gep == matched$event_clinical)
)

matched <- matched |>
  mutate(
    time = time_gep,
    event = event_gep,
    
    Clinical_Risk = case_when(
      Clinical_Group == "CN2/unIGH High" ~ "Favorable",
      Clinical_Group == "CN2/unIGH Low" ~ "Intermediate",
      Clinical_Group == "CN≥3/unIGH ALL" ~ "Adverse",
      TRUE ~ NA_character_
    ),
    
    GEP_Risk = case_when(
      GEP_Group == "CN2/PolyIG High" ~ "Favorable",
      GEP_Group == "CN2/PolyIG Low" ~ "Intermediate",
      GEP_Group == "CN≥3/PolyIG ALL" ~ "Adverse",
      TRUE ~ NA_character_
    ),
    
    Clinical_Risk = factor(Clinical_Risk, levels = RISK_LEVELS),
    GEP_Risk = factor(GEP_Risk, levels = RISK_LEVELS)
  ) |>
  filter(!is.na(Clinical_Risk), !is.na(GEP_Risk))

N_MATCHED <- nrow(matched)

################################################################################
# 4. LOG-RANK P VALUES
################################################################################

get_logrank_p <- function(group) {
  fit <- survdiff(Surv(time, event) ~ group, data = matched)
  pchisq(fit$chisq, df = length(fit$n) - 1, lower.tail = FALSE)
}

fmt_p <- function(p) format(p, scientific = TRUE, digits = 2)

clinical_p <- get_logrank_p(matched$Clinical_Risk)
gep_p <- get_logrank_p(matched$GEP_Risk)

################################################################################
# 5. FIT CLINICAL AND GEP CURVES SEPARATELY
################################################################################

clinical_fit <- survfit(Surv(time, event) ~ Clinical_Risk, data = matched)
gep_fit <- survfit(Surv(time, event) ~ GEP_Risk, data = matched)

extract_curve <- function(fit, prefix, measurement) {
  s <- summary(fit, censored = TRUE)
  
  tibble(
    time = s$time,
    surv = s$surv,
    n.censor = s$n.censor,
    Risk_Group = sub(paste0("^", prefix, "="), "", as.character(s$strata)),
    Measurement = measurement,
    CumProgression = 1 - s$surv
  ) |>
    mutate(Risk_Group = factor(Risk_Group, levels = RISK_LEVELS))
}

clinical_curve <- extract_curve(clinical_fit, "Clinical_Risk", "Clinical")
gep_curve <- extract_curve(gep_fit, "GEP_Risk", "GEP")

zero_rows <- function(measurement) {
  tibble(
    time = 0,
    surv = 1,
    n.censor = 0,
    Risk_Group = factor(RISK_LEVELS, levels = RISK_LEVELS),
    Measurement = measurement,
    CumProgression = 0
  )
}

clinical_curve <- bind_rows(zero_rows("Clinical"), clinical_curve)
gep_curve <- bind_rows(zero_rows("GEP"), gep_curve)

################################################################################
# 6. AXIS SETTINGS
################################################################################

x_max <- ceiling(max(matched$time, na.rm = TRUE) / BREAK_TIME_BY) * BREAK_TIME_BY
risk_times <- seq(0, x_max, by = BREAK_TIME_BY)

################################################################################
# 7. MAIN OVERLAY PLOT
################################################################################

main_plot <- ggplot() +
  
  # Clinical = solid
  geom_step(
    data = clinical_curve,
    aes(x = time, y = CumProgression, color = Risk_Group, group = Risk_Group),
    linewidth = 1.05, linetype = "solid", direction = "hv"
  ) +
  
  # GEP = dashed
  geom_step(
    data = gep_curve,
    aes(x = time, y = CumProgression, color = Risk_Group, group = Risk_Group),
    linewidth = 1.05, linetype = "22", direction = "hv"
  ) +
  
  # Clinical censor marks
  geom_point(
    data = filter(clinical_curve, n.censor > 0),
    aes(x = time, y = CumProgression, color = Risk_Group),
    shape = 3, size = 1.8, stroke = 0.55
  ) +
  
  # GEP censor marks
  geom_point(
    data = filter(gep_curve, n.censor > 0),
    aes(x = time, y = CumProgression, color = Risk_Group),
    shape = 3, size = 1.8, stroke = 0.55
  ) +
  
  scale_color_manual(
    values = RISK_COLORS,
    breaks = RISK_LEVELS,
    name = "Risk Group"
  ) +
  
  scale_x_continuous(breaks = risk_times) +
  
  scale_y_continuous(
    breaks = seq(0, 1, 0.2),
    labels = function(x) sprintf("%.1f", x)
  ) +
  
  coord_cartesian(
    xlim = c(-0.20, x_max + 0.20),
    ylim = c(-0.015, 1.015),
    clip = "off"
  ) +
  
  labs(
    title = "Clinical and GEP-Derived Risk Groups",
    subtitle = paste0(
      "Same matched patients (n = ", N_MATCHED, ")\n",
      "Clinical log-rank p = ", fmt_p(clinical_p),
      "; GEP log-rank p = ", fmt_p(gep_p)
    ),
    x = "Time from First Transplant (years)",
    y = "Cumulative Progression Probability"
  ) +
  
  theme_bw(base_size = 16) +
  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 22),
    plot.subtitle = element_text(hjust = 0.5, size = 17),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.box.just = "center",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 12),
    legend.spacing.y = unit(4, "pt"),
    plot.margin = margin(15, 20, 10, 20)
  )

################################################################################
# 8. SECOND LEGEND ROW FOR MEASUREMENT
################################################################################

legend_df <- tibble(
  x = c(NA_real_, NA_real_),
  y = c(NA_real_, NA_real_),
  Measurement = factor(
    c("Clinical FISH 1q + unIGH", "GEP1q + PolyIG"),
    levels = c("Clinical FISH 1q + unIGH", "GEP1q + PolyIG")
  )
)

main_plot <- main_plot +
  geom_line(
    data = legend_df,
    aes(x = x, y = y, linetype = Measurement),
    color = "black", linewidth = 1.2, na.rm = TRUE
  ) +
  scale_linetype_manual(
    values = c(
      "Clinical FISH 1q + unIGH" = "solid",
      "GEP1q + PolyIG" = "22"
    ),
    name = "Measurement"
  ) +
  guides(
    color = guide_legend(
      order = 1, nrow = 1, byrow = TRUE,
      title.position = "left", title.hjust = 0,
      override.aes = list(linetype = "solid", linewidth = 1.2)
    ),
    linetype = guide_legend(
      order = 2, nrow = 1, byrow = TRUE,
      title.position = "left", title.hjust = 0,
      override.aes = list(color = "black", linewidth = 1.2)
    )
  )

################################################################################
# 9. RISK TABLE
################################################################################

risk_counts <- function(group, measurement) {
  bind_rows(lapply(RISK_LEVELS, function(g) {
    d <- matched[group == g, ]
    
    s <- summary(
      survfit(Surv(time, event) ~ 1, data = d),
      times = risk_times,
      extend = TRUE
    )
    
    tibble(
      time = risk_times,
      n = s$n.risk,
      Risk_Group = g,
      Measurement = measurement
    )
  }))
}

risk_df <- bind_rows(
  risk_counts(matched$Clinical_Risk, "Clinical"),
  risk_counts(matched$GEP_Risk, "GEP")
) |>
  mutate(
    row_label = paste(Risk_Group, Measurement, sep = " — "),
    y = case_when(
      row_label == "Favorable — Clinical" ~ 6,
      row_label == "Favorable — GEP" ~ 5,
      row_label == "Intermediate — Clinical" ~ 4,
      row_label == "Intermediate — GEP" ~ 3,
      row_label == "Adverse — Clinical" ~ 2,
      row_label == "Adverse — GEP" ~ 1
    ),
    Risk_Group = factor(Risk_Group, levels = RISK_LEVELS)
  )

################################################################################
# 10. RISK TABLE LABELS
################################################################################

label_df <- distinct(risk_df, row_label, y)

label_plot <- ggplot(label_df, aes(x = 1, y = y, label = row_label)) +
  geom_text(hjust = 1, size = 4) +
  annotate(
    "text", x = 1, y = 7.1,
    label = "No. at risk",
    hjust = 1, fontface = "bold", size = 4.4
  ) +
  xlim(0, 1.05) +
  ylim(0.3, 7.4) +
  coord_cartesian(clip = "off") +
  theme_void()

################################################################################
# 11. RISK TABLE NUMBERS
################################################################################

number_plot <- ggplot(
  risk_df,
  aes(x = time, y = y, label = n, color = Risk_Group)
) +
  geom_text(size = 4.6) +
  scale_color_manual(values = RISK_COLORS, guide = "none") +
  scale_x_continuous(
    breaks = risk_times,
    limits = c(-0.85, x_max + 0.20),
    expand = c(0, 0)
  ) +
  scale_y_continuous(limits = c(0.3, 7.4), breaks = NULL) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 13) +
  theme(
    axis.line = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title = element_blank()
  )

################################################################################
# 12. COMBINE
################################################################################

risk_grob <- arrangeGrob(
  ggplotGrob(label_plot),
  ggplotGrob(number_plot),
  ncol = 2,
  widths = c(2.15, 6)
)

final_plot <- arrangeGrob(
  ggplotGrob(main_plot),
  risk_grob,
  ncol = 1,
  heights = c(3.3, 1.3)
)

grid.newpage()
grid.draw(final_plot)

################################################################################
# 13. SAVE FIGURE
################################################################################

PNG_FILE <- file.path(OUT_DIR, "Fig5C_Clinical_vs_GEP_TTP_Overlay.png")
PDF_FILE <- file.path(OUT_DIR, "Fig5C_Clinical_vs_GEP_TTP_Overlay.pdf")

png(PNG_FILE, width = 11, height = 9.5, units = "in", res = 300, bg = "white")
grid.draw(final_plot)
dev.off()

pdf(PDF_FILE, width = 11, height = 9.5)
grid.draw(final_plot)
dev.off()

################################################################################
# 14. SAVE MATCHED COHORT
################################################################################

MATCHED_FILE <- file.path(OUT_DIR, "Fig5C_matched-analysis-cohort.csv")
write_csv(matched, MATCHED_FILE)

################################################################################
# 15. FINAL CHECK
################################################################################

cat(
  "\nStep 30 complete\n",
  "Matched n = ", N_MATCHED, "\n",
  "Clinical log-rank p = ", clinical_p, "\n",
  "GEP log-rank p = ", gep_p, "\n",
  "PNG: ", PNG_FILE, "\n",
  "PDF: ", PDF_FILE, "\n",
  "Matched cohort: ", MATCHED_FILE, "\n",
  sep = ""
)

cat("\nClinical groups:\n")
print(table(matched$Clinical_Risk))

cat("\nGEP groups:\n")
print(table(matched$GEP_Risk))