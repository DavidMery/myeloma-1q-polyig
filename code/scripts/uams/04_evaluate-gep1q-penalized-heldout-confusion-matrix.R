################################################################################
# STEP 04 — EVALUATE GEP1q HELD-OUT TEST SET
#
# Purpose:
#   1. Load the frozen held-out predictions from Step 03.
#   2. Compare observed FISH class with predicted GEP1q class.
#   3. Calculate classification performance.
#   4. Create and save the confusion-matrix figure.
#
# Required first:
#   Run through code/run_pipeline.R so MODEL_ROOT is defined.
################################################################################

library(ggplot2)
library(dplyr)


################################################################################
# 1. Verify model root
################################################################################

if (!exists("MODEL_ROOT")) {
  stop(
    "MODEL_ROOT is not defined. ",
    "Run this script through code/run_pipeline.R."
  )
}

dir.create(
  MODEL_ROOT,
  recursive = TRUE,
  showWarnings = FALSE
)


################################################################################
# 2. Files
################################################################################

MODEL_DIR <- file.path(
  MODEL_ROOT,
  "gep1q-penalized",
  "elastic-net-alpha-0p50"
)

TEST_FILE <- file.path(
  MODEL_DIR,
  "gep1q-elastic-net-alpha-0p50-heldout-test-predictions.csv"
)

OUT_DIR <- MODEL_DIR

if (!file.exists(TEST_FILE)) {
  stop(
    "Held-out prediction file not found: ",
    TEST_FILE
  )
}

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


################################################################################
# 3. Load held-out predictions
################################################################################

test <- fread(
  TEST_FILE,
  data.table = FALSE
)


################################################################################
# 4. Standardize class labels
################################################################################

test$Observed <- ifelse(
  test$Test_Class == "2",
  "1q2",
  "1q≥3"
)

test$Predicted <- ifelse(
  test$Predicted_Class == "2",
  "1q2",
  "1q≥3"
)

test$Observed <- factor(
  test$Observed,
  levels = c(
    "1q2",
    "1q≥3"
  )
)

test$Predicted <- factor(
  test$Predicted,
  levels = c(
    "1q2",
    "1q≥3"
  )
)


################################################################################
# 5. Confusion matrix
################################################################################

confusion <- table(
  Observed = test$Observed,
  Predicted = test$Predicted
)

print(
  confusion
)


################################################################################
# 6. Classification performance
################################################################################

TN <- confusion["1q2", "1q2"]
FP <- confusion["1q2", "1q≥3"]
FN <- confusion["1q≥3", "1q2"]
TP <- confusion["1q≥3", "1q≥3"]

accuracy <- (
  TP + TN
) / sum(confusion)

sensitivity <- TP / (
  TP + FN
)

specificity <- TN / (
  TN + FP
)

PPV <- TP / (
  TP + FP
)

NPV <- TN / (
  TN + FN
)

performance <- data.frame(
  Metric = c(
    "Accuracy",
    "Sensitivity",
    "Specificity",
    "Positive predictive value",
    "Negative predictive value"
  ),
  Value = c(
    accuracy,
    sensitivity,
    specificity,
    PPV,
    NPV
  )
)

print(
  performance
)


################################################################################
# 7. Prepare confusion-matrix figure
################################################################################

confusion_percent <- prop.table(
  confusion,
  margin = 1
) * 100

plot_data <- as.data.frame(
  confusion
)

plot_data$Percent <- as.vector(
  confusion_percent
)

plot_data$Label <- paste0(
  "n = ",
  plot_data$Freq,
  "\n(",
  sprintf(
    "%.1f",
    plot_data$Percent
  ),
  "%)"
)


################################################################################
# 8. Create confusion-matrix figure
################################################################################

p <- ggplot(
  plot_data,
  aes(
    x = Predicted,
    y = Observed,
    fill = Percent
  )
) +
  geom_tile(
    color = "white",
    linewidth = 1.5
  ) +
  geom_text(
    aes(
      label = Label
    ),
    size = 6,
    fontface = "bold"
  ) +
  scale_fill_gradient(
    low = "white",
    high = "steelblue",
    limits = c(
      0,
      100
    ),
    name = "% within\nFISH class"
  ) +
  coord_equal() +
  labs(
    title = "Strict Held-Out GEP1q Classification",
    subtitle = paste0(
      "Accuracy = ",
      sprintf(
        "%.1f",
        100 * accuracy
      ),
      "%"
    ),
    x = "Predicted GEP1q class",
    y = "Observed FISH 1q class"
  ) +
  theme_classic(
    base_size = 16
  ) +
  theme(
    axis.title = element_text(
      face = "bold"
    ),
    axis.text = element_text(
      color = "black"
    )
  )

print(
  p
)


################################################################################
# 9. Save results
################################################################################

write.csv(
  performance,
  file.path(
    OUT_DIR,
    "heldout-classification-performance.csv"
  ),
  row.names = FALSE
)

write.csv(
  as.data.frame(confusion),
  file.path(
    OUT_DIR,
    "heldout-confusion-matrix.csv"
  ),
  row.names = FALSE
)

ggsave(
  file.path(
    OUT_DIR,
    "heldout-confusion-matrix.png"
  ),
  p,
  width = 7,
  height = 6,
  dpi = 600
)

ggsave(
  file.path(
    OUT_DIR,
    "heldout-confusion-matrix.pdf"
  ),
  p,
  width = 7,
  height = 6
)


################################################################################
# 10. Complete
################################################################################

cat(
  "Held-out samples:",
  nrow(test),
  "\n"
)

cat(
  "Accuracy:",
  round(
    100 * accuracy,
    1
  ),
  "%\n"
)

cat(
  "Sensitivity:",
  round(
    100 * sensitivity,
    1
  ),
  "%\n"
)

cat(
  "Specificity:",
  round(
    100 * specificity,
    1
  ),
  "%\n"
)