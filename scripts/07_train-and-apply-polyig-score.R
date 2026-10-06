################################################################################
# STEP 07 — TRAIN AND APPLY POLYIG SCORE
#
# Purpose:
#   1. Select RNAS CD138+ NDMM samples for PolyIG discovery.
#   2. Use IG_ARD as the supervised target.
#   3. Split eligible samples into training and held-out test sets.
#   4. Select probes with PCC > 0.40 and FDR < 0.01.
#   5. Retain one highest-PCC probe per gene.
#   6. Define PolyIG_Score as mean expression across selected probes.
#   7. Evaluate the frozen score in the held-out test set.
#   8. Save the PolyIG model and selected probe/gene lists.
#   9. Apply PolyIG_Score to the complete ExpressionSet.
################################################################################

library(dplyr)
library(caret)

set.seed(123)


################################################################################
# 1. Settings
################################################################################

TRAIN_FRAC <- 2 / 3

PCC_THRESHOLD <- 0.40
FDR_THRESHOLD <- 0.01

REMOVE_CHIPID <- "CHIPID5259625"

SCORE_METHOD <- "mean"


################################################################################
# 2. Files
################################################################################

ESET_FILE <- file.path(
  "data",
  "processed",
  "uams",
  "ESET_uams_all_gep1q_fish.rds"
)

MODEL_DIR <- file.path(
  "models",
  "uams",
  "polyig"
)

MODEL_FILE <- file.path(
  MODEL_DIR,
  "polyig-model.rds"
)

OUTPUT_FILE <- file.path(
  "data",
  "processed",
  "uams",
  "ESET_uams_all_gep1q_fish_polyig.rds"
)

dir.create(
  MODEL_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


################################################################################
# 3. Load ExpressionSet
################################################################################

eset <- readRDS(
  ESET_FILE
)

expression <- exprs(
  eset
)

phenotype <- pData(
  eset
)

features <- fData(
  eset
)


################################################################################
# 4. Select RNAS CD138+ NDMM derivation cohort
################################################################################

keep <- phenotype$sample_group == "RNAS_CD138_NDMM"

keep[is.na(keep)] <- FALSE

derivation_ids <- rownames(
  phenotype
)[
  keep
]

expression_derivation <- expression[
  ,
  derivation_ids,
  drop = FALSE
]

phenotype_derivation <- phenotype[
  derivation_ids,
  ,
  drop = FALSE
]

cat(
  "RNAS CD138+ NDMM samples:",
  length(derivation_ids),
  "\n"
)


################################################################################
# 5. Remove specified sample
################################################################################

keep <- rownames(
  phenotype_derivation
) != REMOVE_CHIPID

phenotype_derivation <- phenotype_derivation[
  keep,
  ,
  drop = FALSE
]

expression_derivation <- expression_derivation[
  ,
  rownames(phenotype_derivation),
  drop = FALSE
]


################################################################################
# 6. Keep samples with valid IG_ARD
################################################################################

phenotype_derivation$IG_ARD <- as.numeric(
  phenotype_derivation$IG_ARD
)

keep <- is.finite(
  phenotype_derivation$IG_ARD
)

phenotype_derivation <- phenotype_derivation[
  keep,
  ,
  drop = FALSE
]

expression_derivation <- expression_derivation[
  ,
  rownames(phenotype_derivation),
  drop = FALSE
]


################################################################################
# 7. Create training and held-out test sets
################################################################################

training_index <- createDataPartition(
  phenotype_derivation$IG_ARD,
  p = TRAIN_FRAC,
  list = FALSE
)

training_index <- sort(
  unique(
    as.integer(
      training_index
    )
  )
)

test_index <- setdiff(
  seq_len(
    nrow(phenotype_derivation)
  ),
  training_index
)

train_ids <- rownames(
  phenotype_derivation
)[
  training_index
]

test_ids <- rownames(
  phenotype_derivation
)[
  test_index
]

x_train <- expression_derivation[
  ,
  train_ids,
  drop = FALSE
]

y_train <- phenotype_derivation[
  train_ids,
  "IG_ARD"
]

cat(
  "Training samples:",
  length(train_ids),
  "\n"
)

cat(
  "Held-out samples:",
  length(test_ids),
  "\n"
)


################################################################################
# 8. Calculate probe-level correlation with IG_ARD
################################################################################

pcc <- apply(
  x_train,
  1,
  function(x) {
    cor(
      x,
      y_train,
      method = "pearson"
    )
  }
)

p_value <- apply(
  x_train,
  1,
  function(x) {
    cor.test(
      x,
      y_train,
      method = "pearson"
    )$p.value
  }
)

fdr <- p.adjust(
  p_value,
  method = "BH"
)


################################################################################
# 9. Select positively correlated probes
################################################################################

probe_results <- data.frame(
  Probe_ID = rownames(x_train),
  GENE = trimws(
    as.character(
      features[
        rownames(x_train),
        "GENE"
      ]
    )
  ),
  PCC = pcc,
  FDR = fdr
)

probe_results <- probe_results |>
  filter(
    !is.na(GENE),
    GENE != "",
    tolower(GENE) != "no-gene",
    PCC > PCC_THRESHOLD,
    FDR < FDR_THRESHOLD
  )


################################################################################
# 10. Keep one highest-PCC probe per gene
################################################################################

selected_probes <- probe_results |>
  arrange(
    GENE,
    desc(PCC),
    FDR,
    Probe_ID
  ) |>
  group_by(
    GENE
  ) |>
  slice_head(
    n = 1
  ) |>
  ungroup() |>
  arrange(
    desc(PCC)
  )

cat(
  "Selected PolyIG probes:",
  nrow(selected_probes),
  "\n"
)


################################################################################
# 11. Calculate PolyIG_Score for all samples
################################################################################

polyig_expression <- expression[
  selected_probes$Probe_ID,
  ,
  drop = FALSE
]

PolyIG_Score <- colMeans(
  polyig_expression,
  na.rm = FALSE
)

names(
  PolyIG_Score
) <- colnames(
  polyig_expression
)


################################################################################
# 12. Evaluate held-out test set
################################################################################

heldout_r <- cor(
  phenotype_derivation[
    test_ids,
    "IG_ARD"
  ],
  PolyIG_Score[
    test_ids
  ],
  method = "pearson"
)

cat(
  "Held-out Pearson r:",
  round(
    heldout_r,
    3
  ),
  "\n"
)


################################################################################
# 13. Save frozen PolyIG model
################################################################################

polyig_model <- list(
  selected_probes = selected_probes$Probe_ID,
  selected_genes = selected_probes$GENE,
  score_method = SCORE_METHOD,
  pcc_threshold = PCC_THRESHOLD,
  fdr_threshold = FDR_THRESHOLD,
  training_chipids = train_ids,
  heldout_test_chipids = test_ids
)

saveRDS(
  polyig_model,
  MODEL_FILE
)


################################################################################
# 14. Save selected probes and genes
################################################################################

write.csv(
  selected_probes,
  file.path(
    MODEL_DIR,
    "polyig-selected-probes.csv"
  ),
  row.names = FALSE
)

writeLines(
  selected_probes$Probe_ID,
  file.path(
    MODEL_DIR,
    "polyig-selected-probe-ids.txt"
  )
)

writeLines(
  selected_probes$GENE,
  file.path(
    MODEL_DIR,
    "polyig-selected-gene-symbols.txt"
  )
)


################################################################################
# 15. Add PolyIG_Score to phenotype data
################################################################################

phenotype$PolyIG_Score <- PolyIG_Score[
  rownames(phenotype)
]


################################################################################
# 16. Update ExpressionSet
################################################################################

phenoData(
  eset
) <- AnnotatedDataFrame(
  phenotype
)


################################################################################
# 17. Save PolyIG-scored ExpressionSet
################################################################################

saveRDS(
  eset,
  OUTPUT_FILE
)


################################################################################
# 18. Complete
################################################################################

cat(
  "Selected PolyIG genes:",
  nrow(selected_probes),
  "\n"
)

cat(
  "Samples with PolyIG_Score:",
  sum(
    is.finite(
      phenotype$PolyIG_Score
    )
  ),
  "\n"
)

cat(
  "PolyIG model saved:",
  MODEL_FILE,
  "\n"
)

cat(
  "Selected probe IDs saved:",
  file.path(
    MODEL_DIR,
    "polyig-selected-probe-ids.txt"
  ),
  "\n"
)

cat(
  "Selected gene symbols saved:",
  file.path(
    MODEL_DIR,
    "polyig-selected-gene-symbols.txt"
  ),
  "\n"
)

cat(
  "ExpressionSet saved:",
  OUTPUT_FILE,
  "\n"
)