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
#
# Required first:
#   Run through code/run_pipeline.R so DATA_ROOT and MODEL_ROOT are defined.
################################################################################

library(dplyr)
library(caret)

set.seed(123)


################################################################################
# 1. VERIFY PIPELINE ROOTS
################################################################################

if (!exists("DATA_ROOT")) {
  
  stop(
    "DATA_ROOT is not defined. ",
    "Run this script through code/run_pipeline.R."
  )
}

if (!exists("MODEL_ROOT")) {
  
  stop(
    "MODEL_ROOT is not defined. ",
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
  MODEL_ROOT,
  recursive = TRUE,
  showWarnings = FALSE
)


################################################################################
# 2. SETTINGS
################################################################################

TRAIN_FRAC <- 2 / 3

PCC_THRESHOLD <- 0.40
FDR_THRESHOLD <- 0.01

REMOVE_CHIPID <- "CHIPID5259625"

SCORE_METHOD <- "mean"


################################################################################
# 3. FILES
################################################################################

ESET_FILE <- file.path(
  DATA_ROOT,
  "processed",
  "ESET_uams_all_gep1q_fish.rds"
)

MODEL_DIR <- file.path(
  MODEL_ROOT,
  "polyig"
)

MODEL_FILE <- file.path(
  MODEL_DIR,
  "polyig-model.rds"
)

OUTPUT_FILE <- file.path(
  DATA_ROOT,
  "processed",
  "ESET_uams_all_gep1q_fish_polyig.rds"
)

if (!file.exists(ESET_FILE)) {
  
  stop(
    "FISH-annotated ExpressionSet not found: ",
    ESET_FILE
  )
}

dir.create(
  MODEL_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  dirname(OUTPUT_FILE),
  recursive = TRUE,
  showWarnings = FALSE
)


################################################################################
# 4. LOAD EXPRESSIONSET
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
# 5. SELECT RNAS CD138+ NDMM DERIVATION COHORT
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
# 6. REMOVE SPECIFIED SAMPLE
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
# 7. KEEP SAMPLES WITH VALID IG_ARD
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
# 8. CREATE TRAINING AND HELD-OUT TEST SETS
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
# 9. CALCULATE PROBE-LEVEL CORRELATION WITH IG_ARD
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
# 10. SELECT POSITIVELY CORRELATED PROBES
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
# 11. KEEP ONE HIGHEST-PCC PROBE PER GENE
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
# 12. CALCULATE POLYIG_SCORE FOR ALL SAMPLES
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
# 13. EVALUATE HELD-OUT TEST SET
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
# 14. SAVE FROZEN POLYIG MODEL
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
# 15. SAVE SELECTED PROBES AND GENES
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
# 16. ADD POLYIG_SCORE TO PHENOTYPE DATA
################################################################################

phenotype$PolyIG_Score <- PolyIG_Score[
  rownames(phenotype)
]


################################################################################
# 17. UPDATE EXPRESSIONSET
################################################################################

phenoData(
  eset
) <- AnnotatedDataFrame(
  phenotype
)


################################################################################
# 18. SAVE POLYIG-SCORED EXPRESSIONSET
################################################################################

saveRDS(
  eset,
  OUTPUT_FILE
)


################################################################################
# 19. COMPLETE
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