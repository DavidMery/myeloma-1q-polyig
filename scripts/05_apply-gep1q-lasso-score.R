################################################################################
# STEP 05 — APPLY GEP1q ELASTIC-NET SCORE
#
# Purpose:
#   1. Load the complete UAMS ExpressionSet.
#   2. Reconstruct chromosome 1q gene-level expression as used in training.
#   3. Apply the frozen elastic-net model to all RNAS samples.
#   4. Add:
#        GEP1qscore_lasso
#        GEP1qcopy_lasso
#        GEP1q_binary_lasso
#   5. Save the scored ExpressionSet.
################################################################################

library(glmnet)


################################################################################
# 1. Files
################################################################################

ESET_FILE <- file.path(
  "data",
  "processed",
  "uams",
  "ESET_uams_all.rds"
)

MODEL_FILE <- file.path(
  "models",
  "uams",
  "gep1q-penalized",
  "elastic-net-alpha-0p50",
  "gep1q-elastic-net-alpha-0p50-model.rds"
)

OUTPUT_FILE <- file.path(
  "data",
  "processed",
  "uams",
  "ESET_uams_all_gep1q_lasso.rds"
)


################################################################################
# 2. Load ExpressionSet and frozen model
################################################################################

eset <- readRDS(
  ESET_FILE
)

model <- readRDS(
  MODEL_FILE
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
# 3. Select all RNAS samples
################################################################################

score_samples <- grepl(
  "^RNAS_",
  phenotype$sample_group
)

score_samples[is.na(score_samples)] <- FALSE

score_ids <- rownames(
  phenotype
)[
  score_samples
]

cat(
  "RNAS samples to score:",
  length(score_ids),
  "\n"
)


################################################################################
# 4. Restrict expression to chromosome 1q
################################################################################

probe_annotation <- data.frame(
  Probe_ID = rownames(features),
  GENE = trimws(
    as.character(features$GENE)
  ),
  CHR_ARM = tolower(
    gsub(
      "[[:space:]]+",
      "",
      trimws(
        as.character(features$CHR_ARM)
      )
    )
  )
)

probe_annotation <- probe_annotation[
  probe_annotation$CHR_ARM == "1q",
  ,
  drop = FALSE
]

expression_1q <- expression[
  probe_annotation$Probe_ID,
  score_ids,
  drop = FALSE
]


################################################################################
# 5. Collapse probes to gene-level mean expression
################################################################################

keep_genes <- !is.na(probe_annotation$GENE) &
  probe_annotation$GENE != "" &
  tolower(probe_annotation$GENE) != "no-gene"

probe_annotation <- probe_annotation[
  keep_genes,
  ,
  drop = FALSE
]

expression_1q <- expression_1q[
  keep_genes,
  ,
  drop = FALSE
]

gene_sum <- rowsum(
  expression_1q,
  group = probe_annotation$GENE
)

gene_count <- table(
  probe_annotation$GENE
)[
  rownames(gene_sum)
]

gene_expression <- sweep(
  gene_sum,
  1,
  as.numeric(gene_count),
  "/"
)


################################################################################
# 6. Keep the exact candidate genes used during training
################################################################################

gene_expression <- gene_expression[
  model$candidate_genes,
  ,
  drop = FALSE
]

x_apply <- t(
  gene_expression
)


################################################################################
# 7. Apply frozen GEP1q model
################################################################################

gep1q_score <- as.numeric(
  predict(
    model$fit,
    newx = x_apply,
    s = model$lambda,
    type = "response"
  )
)

names(
  gep1q_score
) <- score_ids


################################################################################
# 8. Convert continuous score to predicted 1q class
################################################################################

gep1q_copy <- ifelse(
  gep1q_score >= model$cutoff,
  "1q3plus",
  "1q2"
)

gep1q_binary <- ifelse(
  gep1q_score >= model$cutoff,
  1L,
  0L
)


################################################################################
# 9. Add GEP1q fields to phenotype data
################################################################################

phenotype$GEP1qscore_lasso <- NA_real_
phenotype$GEP1qcopy_lasso <- NA_character_
phenotype$GEP1q_binary_lasso <- NA_integer_

phenotype[
  score_ids,
  "GEP1qscore_lasso"
] <- gep1q_score

phenotype[
  score_ids,
  "GEP1qcopy_lasso"
] <- gep1q_copy

phenotype[
  score_ids,
  "GEP1q_binary_lasso"
] <- gep1q_binary


################################################################################
# 10. Update ExpressionSet
################################################################################

phenoData(
  eset
) <- AnnotatedDataFrame(
  phenotype
)


################################################################################
# 11. Save scored ExpressionSet
################################################################################

saveRDS(
  eset,
  OUTPUT_FILE
)


################################################################################
# 12. Complete
################################################################################

cat(
  "Samples scored:",
  length(score_ids),
  "\n"
)

cat(
  "Predicted 1q2:",
  sum(
    gep1q_copy == "1q2"
  ),
  "\n"
)

cat(
  "Predicted 1q3plus:",
  sum(
    gep1q_copy == "1q3plus"
  ),
  "\n"
)

cat(
  "Saved:",
  OUTPUT_FILE,
  "\n"
)