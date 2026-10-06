################################################################################
# STEP 03 — TRAIN GEP1q ELASTIC-NET MODEL
#
# Purpose:
#   1. Select baseline CD138+ NDMM samples from TT2-TT5.
#   2. Restrict expression to chromosome 1q genes.
#   3. Define nonoverlapping training and strict held-out FISH cohorts.
#   4. Average multiple probes belonging to the same gene.
#   5. Fit an alpha = 0.50 elastic-net model using 10-fold CV.
#   6. Apply the frozen model to the strict held-out cohort.
#   7. Save the frozen model, selected genes, and held-out predictions.
#
# Required first:
#   source("scripts/00_prepare-environment.R")
#   source("scripts/02_prepare-uams-data.R")
################################################################################

library(glmnet)
library(dplyr)


################################################################################
# 1. Settings
################################################################################

set.seed(123)

TRAIN_PURITY_MIN <- 80
TEST_PURITY_MIN <- 95
TEST_DOMINANCE_PCT <- 90

MODEL_ALPHA <- 0.50
CV_FOLDS <- 10

FISH_CLASSIFICATION_CUTOFF <- 2.5

TT_KEEP <- c(
  "TT2_NoThal",
  "TT2_Thal",
  "TT3a",
  "TT3b",
  "TT4_S-TT3",
  "TT5"
)


################################################################################
# 2. Files
################################################################################

ESET_PATH <- file.path(
  "data",
  "processed",
  "uams",
  "ESET_uams_all.rds"
)

FISH_FILE <- file.path(
  "data",
  "raw",
  "uams",
  "fish-1q-cell-percentages.csv"
)

MODEL_DIR <- file.path(
  "models",
  "uams",
  "gep1q-penalized",
  "elastic-net-alpha-0p50"
)

dir.create(
  MODEL_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


################################################################################
# 3. Helper functions
################################################################################

normalize_chipid <- function(x) {
  
  x <- toupper(
    trimws(
      as.character(x)
    )
  )
  
  numeric_only <- grepl(
    "^[0-9]+$",
    x
  )
  
  x[numeric_only] <- paste0(
    "CHIPID",
    x[numeric_only]
  )
  
  x
}


make_stratified_folds <- function(group, nfolds) {
  
  foldid <- integer(
    length(group)
  )
  
  for (level in unique(group)) {
    
    index <- which(
      group == level
    )
    
    index <- sample(
      index
    )
    
    foldid[index] <- rep(
      1:nfolds,
      length.out = length(index)
    )
  }
  
  foldid
}


################################################################################
# 4. Load UAMS ExpressionSet
################################################################################

eset <- readRDS(
  ESET_PATH
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
# 5. Select baseline CD138+ NDMM from TT2-TT5
################################################################################

keep_samples <- phenotype$sample_group == "RNAS_CD138_NDMM" &
  phenotype$tt_on_for_tx1 %in% TT_KEEP

keep_samples[is.na(keep_samples)] <- FALSE

phenotype <- phenotype[
  keep_samples,
  ,
  drop = FALSE
]

expression <- expression[
  ,
  rownames(phenotype),
  drop = FALSE
]

cat(
  "Baseline CD138+ NDMM samples:",
  ncol(expression),
  "\n"
)


################################################################################
# 6. Restrict expression to chromosome 1q
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
  ,
  drop = FALSE
]

cat(
  "Chromosome 1q probes:",
  nrow(expression_1q),
  "\n"
)


################################################################################
# 7. Read FISH 1q cell-percentage data
################################################################################

fish_raw <- fread(
  FISH_FILE,
  data.table = FALSE
)

fish_raw$CHIPID <- normalize_chipid(
  fish_raw$CHIPID
)

get_copy_percent <- function(copy_number) {
  
  as.numeric(
    fish_raw[[
      paste0(
        "1q21_Percent_Cells_Copy_",
        copy_number
      )
    ]]
  )
}

fish <- data.frame(
  CHIPID = fish_raw$CHIPID,
  Copy0 = get_copy_percent(0),
  Copy1 = get_copy_percent(1),
  Copy2 = get_copy_percent(2),
  Copy3 = get_copy_percent(3),
  Copy4 = get_copy_percent(4),
  Copy5 = get_copy_percent(5),
  Copy6 = get_copy_percent(6)
)

fish$Copy3plus <- (
  fish$Copy3 +
    fish$Copy4 +
    fish$Copy5 +
    fish$Copy6
)

fish$Weighted_1q_Copy <- (
  0 * fish$Copy0 +
    1 * fish$Copy1 +
    2 * fish$Copy2 +
    3 * fish$Copy3 +
    4 * fish$Copy4 +
    5 * fish$Copy5 +
    6 * fish$Copy6
) / 100


################################################################################
# 8. Add phenotype information to FISH data
################################################################################

fish$POST_SORT_PURITY <- as.numeric(
  phenotype[
    match(
      fish$CHIPID,
      rownames(phenotype)
    ),
    "POST_SORT_PURITY"
  ]
)

fish$TT_Protocol <- phenotype[
  match(
    fish$CHIPID,
    rownames(phenotype)
  ),
  "tt_on_for_tx1"
]

fish <- fish[
  fish$CHIPID %in% colnames(expression),
  ,
  drop = FALSE
]


################################################################################
# 9. Define strict held-out test cohort
################################################################################

fish$Test_Copy_State <- case_when(
  
  fish$POST_SORT_PURITY >= TEST_PURITY_MIN &
    fish$Copy2 >= TEST_DOMINANCE_PCT ~ "2",
  
  fish$POST_SORT_PURITY >= TEST_PURITY_MIN &
    fish$Copy3 >= TEST_DOMINANCE_PCT ~ "3",
  
  fish$POST_SORT_PURITY >= TEST_PURITY_MIN &
    fish$Copy4 >= TEST_DOMINANCE_PCT ~ "4",
  
  fish$POST_SORT_PURITY >= TEST_PURITY_MIN &
    fish$Copy5 >= TEST_DOMINANCE_PCT ~ "5",
  
  fish$POST_SORT_PURITY >= TEST_PURITY_MIN &
    fish$Copy6 >= TEST_DOMINANCE_PCT ~ "6",
  
  TRUE ~ NA_character_
)

fish$Test_Class <- case_when(
  fish$Test_Copy_State == "2" ~ "2",
  fish$Test_Copy_State %in% c(
    "3",
    "4",
    "5",
    "6"
  ) ~ "3plus",
  TRUE ~ NA_character_
)

test_fish <- fish |>
  filter(
    !is.na(Test_Class),
    is.finite(Weighted_1q_Copy)
  )

test_ids <- test_fish$CHIPID


################################################################################
# 10. Define training cohort
################################################################################

train_fish <- fish |>
  filter(
    POST_SORT_PURITY >= TRAIN_PURITY_MIN,
    !CHIPID %in% test_ids,
    is.finite(Weighted_1q_Copy)
  ) |>
  mutate(
    Train_Class = ifelse(
      Copy3plus >= 20,
      "1q3plus",
      "1q2"
    )
  )

train_ids <- train_fish$CHIPID

cat(
  "Training samples:",
  length(train_ids),
  "\n"
)

cat(
  "Strict held-out samples:",
  length(test_ids),
  "\n"
)


################################################################################
# 11. Collapse 1q probes to gene-level mean expression
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
# 11A. Define one representative probe per gene for downstream visualization
################################################################################
#
# The elastic-net model remains gene-level: all probes mapping to a gene are
# averaged before model fitting. This representative-probe step therefore
# does NOT change the fitted model, lambda, coefficients, or predictions.
#
# For downstream heatmaps, one probe per gene is selected using:
#   1. highest Pearson correlation with the gene-level mean across training
#      samples;
#   2. highest training-sample SD if tied;
#   3. smallest Probe_ID if still tied.
#
# This creates a reproducible one-probe-per-gene file for visualization.

representative_probe_table <- lapply(
  rownames(gene_expression),
  function(gene_symbol) {
    
    gene_probe_rows <- which(
      probe_annotation$GENE == gene_symbol
    )
    
    if (length(gene_probe_rows) == 0) {
      return(NULL)
    }
    
    gene_probe_ids <- probe_annotation$Probe_ID[
      gene_probe_rows
    ]
    
    probe_training_expression <- expression_1q[
      gene_probe_ids,
      train_ids,
      drop = FALSE
    ]
    
    gene_mean_training <- gene_expression[
      gene_symbol,
      train_ids
    ]
    
    probe_correlations <- apply(
      probe_training_expression,
      1,
      function(x) {
        
        ok <- is.finite(x) &
          is.finite(gene_mean_training)
        
        if (sum(ok) < 3) {
          return(NA_real_)
        }
        
        suppressWarnings(
          cor(
            x[ok],
            gene_mean_training[ok],
            method = "pearson"
          )
        )
      }
    )
    
    probe_sd <- apply(
      probe_training_expression,
      1,
      sd,
      na.rm = TRUE
    )
    
    candidate_table <- data.frame(
      GENE = gene_symbol,
      Probe_ID = gene_probe_ids,
      Correlation_to_GeneMean = as.numeric(
        probe_correlations
      ),
      Training_SD = as.numeric(
        probe_sd
      ),
      stringsAsFactors = FALSE
    )
    
    candidate_table <- candidate_table[
      is.finite(
        candidate_table$Correlation_to_GeneMean
      ),
      ,
      drop = FALSE
    ]
    
    if (nrow(candidate_table) == 0) {
      return(NULL)
    }
    
    candidate_table <- candidate_table[
      order(
        -candidate_table$Correlation_to_GeneMean,
        -candidate_table$Training_SD,
        candidate_table$Probe_ID
      ),
      ,
      drop = FALSE
    ]
    
    candidate_table[1, , drop = FALSE]
  }
)

representative_probe_table <- bind_rows(
  representative_probe_table
)

################################################################################
# 12. Remove genes with no variation in training samples
################################################################################

training_sd <- apply(
  gene_expression[
    ,
    train_ids,
    drop = FALSE
  ],
  1,
  sd
)

gene_expression <- gene_expression[
  is.finite(training_sd) &
    training_sd > 0,
  ,
  drop = FALSE
]

x_train <- t(
  gene_expression[
    ,
    train_ids,
    drop = FALSE
  ]
)

x_test <- t(
  gene_expression[
    ,
    test_ids,
    drop = FALSE
  ]
)

y_train <- train_fish$Weighted_1q_Copy

cat(
  "Candidate 1q genes:",
  ncol(x_train),
  "\n"
)


################################################################################
# 13. Create stratified cross-validation folds
################################################################################

foldid <- make_stratified_folds(
  train_fish$Train_Class,
  CV_FOLDS
)


################################################################################
# 14. Fit elastic-net model
################################################################################

cv_fit <- cv.glmnet(
  x = x_train,
  y = y_train,
  family = "gaussian",
  alpha = MODEL_ALPHA,
  nfolds = CV_FOLDS,
  foldid = foldid,
  type.measure = "mse",
  standardize = TRUE,
  intercept = TRUE,
  keep = TRUE
)

selected_lambda <- cv_fit$lambda.1se
lambda_used <- "lambda.1se"

coefficient_matrix <- as.matrix(
  coef(
    cv_fit,
    s = selected_lambda
  )
)

coefficient_table <- data.frame(
  term = rownames(coefficient_matrix),
  coefficient = coefficient_matrix[, 1]
)

selected_genes <- coefficient_table |>
  filter(
    term != "(Intercept)",
    coefficient != 0
  ) |>
  arrange(
    desc(
      abs(coefficient)
    )
  )

names(selected_genes)[1] <- "GENE"


################################################################################
# 15. Fall back to lambda.min if lambda.1se selects no genes
################################################################################

if (nrow(selected_genes) == 0) {
  
  selected_lambda <- cv_fit$lambda.min
  lambda_used <- "lambda.min"
  
  coefficient_matrix <- as.matrix(
    coef(
      cv_fit,
      s = selected_lambda
    )
  )
  
  coefficient_table <- data.frame(
    term = rownames(coefficient_matrix),
    coefficient = coefficient_matrix[, 1]
  )
  
  selected_genes <- coefficient_table |>
    filter(
      term != "(Intercept)",
      coefficient != 0
    ) |>
    arrange(
      desc(
        abs(coefficient)
      )
    )
  
  names(selected_genes)[1] <- "GENE"
}


################################################################################
# 16. Apply frozen model to held-out test cohort
################################################################################

test_prediction <- as.numeric(
  predict(
    cv_fit,
    newx = x_test,
    s = selected_lambda,
    type = "response"
  )
)

test_fish <- test_fish[
  match(
    rownames(x_test),
    test_fish$CHIPID
  ),
  ,
  drop = FALSE
]

test_fish$Predicted_Weighted_1q_Copy <- test_prediction

test_fish$Predicted_Class <- ifelse(
  test_fish$Predicted_Weighted_1q_Copy >=
    FISH_CLASSIFICATION_CUTOFF,
  "3plus",
  "2"
)


################################################################################
# 17. Held-out performance
################################################################################

heldout_r <- cor(
  test_fish$Weighted_1q_Copy,
  test_fish$Predicted_Weighted_1q_Copy,
  method = "pearson"
)

heldout_accuracy <- mean(
  test_fish$Test_Class ==
    test_fish$Predicted_Class
)

cat(
  "Held-out Pearson r:",
  round(
    heldout_r,
    3
  ),
  "\n"
)

cat(
  "Held-out accuracy:",
  round(
    heldout_accuracy * 100,
    1
  ),
  "%\n"
)


################################################################################
# 18. Save selected genes and held-out predictions
################################################################################

write.csv(
  selected_genes,
  file.path(
    MODEL_DIR,
    "gep1q-elastic-net-alpha-0p50-selected-genes-and-coefficients.csv"
  ),
  row.names = FALSE
)

writeLines(
  selected_genes$GENE,
  file.path(
    MODEL_DIR,
    "gep1q-elastic-net-alpha-0p50-selected-gene-symbols.txt"
  )
)

selected_probe_map <- representative_probe_table |>
  filter(
    GENE %in% selected_genes$GENE
  ) |>
  arrange(
    match(
      GENE,
      selected_genes$GENE
    )
  )

if (
  nrow(selected_probe_map) !=
  nrow(selected_genes)
) {
  missing_probe_genes <- setdiff(
    selected_genes$GENE,
    selected_probe_map$GENE
  )
  
  stop(
    "Could not assign one representative probe to every selected GEP1q gene: ",
    paste(
      missing_probe_genes,
      collapse = ", "
    )
  )
}

writeLines(
  selected_probe_map$Probe_ID,
  file.path(
    MODEL_DIR,
    "gep1q-elastic-net-alpha-0p50-selected-probe-ids.txt"
  )
)

write.csv(
  selected_probe_map,
  file.path(
    MODEL_DIR,
    "gep1q-elastic-net-alpha-0p50-selected-probe-map.csv"
  ),
  row.names = FALSE
)

write.csv(
  test_fish,
  file.path(
    MODEL_DIR,
    "gep1q-elastic-net-alpha-0p50-heldout-test-predictions.csv"
  ),
  row.names = FALSE
)


################################################################################
# 19. Save frozen model
################################################################################

gep1q_model <- list(
  fit = cv_fit,
  alpha = MODEL_ALPHA,
  lambda = selected_lambda,
  lambda_used = lambda_used,
  cutoff = FISH_CLASSIFICATION_CUTOFF,
  candidate_genes = colnames(x_train),
  selected_genes = selected_genes$GENE,
  selected_probe_map = selected_probe_map
)

saveRDS(
  gep1q_model,
  file.path(
    MODEL_DIR,
    "gep1q-elastic-net-alpha-0p50-model.rds"
  )
)


################################################################################
# 20. Complete
################################################################################

cat(
  "Lambda used:",
  lambda_used,
  "\n"
)

cat(
  "Selected genes:",
  nrow(selected_genes),
  "\n"
)

cat(
  "GEP1q elastic-net model complete.\n"
)