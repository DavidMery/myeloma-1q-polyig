################################################################################
# PREPARE COMPLETE UAMS EXPRESSIONSET
#
# Purpose:
#   1. Read the raw UAMS expression files.
#   2. Combine all expression datasets.
#   3. Match feature and phenotype data.
#   4. Create and save the complete UAMS ExpressionSet.
#
# Missing expression values are retained as NA.
# Raw input files are not modified.
#
# Required first:
#   source("scripts/00_prepare-environment.R")
################################################################################


################################################################################
# 1. Input files
################################################################################

expression_files <- c(
  "expression-rnas-mgus.csv",
  "expression-rnas-ndmm.csv",
  "expression-rnbx-mgus.csv",
  "expression-rnbx-ndmm.csv"
)

raw_dir <- file.path(
  "data",
  "raw",
  "uams"
)


################################################################################
# 2. Standardize CHIPID
################################################################################

standardize_chipid <- function(x) {
  
  x <- trimws(
    as.character(x)
  )
  
  # Standardize existing CHIPID prefix.
  x <- sub(
    "^chipid",
    "CHIPID",
    x,
    ignore.case = TRUE
  )
  
  # Add CHIPID prefix to numeric-only identifiers.
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


################################################################################
# 3. Read expression files
################################################################################

expression_matrices <- list()

for (file in expression_files) {
  
  expression <- fread(
    file.path(
      raw_dir,
      file
    ),
    na.strings = c(
      "",
      "NA",
      "N/A",
      "NaN",
      "NULL",
      ".",
      "#NUM!"
    )
  )
  
  probe_id <- expression$Probe_ID
  
  expression$Probe_ID <- NULL
  
  setnames(
    expression,
    standardize_chipid(
      names(expression)
    )
  )
  
  expression_matrix <- as.matrix(
    expression
  )
  
  storage.mode(
    expression_matrix
  ) <- "double"
  
  rownames(
    expression_matrix
  ) <- probe_id
  
  expression_matrices[[file]] <- expression_matrix
}


################################################################################
# 4. Align probes and combine expression matrices
################################################################################

reference_probes <- rownames(
  expression_matrices[[1]]
)

expression_matrices <- lapply(
  expression_matrices,
  function(x) {
    
    x[
      reference_probes,
      ,
      drop = FALSE
    ]
    
  }
)

expression_matrix <- do.call(
  cbind,
  expression_matrices
)


################################################################################
# 5. Read and match feature data
################################################################################

features <- fread(
  file.path(
    raw_dir,
    "features.csv"
  ),
  data.table = FALSE
)

features <- features[
  match(
    rownames(expression_matrix),
    features$Probe_ID
  ),
  ,
  drop = FALSE
]

rownames(features) <- features$Probe_ID


################################################################################
# 6. Read and match phenotype data
################################################################################

phenotype <- fread(
  file.path(
    raw_dir,
    "phenotype.csv"
  ),
  data.table = FALSE
)

phenotype$chipid <- standardize_chipid(
  phenotype$chipid
)

phenotype <- phenotype[
  match(
    colnames(expression_matrix),
    phenotype$chipid
  ),
  ,
  drop = FALSE
]

rownames(phenotype) <- phenotype$chipid


################################################################################
# 7. Create ExpressionSet
################################################################################

eset <- ExpressionSet(
  assayData = expression_matrix,
  phenoData = AnnotatedDataFrame(
    phenotype
  ),
  featureData = AnnotatedDataFrame(
    features
  )
)


################################################################################
# 8. Save ExpressionSet
################################################################################

dir.create(
  file.path(
    "data",
    "processed",
    "uams"
  ),
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  eset,
  file.path(
    "data",
    "processed",
    "uams",
    "ESET_uams_all.rds"
  )
)

cat(
  "ESET_uams_all.rds created:",
  nrow(eset),
  "probes x",
  ncol(eset),
  "samples\n"
)