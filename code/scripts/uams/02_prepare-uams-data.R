################################################################################
# STEP 02 — PREPARE COMPLETE UAMS EXPRESSIONSET
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
# Data location:
#   Controlled by DATA_ROOT in code/run_pipeline.R.
#
#   PIPELINE_MODE = "EXAMPLE"
#     -> GitHub/myeloma-1q-polyig/data/simulated/uams
#
#   PIPELINE_MODE = "UAMS"
#     -> Box/PROJECTS/myeloma-1q-polyig/data
################################################################################


################################################################################
# 1. VERIFY DATA ROOT
################################################################################

if (!exists("DATA_ROOT")) {
  
  stop(
    "DATA_ROOT is not defined. ",
    "Run this script through code/run_pipeline.R."
  )
  
}

if (!dir.exists(DATA_ROOT)) {
  
  stop(
    "DATA_ROOT does not exist: ",
    DATA_ROOT
  )
  
}


################################################################################
# 2. DEFINE INPUT AND OUTPUT DIRECTORIES
################################################################################

raw_dir <- file.path(
  DATA_ROOT,
  "raw"
)

processed_dir <- file.path(
  DATA_ROOT,
  "processed"
)

if (!dir.exists(raw_dir)) {
  
  stop(
    "Raw UAMS data directory not found: ",
    raw_dir
  )
  
}

dir.create(
  processed_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


################################################################################
# 3. INPUT FILES
################################################################################

expression_files <- c(
  "expression-rnas-mgus.csv",
  "expression-rnas-ndmm.csv",
  "expression-rnbx-mgus.csv",
  "expression-rnbx-ndmm.csv"
)


################################################################################
# 4. VERIFY REQUIRED INPUT FILES
################################################################################

required_files <- c(
  expression_files,
  "features.csv",
  "phenotype.csv"
)

required_paths <- file.path(
  raw_dir,
  required_files
)

missing_files <- required_paths[
  !file.exists(required_paths)
]

if (length(missing_files) > 0L) {
  
  stop(
    "Required input file(s) not found:\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )
  
}


################################################################################
# 5. STANDARDIZE CHIPID
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
# 6. READ EXPRESSION FILES
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
# 7. ALIGN PROBES AND COMBINE EXPRESSION MATRICES
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
# 8. READ AND MATCH FEATURE DATA
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
# 9. READ AND MATCH PHENOTYPE DATA
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
# 10. CREATE EXPRESSIONSET
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
# 11. SAVE EXPRESSIONSET
################################################################################

output_file <- file.path(
  processed_dir,
  "ESET_uams_all.rds"
)

saveRDS(
  eset,
  output_file
)


################################################################################
# 12. COMPLETE
################################################################################

cat(
  "\n",
  "Step 02 completed successfully.\n",
  "ESET_uams_all.rds created:\n",
  "  ",
  nrow(eset),
  " probes x ",
  ncol(eset),
  " samples\n",
  "Saved to:\n",
  "  ",
  output_file,
  "\n",
  sep = ""
)