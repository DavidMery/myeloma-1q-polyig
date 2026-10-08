################################################################################
# STEP 01 — CHECK UAMS EXPRESSION FILES FOR MISSING VALUES
#
# Purpose:
#   Search each raw UAMS expression file for common missing-value entries.
#
# Missing values checked:
#   "#NUM!", "NA", "N/A", "NaN", "NULL", ".", and blank values.
#
# Only expression measurements are checked.
# Probe_ID is excluded.
# No data are modified.
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
# 2. INPUT FILES
################################################################################

files <- c(
  "expression-rnas-mgus.csv",
  "expression-rnas-ndmm.csv",
  "expression-rnbx-mgus.csv",
  "expression-rnbx-ndmm.csv"
)


################################################################################
# 3. CHECK EACH EXPRESSION FILE
################################################################################

for (file in files) {
  
  input_file <- file.path(
    DATA_ROOT,
    "raw",
    file
  )
  
  if (!file.exists(input_file)) {
    stop(
      "Expression file not found: ",
      input_file
    )
  }
  
  x <- fread(
    input_file,
    colClasses = "character",
    na.strings = NULL
  )
  
  x <- as.matrix(
    x[, setdiff(names(x), "Probe_ID"), with = FALSE]
  )
  
  x <- trimws(x)
  
  n_missing <- sum(
    is.na(x) |
      x %in% c(
        "#NUM!",
        "NA",
        "N/A",
        "NaN",
        "NULL",
        ".",
        ""
      )
  )
  
  cat(
    file,
    ": ",
    n_missing,
    " missing expression values\n",
    sep = ""
  )
}


################################################################################
# 4. COMPLETE
################################################################################

cat(
  "\nStep 01 completed successfully.\n"
)