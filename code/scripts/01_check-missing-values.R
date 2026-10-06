################################################################################
# CHECK UAMS EXPRESSION FILES FOR MISSING VALUES
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
################################################################################


files <- c(
  "expression-rnas-mgus.csv",
  "expression-rnas-ndmm.csv",
  "expression-rnbx-mgus.csv",
  "expression-rnbx-ndmm.csv"
)


for (file in files) {
  
  x <- fread(
    file.path("data", "raw", "uams", file),
    colClasses = "character",
    na.strings = NULL
  )
  
  x <- as.matrix(
    x[, setdiff(names(x), "Probe_ID"), with = FALSE]
  )
  
  x <- trimws(x)
  
  n_missing <- sum(
    is.na(x) |
      x %in% c("#NUM!", "NA", "N/A", "NaN", "NULL", ".", "")
  )
  
  cat(
    file,
    ":",
    n_missing,
    "missing expression values\n"
  )
}