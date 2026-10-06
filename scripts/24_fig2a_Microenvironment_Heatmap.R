################################################################################
# STEP 24C — FIGURE 2A
# PAIRED-RNAS HEATMAP OF ELASTIC-NET GEP1Q AND POLYIG GENES
#
# Pairing:
#   RNAS_CD138_NDMM -> RNBX_WB_NDMM by PATID + disease type
#   RNAS_CD138_MGUS -> RNBX_WB_MGUS by PATID + disease type
#
# Heatmap columns:
#   Paired RNAS samples
#
# Heatmap rows:
#   1. PolyIG
#   2. GEP1q
#
# Top annotations:
#   GEP Groups
#   PolyIG
#   Supervised ME score
#
# No continuous GEP1q or GEP70 annotation bars.
################################################################################

suppressPackageStartupMessages({
  library(Biobase)
  library(limma)
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(tibble)
  library(grid)
})

set.seed(123)

################################################################################
# 1. CONFIGURATION
################################################################################

FIGURE_ID <- "Fig2A_Lasso_BestProbe"

ESET_PATH <- file.path(
  "data", "processed", "uams",
  "ESET_uams_all_gep1q_fish_polyig_pcme_binaries_lasso_groups_rnbx-matched.rds"
)

GEP1Q_PROBE_FILE <- file.path(
  "models", "uams", "gep1q-penalized", "elastic-net-alpha-0p50",
  "gep1q-elastic-net-alpha-0p50-selected-probe-ids.txt"
)

POLYIG_PROBE_FILE <- file.path(
  "models", "uams", "polyig",
  "polyig-selected-probe-ids.txt"
)

FIGURE_DIR <- file.path("results", "figures", "Fig2", FIGURE_ID)
SOURCE_DATA_DIR <- file.path(FIGURE_DIR, "source-data")
OUT_DIR <- file.path(FIGURE_DIR, paste0(FIGURE_ID, "_RNAS_GEP1q_PolyIG_Heatmap"))

dir.create(SOURCE_DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SAMPLE_ID_COL <- "chipid"
PATID_COL <- "patid"
SAMPLE_GROUP_COL <- "sample_group"

RNBX_GROUP_COL <- "GEP_RNBX_Groups_MGUS_Lasso"
RNAS_GROUP_COL <- "GEP_RNAS_Groups_MGUS_Lasso"

MATCH_STATUS_COL <- "RNAS_to_RNBX_PATID_match_MGUS_Lasso"
MATCH_TYPE_COL <- "RNAS_to_RNBX_match_type_Lasso"
MATCHED_STATUS_VALUE <- "matched_one_to_one"

PCBM_BX_COL <- "PCBmBx_All"
POLYIG_COL <- "PolyIG_Score"

GENE_COL <- "GENE"
PC_ME_OTHER_COL <- "PC_ME_Other"

RNAS_NDMM <- "RNAS_CD138_NDMM"
RNBX_NDMM <- "RNBX_WB_NDMM"

RNAS_MGUS <- "RNAS_CD138_MGUS"
RNBX_MGUS <- "RNBX_WB_MGUS"

MGUS_GROUP_VALUE <- "MGUS"
GEP70_HIGH_GROUP_VALUE <- "GEP70>0.66"

MGUS_PCBMBX_MAX <- 100
GEP70_PCBMBX_MIN <- 0
GEP70_PCBMBX_MAX_EXCLUSIVE <- 70

# Must match Step 23C
MGUS_ME_FDR_CUTOFF <- 5e-6
GEP70_ME_FDR_CUTOFF <- 1e-2

USE_SUPERVISED_ME_HIGH_MINUS_LOW_SCORE <- TRUE
ORDER_COLUMNS_BY_SUPERVISED_ME_SCORE <- TRUE
SPLIT_SUPERVISED_ORDER_BY_RNBX_GROUP <- TRUE
SUPERVISED_ME_SCORE_LEFT_TO_RIGHT <- "LowToHigh"

HEATMAP_INCLUDE_ALL_NONMISSING_RNBX_GROUPS <- TRUE

HEATMAP_GROUP_ORDER <- c(
  "MGUS",
  "GEP1q-CN2/PolyIG-high",
  "GEP1q-CN>=3/PolyIG-high",
  "GEP1q-CN2/PolyIG-low",
  "GEP1q-CN>=3/PolyIG-low",
  "GEP70>0.66"
)

CUSTOM_RNBX_GROUP_COLORS <- c(
  "MGUS" = "yellow",
  "GEP1q-CN2/PolyIG-high" = "steelblue",
  "GEP1q-CN2/PolyIG-low" = "grey50",
  "GEP1q-CN>=3/PolyIG-high" = "#8B4513",
  "GEP1q-CN>=3/PolyIG-low" = "red",
  "GEP70>0.66" = "black"
)

SHOW_ROW_NAMES <- TRUE
SHOW_COLUMN_NAMES <- FALSE

HEATMAP_SCALE_METHOD <- "Zscore"
Z_SCORE_COLOR_LIMIT <- 2

POLYIG_COLOR_QUANTILES <- c(0.20, 0.80)

FIGURE_WIDTH_IN <- 20
FIGURE_HEIGHT_IN <- 10
PNG_RESOLUTION <- 600

################################################################################
# 2. HELPERS
################################################################################

clean_text <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x %in% c("", "NA", "N/A", "<NA>", "NaN")] <- NA_character_
  x
}

clean_key <- function(x) {
  toupper(clean_text(x))
}

clean_gene <- function(x) {
  x <- clean_text(x)
  x[tolower(x) == "no-gene"] <- NA_character_
  toupper(x)
}

write_source <- function(x, filename) {
  utils::write.csv(
    x,
    file.path(SOURCE_DATA_DIR, filename),
    row.names = FALSE
  )
}

format_group <- function(x) {
  x <- sub("^GEP1q-", "", as.character(x))
  gsub(">=", "≥", x, fixed = TRUE)
}

group_order <- function(x) {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & x != ""]
  
  if (!length(x)) {
    stop("No nonmissing RNBX groups available.")
  }
  
  c(
    HEATMAP_GROUP_ORDER[HEATMAP_GROUP_ORDER %in% x],
    setdiff(sort(x), HEATMAP_GROUP_ORDER)
  )
}

group_colors <- function(levels, custom) {
  missing <- setdiff(levels, names(custom))
  
  if (length(missing)) {
    stop(
      "Missing custom colors for: ",
      paste(missing, collapse = ", ")
    )
  }
  
  stats::setNames(
    unname(custom[levels]),
    levels
  )
}

quantile_scale <- function(x, probs, low, mid, high, digits = 3L) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  
  if (!length(x)) {
    stop("No finite annotation values available.")
  }
  
  lo <- unname(
    quantile(
      x,
      probs[1],
      na.rm = TRUE,
      names = FALSE
    )
  )
  
  md <- median(
    x,
    na.rm = TRUE
  )
  
  hi <- unname(
    quantile(
      x,
      probs[2],
      na.rm = TRUE,
      names = FALSE
    )
  )
  
  if (lo >= hi) {
    d <- max(abs(lo) * 0.01, 1e-6)
    lo <- lo - d
    hi <- hi + d
  }
  
  if (md <= lo || md >= hi) {
    md <- (lo + hi) / 2
  }
  
  at <- c(lo, md, hi)
  
  list(
    color_function = circlize::colorRamp2(
      at,
      c(low, mid, high)
    ),
    legend_at = at,
    legend_labels = formatC(
      at,
      format = "fg",
      digits = digits
    )
  )
}

transform_rows <- function(x, method) {
  means <- rowMeans(x)
  centered <- sweep(x, 1, means, "-")
  
  if (method == "Center") {
    return(
      list(
        matrix = centered,
        row_sds = rep(NA_real_, nrow(x)),
        zero_variance = rep(FALSE, nrow(x))
      )
    )
  }
  
  sds <- apply(x, 1, sd)
  zero <- !is.finite(sds) | sds <= .Machine$double.eps
  scaled <- centered
  
  if (any(!zero)) {
    scaled[!zero, ] <- sweep(
      centered[!zero, , drop = FALSE],
      1,
      sds[!zero],
      "/"
    )
  }
  
  if (any(zero)) {
    scaled[zero, ] <- 0
  }
  
  list(
    matrix = scaled,
    row_sds = sds,
    zero_variance = zero
  )
}

read_probe_ids <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " probe file not found:\n", path)
  }
  
  x <- readLines(
    path,
    warn = FALSE,
    encoding = "UTF-8"
  )
  
  if (length(x)) {
    x[1] <- sub("^\\ufeff", "", x[1])
  }
  
  x <- trimws(
    gsub(
      '"',
      "",
      x,
      fixed = TRUE
    )
  )
  
  x <- sub("[\\t,].*$", "", x)
  x <- clean_text(x)
  x <- x[!is.na(x)]
  
  x <- x[
    !tolower(x) %in% c(
      "probe",
      "probe_id",
      "probeid",
      "probeset",
      "probeset_id",
      "probesetid",
      "selected_probe_id",
      "selected_probe_ids"
    )
  ]
  
  x <- unique(x)
  
  if (!length(x)) {
    stop(label, " probe file contained no usable Probe_IDs.")
  }
  
  x
}

build_gene_matrix <- function(
    probe_ids,
    set_label,
    expr,
    fdat,
    sample_ids
) {
  expr_key <- clean_key(
    rownames(expr)
  )
  
  input_key <- clean_key(
    probe_ids
  )
  
  idx <- match(
    input_key,
    expr_key
  )
  
  resolution <- rep(
    "Exact",
    length(probe_ids)
  )
  
  unresolved <- is.na(idx)
  
  idx_t <- match(
    paste0(input_key, "T"),
    expr_key
  )
  
  repaired <- unresolved &
    !is.na(idx_t)
  
  idx[repaired] <- idx_t[repaired]
  
  resolution[repaired] <-
    "Appended terminal t"
  
  resolution[
    is.na(idx)
  ] <- "Unresolved"
  
  found <- !is.na(idx)
  
  matched_probe <- rep(
    NA_character_,
    length(probe_ids)
  )
  
  matched_gene <- rep(
    NA_character_,
    length(probe_ids)
  )
  
  matched_probe[found] <-
    rownames(expr)[idx[found]]
  
  matched_gene[found] <-
    clean_gene(
      fdat[[GENE_COL]][idx[found]]
    )
  
  audit <- tibble(
    Gene_Set = set_label,
    Input_Probe_ID = probe_ids,
    Probe_ID_Resolution = resolution,
    Found_In_ESET = found,
    Probe_ID = matched_probe,
    Gene = matched_gene
  )
  
  if (!any(found)) {
    return(
      list(
        matrix = NULL,
        row_labels = character(0),
        audit = audit,
        retained = tibble()
      )
    )
  }
  
  found_idx <- idx[found]
  found_probe <- rownames(expr)[found_idx]
  found_gene <- clean_gene(
    fdat[[GENE_COL]][found_idx]
  )
  
  found_expr <- expr[
    found_probe,
    sample_ids,
    drop = FALSE
  ]
  
  keep <- !is.na(found_gene) &
    rowSums(!is.finite(found_expr)) == 0
  
  retained <- tibble(
    Gene_Set = set_label,
    Input_Probe_ID = probe_ids[found][keep],
    Probe_ID_Resolution = resolution[found][keep],
    Probe_ID = found_probe[keep],
    Gene = found_gene[keep]
  )
  
  if (!nrow(retained)) {
    return(
      list(
        matrix = NULL,
        row_labels = character(0),
        audit = audit,
        retained = retained
      )
    )
  }
  
  retained_expr <- found_expr[
    keep,
    ,
    drop = FALSE
  ]
  
  retained_gene <- found_gene[
    keep
  ]
  
  gene_order <- unique(
    retained_gene
  )
  
  gene_sum <- rowsum(
    retained_expr,
    group = retained_gene,
    reorder = FALSE
  )
  
  n_probe <- table(
    retained_gene
  )
  
  gene_expr <- sweep(
    gene_sum,
    1,
    as.numeric(
      n_probe[
        rownames(gene_sum)
      ]
    ),
    "/"
  )
  
  gene_expr <- gene_expr[
    gene_order,
    ,
    drop = FALSE
  ]
  
  rownames(gene_expr) <- paste0(
    set_label,
    "__",
    gene_order
  )
  
  list(
    matrix = gene_expr,
    row_labels = gene_order,
    audit = audit,
    retained = retained
  )
}

save_ht_png <- function(ht, file) {
  if (file.exists(file)) {
    unlink(file)
  }
  
  png(
    file,
    width = FIGURE_WIDTH_IN,
    height = FIGURE_HEIGHT_IN,
    units = "in",
    res = PNG_RESOLUTION
  )
  
  on.exit(
    dev.off(),
    add = TRUE
  )
  
  ComplexHeatmap::draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    merge_legends = FALSE
  )
}

save_ht_pdf <- function(ht, file) {
  if (file.exists(file)) {
    unlink(file)
  }
  
  pdf(
    file,
    width = FIGURE_WIDTH_IN,
    height = FIGURE_HEIGHT_IN,
    useDingbats = FALSE
  )
  
  on.exit(
    dev.off(),
    add = TRUE
  )
  
  ComplexHeatmap::draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    merge_legends = FALSE
  )
}

################################################################################
# 3. LOAD EXPRESSIONSET
################################################################################

if (!file.exists(ESET_PATH)) {
  stop("Cannot find ESET:\n", ESET_PATH)
}

eset <- readRDS(ESET_PATH)

if (!methods::is(eset, "ExpressionSet")) {
  stop("Input is not an ExpressionSet.")
}

expr <- as.matrix(
  Biobase::exprs(eset)
)

storage.mode(expr) <- "numeric"

pdat <- as.data.frame(
  Biobase::pData(eset),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

fdat <- as.data.frame(
  Biobase::fData(eset),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_pdat <- c(
  SAMPLE_ID_COL,
  PATID_COL,
  SAMPLE_GROUP_COL,
  RNBX_GROUP_COL,
  RNAS_GROUP_COL,
  MATCH_STATUS_COL,
  MATCH_TYPE_COL,
  PCBM_BX_COL,
  POLYIG_COL
)

missing <- setdiff(
  required_pdat,
  colnames(pdat)
)

if (length(missing)) {
  stop(
    "Missing required pData columns: ",
    paste(missing, collapse = ", ")
  )
}

required_fdat <- c(
  GENE_COL,
  PC_ME_OTHER_COL
)

missing <- setdiff(
  required_fdat,
  colnames(fdat)
)

if (length(missing)) {
  stop(
    "Missing required fData columns: ",
    paste(missing, collapse = ", ")
  )
}

################################################################################
# 4. ALIGN EXPRESSION, pData, fData
################################################################################

expr_key <- clean_key(
  colnames(expr)
)

pdat_key <- clean_key(
  pdat[[SAMPLE_ID_COL]]
)

if (
  anyNA(expr_key) ||
  anyNA(pdat_key) ||
  anyDuplicated(expr_key) ||
  anyDuplicated(pdat_key)
) {
  stop(
    "Expression and pData sample identifiers must be nonmissing and unique."
  )
}

keep <- expr_key %in%
  pdat_key

if (!any(keep)) {
  stop(
    "No expression samples matched pData."
  )
}

expr <- expr[
  ,
  keep,
  drop = FALSE
]

idx <- match(
  clean_key(
    colnames(expr)
  ),
  pdat_key
)

pdat <- pdat[
  idx,
  ,
  drop = FALSE
]

rownames(pdat) <- colnames(expr)

fdat <- fdat[
  match(
    rownames(expr),
    rownames(fdat)
  ),
  ,
  drop = FALSE
]

rownames(fdat) <- rownames(expr)

if (
  !identical(
    colnames(expr),
    rownames(pdat)
  )
) {
  stop(
    "Expression/pData alignment failed."
  )
}

################################################################################
# 5. SAMPLE METADATA
################################################################################

meta <- tibble(
  Sample_ID = colnames(expr),
  PATID = pdat[[PATID_COL]],
  sample_group = pdat[[SAMPLE_GROUP_COL]],
  RNBX_Group = as.character(pdat[[RNBX_GROUP_COL]]),
  RNAS_Group = as.character(pdat[[RNAS_GROUP_COL]]),
  Match_Status = as.character(pdat[[MATCH_STATUS_COL]]),
  Match_Type = as.character(pdat[[MATCH_TYPE_COL]]),
  PCBmBx_All = pdat[[PCBM_BX_COL]],
  Row_PolyIG = pdat[[POLYIG_COL]]
)

meta$Is_RNAS_NDMM <-
  meta$sample_group ==
  RNAS_NDMM

meta$Is_RNAS_MGUS <-
  meta$sample_group ==
  RNAS_MGUS

meta$Is_RNBX_NDMM <-
  meta$sample_group ==
  RNBX_NDMM

meta$Is_RNBX_MGUS <-
  meta$sample_group ==
  RNBX_MGUS

meta$Is_RNBX_Target <-
  meta$Is_RNBX_NDMM |
  meta$Is_RNBX_MGUS

meta$Is_Strict_Match <-
  meta$Is_RNBX_Target &
  meta$Match_Status ==
  MATCHED_STATUS_VALUE

meta$Disease_Type_Consistent <- dplyr::case_when(
  meta$Is_RNBX_NDMM &
    meta$Match_Type ==
    "NDMM" ~ TRUE,
  
  meta$Is_RNBX_MGUS &
    meta$Match_Type ==
    "MGUS" ~ TRUE,
  
  meta$Is_RNBX_Target ~ FALSE,
  
  TRUE ~ NA
)

bad_type <- meta |>
  dplyr::filter(
    Is_Strict_Match,
    !Disease_Type_Consistent
  )

if (nrow(bad_type)) {
  write_source(
    bad_type,
    "Fig2A_Lasso_BestProbe_Matched_RNBX_DiseaseType_Inconsistencies.csv"
  )
  
  stop(
    "Some matched RNBX rows have inconsistent sample_group and match type."
  )
}

################################################################################
# 6. RNAS SOURCE LOOKUP — PATID + DISEASE TYPE
################################################################################

source_lookup <- meta |>
  dplyr::filter(
    Is_RNAS_NDMM |
      Is_RNAS_MGUS
  ) |>
  dplyr::transmute(
    PATID,
    Match_Type = dplyr::if_else(
      Is_RNAS_NDMM,
      "NDMM",
      "MGUS"
    ),
    Source_Sample_ID = Sample_ID,
    Source_RNAS_Group = RNAS_Group,
    Paired_RNAS_PolyIG = Row_PolyIG
  )

duplicate_sources <- source_lookup |>
  dplyr::filter(
    !is.na(PATID)
  ) |>
  dplyr::group_by(
    PATID,
    Match_Type
  ) |>
  dplyr::summarise(
    n = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::filter(
    n > 1L
  )

if (nrow(duplicate_sources)) {
  write_source(
    duplicate_sources,
    "Fig2A_Lasso_BestProbe_Duplicate_RNAS_PATID_Type.csv"
  )
  
  stop(
    "Multiple RNAS rows were found for the same PATID + disease type."
  )
}

################################################################################
# 7. PAIR RNBX TO RNAS
################################################################################

rnbx <- meta |>
  dplyr::filter(
    Is_RNBX_Target
  ) |>
  dplyr::left_join(
    source_lookup,
    by = c(
      "PATID",
      "Match_Type"
    )
  ) |>
  dplyr::mutate(
    Source_RNAS_Found =
      !is.na(
        Source_Sample_ID
      ),
    
    Source_Group_Matches_Target =
      Source_RNAS_Found &
      !is.na(RNBX_Group) &
      !is.na(Source_RNAS_Group) &
      RNBX_Group ==
      Source_RNAS_Group,
    
    Pairing_Integrity =
      dplyr::case_when(
        Match_Status !=
          MATCHED_STATUS_VALUE ~
          "Not a Step 21B one-to-one match",
        
        !Source_RNAS_Found ~
          "Matched RNAS source not found",
        
        !Source_Group_Matches_Target ~
          "RNAS and RNBX group-label mismatch",
        
        TRUE ~
          "Complete"
      )
  )

write_source(
  rnbx,
  "Fig2A_Lasso_BestProbe_RNBX_to_RNAS_PATID_DiseaseType_Pairing_Audit.csv"
)

strict_rnbx <- rnbx |>
  dplyr::filter(
    Pairing_Integrity ==
      "Complete"
  )

if (!nrow(strict_rnbx)) {
  stop(
    "No RNBX rows passed Step 21B pairing-integrity checks."
  )
}

################################################################################
# 8. LIMMA SAMPLE SELECTION
################################################################################

strict_rnbx <- strict_rnbx |>
  dplyr::mutate(
    MGUS_Selected =
      RNBX_Group ==
      MGUS_GROUP_VALUE &
      is.finite(PCBmBx_All) &
      PCBmBx_All <=
      MGUS_PCBMBX_MAX,
    
    GEP70_Selected =
      RNBX_Group ==
      GEP70_HIGH_GROUP_VALUE &
      is.finite(PCBmBx_All) &
      PCBmBx_All >=
      GEP70_PCBMBX_MIN &
      PCBmBx_All <
      GEP70_PCBMBX_MAX_EXCLUSIVE,
    
    Analysis_Group =
      dplyr::case_when(
        MGUS_Selected ~
          "MGUS",
        
        GEP70_Selected ~
          "GEP70>0.66",
        
        TRUE ~
          NA_character_
      )
  )

selected <- strict_rnbx |>
  dplyr::filter(
    !is.na(
      Analysis_Group
    )
  ) |>
  dplyr::mutate(
    Analysis_Group = factor(
      Analysis_Group,
      levels = c(
        "MGUS",
        "GEP70>0.66"
      )
    )
  ) |>
  dplyr::arrange(
    Analysis_Group,
    Sample_ID
  )

n_mgus <- sum(
  selected$Analysis_Group ==
    "MGUS"
)

n_gep70 <- sum(
  selected$Analysis_Group ==
    "GEP70>0.66"
)

if (
  n_mgus < 2 ||
  n_gep70 < 2
) {
  stop(
    "Each limma group must contain at least two samples."
  )
}

################################################################################
# 9. PAIRED RNAS HEATMAP SAMPLE SELECTION
################################################################################

heatmap_pairs <- if (
  HEATMAP_INCLUDE_ALL_NONMISSING_RNBX_GROUPS
) {
  strict_rnbx |>
    dplyr::filter(
      !is.na(
        RNBX_Group
      )
    )
} else {
  selected
}

heatmap_metadata <- heatmap_pairs |>
  dplyr::transmute(
    RNBX_Sample_ID = Sample_ID,
    Sample_ID = Source_Sample_ID,
    PATID = PATID,
    sample_group = dplyr::if_else(
      Match_Type == "NDMM",
      RNAS_NDMM,
      RNAS_MGUS
    ),
    RNBX_Group = RNBX_Group,
    RNAS_Group = Source_RNAS_Group,
    Match_Type = Match_Type,
    Paired_RNAS_PolyIG = Paired_RNAS_PolyIG
  )

if (nrow(heatmap_metadata) < 2) {
  stop(
    "Fewer than two paired RNAS samples are available."
  )
}

if (
  anyNA(
    heatmap_metadata$Sample_ID
  )
) {
  stop(
    "A retained pair is missing its RNAS sample ID."
  )
}

if (
  anyDuplicated(
    heatmap_metadata$Sample_ID
  )
) {
  stop(
    "A paired RNAS sample was assigned to more than one retained pair."
  )
}

if (
  !all(
    heatmap_metadata$Sample_ID %in%
    colnames(expr)
  )
) {
  stop(
    "A paired RNAS sample is missing from the expression matrix."
  )
}

################################################################################
# 10. PROBE-LEVEL LIMMA ON RNBX
################################################################################

gene <- clean_gene(
  fdat[[GENE_COL]]
)

pc_me <- suppressWarnings(
  as.integer(
    as.character(
      fdat[[PC_ME_OTHER_COL]]
    )
  )
)

valid_probe <-
  !is.na(gene) &
  !is.na(pc_me) &
  pc_me %in% c(
    0L,
    1L,
    2L
  ) &
  rowSums(
    !is.finite(expr)
  ) == 0

probe_expr <- expr[
  valid_probe,
  ,
  drop = FALSE
]

probe_gene <- gene[
  valid_probe
]

probe_pcme <- pc_me[
  valid_probe
]

probe_expr_limma <- probe_expr[
  ,
  selected$Sample_ID,
  drop = FALSE
]

design <- model.matrix(
  ~ 0 +
    selected$Analysis_Group
)

colnames(design) <- c(
  "MGUS",
  "GEP70_High"
)

rownames(design) <-
  selected$Sample_ID

contrast <- limma::makeContrasts(
  GEP70_High_minus_MGUS =
    GEP70_High - MGUS,
  levels = design
)

fit <- limma::lmFit(
  probe_expr_limma,
  design
)

fit <- limma::contrasts.fit(
  fit,
  contrast
)

fit <- limma::eBayes(
  fit
)

tt <- limma::topTable(
  fit,
  coef = "GEP70_High_minus_MGUS",
  number = Inf,
  sort.by = "none",
  p.value = 1,
  lfc = 0
)

mean_mgus <- rowMeans(
  probe_expr_limma[
    ,
    selected$Analysis_Group ==
      "MGUS",
    drop = FALSE
  ]
)

mean_gep70 <- rowMeans(
  probe_expr_limma[
    ,
    selected$Analysis_Group ==
      "GEP70>0.66",
    drop = FALSE
  ]
)

probe_results <- tibble(
  Probe_ID = rownames(probe_expr_limma),
  Gene = probe_gene,
  PC_ME_Other = probe_pcme,
  N_MGUS = n_mgus,
  N_GEP70_High = n_gep70,
  Mean_MGUS = mean_mgus,
  Mean_GEP70_High = mean_gep70,
  Log2FC_GEP70_High_vs_MGUS = tt$logFC,
  PValue = tt$P.Value,
  FDR = tt$adj.P.Val
)

write_source(
  probe_results,
  "Fig2A_Lasso_BestProbe_MGUS_vs_GEP70_LIMMA_AllProbeResults.csv"
)

################################################################################
# 11. ONE REPRESENTATIVE PROBE PER GENE
################################################################################

best_probe <- probe_results |>
  dplyr::filter(
    !is.na(Gene),
    is.finite(FDR),
    is.finite(PValue),
    is.finite(
      Log2FC_GEP70_High_vs_MGUS
    )
  ) |>
  dplyr::group_by(
    Gene
  ) |>
  dplyr::arrange(
    FDR,
    PValue,
    dplyr::desc(
      abs(
        Log2FC_GEP70_High_vs_MGUS
      )
    ),
    Probe_ID,
    .by_group = TRUE
  ) |>
  dplyr::slice_head(
    n = 1L
  ) |>
  dplyr::ungroup()

if (!nrow(best_probe)) {
  stop(
    "No representative probes remained."
  )
}

best_idx <- match(
  best_probe$Probe_ID,
  rownames(probe_expr)
)

best_gene_expr <- probe_expr[
  best_idx,
  ,
  drop = FALSE
]

rownames(best_gene_expr) <-
  best_probe$Gene

################################################################################
# 12. STEP 23C MICROENVIRONMENT GENE SETS
################################################################################

de <- best_probe |>
  dplyr::mutate(
    ME_GeneSet =
      dplyr::case_when(
        PC_ME_Other ==
          2L &
          Log2FC_GEP70_High_vs_MGUS <
          0 &
          FDR <
          MGUS_ME_FDR_CUTOFF ~
          "MGUS_ME",
        
        PC_ME_Other ==
          2L &
          Log2FC_GEP70_High_vs_MGUS >
          0 &
          FDR <
          GEP70_ME_FDR_CUTOFF ~
          "GEP70_ME",
        
        TRUE ~
          "Not_selected"
      )
  )

mgus_me <- de |>
  dplyr::filter(
    ME_GeneSet ==
      "MGUS_ME"
  )

gep70_me <- de |>
  dplyr::filter(
    ME_GeneSet ==
      "GEP70_ME"
  )

if (
  !nrow(mgus_me) ||
  !nrow(gep70_me)
) {
  stop(
    "Both MGUS_ME and GEP70_ME genes are required for the supervised score."
  )
}

write_source(
  best_probe,
  "Fig2A_Lasso_BestProbe_MGUS_vs_GEP70_LIMMA_BestProbePerGene.csv"
)

write_source(
  mgus_me,
  "Fig2A_Lasso_BestProbe_MGUS_ME_Selected.csv"
)

write_source(
  gep70_me,
  "Fig2A_Lasso_BestProbe_GEP70_ME_Selected.csv"
)

################################################################################
# 13. SUPERVISED ME SCORE FROM PAIRED RNBX SAMPLES
################################################################################

high_genes <- intersect(
  gep70_me$Gene,
  rownames(best_gene_expr)
)

low_genes <- intersect(
  mgus_me$Gene,
  rownames(best_gene_expr)
)

high_mean <- colMeans(
  best_gene_expr[
    high_genes,
    strict_rnbx$Sample_ID,
    drop = FALSE
  ]
)

low_mean <- colMeans(
  best_gene_expr[
    low_genes,
    strict_rnbx$Sample_ID,
    drop = FALSE
  ]
)

rnbx_score <- tibble(
  RNBX_Sample_ID =
    strict_rnbx$Sample_ID,
  
  Supervised_ME_HighMinusLow =
    unname(
      high_mean[
        strict_rnbx$Sample_ID
      ]
    ) -
    unname(
      low_mean[
        strict_rnbx$Sample_ID
      ]
    )
)

heatmap_metadata <- heatmap_metadata |>
  dplyr::left_join(
    rnbx_score,
    by = "RNBX_Sample_ID"
  )

if (
  any(
    !is.finite(
      heatmap_metadata$Supervised_ME_HighMinusLow
    )
  )
) {
  stop(
    "A heatmap sample is missing its supervised ME score."
  )
}

################################################################################
# 14. READ GEP1Q AND POLYIG PROBE LISTS
################################################################################

gep1q_probe_ids <- read_probe_ids(
  GEP1Q_PROBE_FILE,
  "GEP1q"
)

polyig_probe_ids <- read_probe_ids(
  POLYIG_PROBE_FILE,
  "PolyIG"
)

################################################################################
# 15. BUILD PAIRED-RNAS GEP1Q AND POLYIG MATRICES
################################################################################

gep1q_set <- build_gene_matrix(
  probe_ids = gep1q_probe_ids,
  set_label = "GEP1q",
  expr = expr,
  fdat = fdat,
  sample_ids = heatmap_metadata$Sample_ID
)

polyig_set <- build_gene_matrix(
  probe_ids = polyig_probe_ids,
  set_label = "PolyIG",
  expr = expr,
  fdat = fdat,
  sample_ids = heatmap_metadata$Sample_ID
)

write_source(
  gep1q_set$audit,
  "Fig2A_Lasso_BestProbe_GEP1q_Selected_Probe_Audit.csv"
)

write_source(
  polyig_set$audit,
  "Fig2A_Lasso_BestProbe_PolyIG_Selected_Probe_Audit.csv"
)

if (is.null(gep1q_set$matrix)) {
  stop(
    "No GEP1q probes remained."
  )
}

if (is.null(polyig_set$matrix)) {
  stop(
    "No PolyIG probes remained."
  )
}

################################################################################
# 16. VERIFY FROZEN GEP1Q FEATURES
################################################################################

gep1q_unresolved <- gep1q_set$audit |>
  dplyr::filter(
    !Found_In_ESET
  )

if (nrow(gep1q_unresolved)) {
  write_source(
    gep1q_unresolved,
    "Fig2A_Lasso_BestProbe_GEP1q_Selected_Probe_Unresolved.csv"
  )
  
  stop(
    "Some frozen GEP1q Probe_IDs could not be matched to the ExpressionSet."
  )
}

if (
  nrow(
    gep1q_set$retained
  ) !=
  length(
    gep1q_probe_ids
  )
) {
  stop(
    "Not every frozen GEP1q Probe_ID was retained. Input: ",
    length(gep1q_probe_ids),
    "; retained: ",
    nrow(gep1q_set$retained)
  )
}

gep1q_gene_counts <- gep1q_set$retained |>
  dplyr::group_by(
    Gene
  ) |>
  dplyr::summarise(
    n = dplyr::n(),
    .groups = "drop"
  )

if (
  any(
    gep1q_gene_counts$n != 1L
  )
) {
  write_source(
    gep1q_gene_counts,
    "Fig2A_Lasso_BestProbe_GEP1q_Selected_Probe_Counts_By_Gene.csv"
  )
  
  stop(
    "The frozen GEP1q probe file is not one-probe-per-gene."
  )
}

################################################################################
# 17. COMBINE POLYIG + GEP1Q HEATMAP ROWS
################################################################################

# Preserve original display order:
#   1. PolyIG
#   2. GEP1q

heatmap_expression <- rbind(
  polyig_set$matrix,
  gep1q_set$matrix
)

heatmap_row_labels <- c(
  polyig_set$row_labels,
  gep1q_set$row_labels
)

row_gene_set <- factor(
  c(
    rep(
      "PolyIG",
      nrow(
        polyig_set$matrix
      )
    ),
    rep(
      "GEP1q",
      nrow(
        gep1q_set$matrix
      )
    )
  ),
  levels = c(
    "PolyIG",
    "GEP1q"
  )
)

if (
  !identical(
    colnames(
      heatmap_expression
    ),
    heatmap_metadata$Sample_ID
  )
) {
  stop(
    "Heatmap expression columns and metadata are not aligned."
  )
}

################################################################################
# 18. GROUP ORDER + SUPERVISED COLUMN ORDER
################################################################################

heatmap_group_levels <- group_order(
  heatmap_metadata$RNBX_Group
)

heatmap_metadata$RNBX_Group <- factor(
  heatmap_metadata$RNBX_Group,
  levels =
    heatmap_group_levels
)

if (
  ORDER_COLUMNS_BY_SUPERVISED_ME_SCORE
) {
  score_order <-
    heatmap_metadata$Supervised_ME_HighMinusLow
  
  if (
    SUPERVISED_ME_SCORE_LEFT_TO_RIGHT ==
    "HighToLow"
  ) {
    score_order <- -score_order
  }
  
  idx <- if (
    SPLIT_SUPERVISED_ORDER_BY_RNBX_GROUP
  ) {
    order(
      heatmap_metadata$RNBX_Group,
      score_order,
      heatmap_metadata$Sample_ID
    )
  } else {
    order(
      score_order,
      heatmap_metadata$Sample_ID
    )
  }
  
  heatmap_metadata <- heatmap_metadata[
    idx,
    ,
    drop = FALSE
  ]
  
  heatmap_expression <- heatmap_expression[
    ,
    heatmap_metadata$Sample_ID,
    drop = FALSE
  ]
}

group_counts <- table(
  heatmap_metadata$RNBX_Group
)

split_labels <- stats::setNames(
  paste0(
    format_group(
      names(
        group_counts
      )
    ),
    "\n(n = ",
    as.integer(
      group_counts
    ),
    ")"
  ),
  names(
    group_counts
  )
)

heatmap_metadata$Column_Split <- factor(
  split_labels[
    as.character(
      heatmap_metadata$RNBX_Group
    )
  ],
  levels =
    unname(
      split_labels[
        heatmap_group_levels
      ]
    )
)

write_source(
  heatmap_metadata,
  "Fig2A_Lasso_BestProbe_Heatmap_Retained_RNAS_SampleMetadata_ColumnOrder.csv"
)

################################################################################
# 19. TRANSFORM HEATMAP ROWS
################################################################################

heatmap_transform <- transform_rows(
  heatmap_expression,
  HEATMAP_SCALE_METHOD
)

scaled_heatmap_expression <-
  heatmap_transform$matrix

if (
  any(
    !is.finite(
      scaled_heatmap_expression
    )
  )
) {
  stop(
    "Transformed heatmap matrix contains nonfinite values."
  )
}

################################################################################
# 20. COLORS
################################################################################

heatmap_color_limit <- if (
  HEATMAP_SCALE_METHOD ==
  "Zscore"
) {
  Z_SCORE_COLOR_LIMIT
} else {
  max(
    abs(
      scaled_heatmap_expression
    ),
    na.rm = TRUE
  )
}

if (
  !is.finite(
    heatmap_color_limit
  ) ||
  heatmap_color_limit <= 0
) {
  heatmap_color_limit <- 1
}

heatmap_color_function <- circlize::colorRamp2(
  c(
    -heatmap_color_limit,
    0,
    heatmap_color_limit
  ),
  c(
    "#2166AC",
    "white",
    "#B2182B"
  )
)

heatmap_legend_at <- seq(
  -heatmap_color_limit,
  heatmap_color_limit,
  length.out = 5
)

row_gene_set_colors <- c(
  "PolyIG" = "steelblue3",
  "GEP1q" = "red3"
)

group_cols <- group_colors(
  heatmap_group_levels,
  CUSTOM_RNBX_GROUP_COLORS
)

polyig_scale <- quantile_scale(
  heatmap_metadata$Paired_RNAS_PolyIG,
  POLYIG_COLOR_QUANTILES,
  "green",
  "black",
  "red",
  4L
)

score_limit <- max(
  abs(
    heatmap_metadata$Supervised_ME_HighMinusLow
  ),
  na.rm = TRUE
)

if (
  !is.finite(
    score_limit
  ) ||
  score_limit <= 0
) {
  score_limit <- 1
}

score_color_function <- circlize::colorRamp2(
  c(
    -score_limit,
    0,
    score_limit
  ),
  c(
    "#2166AC",
    "white",
    "#B2182B"
  )
)

score_legend_at <- seq(
  -score_limit,
  score_limit,
  length.out = 5
)

################################################################################
# 21. TOP ANNOTATION
################################################################################

top_annotation <- ComplexHeatmap::HeatmapAnnotation(
  `GEP Groups` =
    heatmap_metadata$RNBX_Group,
  
  PolyIG =
    heatmap_metadata$Paired_RNAS_PolyIG,
  
  ME_HighMinusLow =
    heatmap_metadata$Supervised_ME_HighMinusLow,
  
  col = list(
    `GEP Groups` =
      group_cols,
    
    PolyIG =
      polyig_scale$color_function,
    
    ME_HighMinusLow =
      score_color_function
  ),
  
  annotation_name_gp =
    grid::gpar(
      fontsize = 10,
      fontface = "bold"
    ),
  
  simple_anno_size =
    grid::unit(
      4,
      "mm"
    ),
  
  na_col =
    "grey85",
  
  annotation_legend_param =
    list(
      `GEP Groups` =
        list(
          title =
            "GEP Groups",
          
          at =
            names(
              group_cols
            ),
          
          labels =
            format_group(
              names(
                group_cols
              )
            )
        ),
      
      PolyIG =
        list(
          title =
            "Paired RNAS\nPolyIG",
          
          at =
            polyig_scale$legend_at,
          
          labels =
            polyig_scale$legend_labels
        ),
      
      ME_HighMinusLow =
        list(
          title =
            "Supervised ME score\nmean high - mean low",
          
          at =
            score_legend_at,
          
          labels =
            formatC(
              score_legend_at,
              format = "fg",
              digits = 3
            )
        )
    )
)

################################################################################
# 22. ROW ANNOTATION
################################################################################

row_annotation <- ComplexHeatmap::rowAnnotation(
  Signature =
    row_gene_set,
  
  col = list(
    Signature =
      row_gene_set_colors
  ),
  
  annotation_name_gp =
    grid::gpar(
      fontsize = 8,
      fontface = "bold"
    ),
  
  simple_anno_size =
    grid::unit(
      4,
      "mm"
    ),
  
  annotation_legend_param =
    list(
      Signature =
        list(
          title =
            "Selected gene set"
        )
    )
)

################################################################################
# 23. CREATE FIGURE 2A
################################################################################

row_name_font_size <- if (
  nrow(
    scaled_heatmap_expression
  ) <= 50
) {
  8
} else {
  6
}

heatmap_title <-
  "GEP 1q Copy Number (CN) and Polyclonal Immunoglobulin Expression Across MGUS and NDMM Risk Groups"

figure_2a_heatmap <- ComplexHeatmap::Heatmap(
  scaled_heatmap_expression,
  
  name = if (
    HEATMAP_SCALE_METHOD ==
    "Zscore"
  ) {
    "Row z-score\nof log2 expression"
  } else {
    "Row mean-centered\nlog2 expression"
  },
  
  col =
    heatmap_color_function,
  
  top_annotation =
    top_annotation,
  
  left_annotation =
    row_annotation,
  
  row_split =
    row_gene_set,
  
  column_split =
    heatmap_metadata$Column_Split,
  
  cluster_rows =
    TRUE,
  
  cluster_row_slices =
    FALSE,
  
  cluster_columns =
    FALSE,
  
  cluster_column_slices =
    FALSE,
  
  show_row_dend =
    TRUE,
  
  show_column_dend =
    FALSE,
  
  show_row_names =
    SHOW_ROW_NAMES,
  
  show_column_names =
    SHOW_COLUMN_NAMES,
  
  row_labels =
    heatmap_row_labels,
  
  row_names_gp =
    grid::gpar(
      fontsize =
        row_name_font_size
    ),
  
  row_title_gp =
    grid::gpar(
      fontsize = 12,
      fontface = "bold"
    ),
  
  column_title =
    heatmap_title,
  
  column_title_gp =
    grid::gpar(
      fontsize = 20,
      fontface = "bold"
    ),
  
  column_gap =
    grid::unit(
      2.5,
      "mm"
    ),
  
  row_gap =
    grid::unit(
      2.5,
      "mm"
    ),
  
  use_raster =
    ncol(
      scaled_heatmap_expression
    ) > 200,
  
  raster_device =
    "png",
  
  heatmap_legend_param =
    list(
      at =
        heatmap_legend_at,
      
      labels =
        formatC(
          heatmap_legend_at,
          format = "fg",
          digits = 3
        ),
      
      title_gp =
        grid::gpar(
          fontface = "bold"
        )
    )
)

ComplexHeatmap::draw(
  figure_2a_heatmap,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legends = FALSE
)

################################################################################
# 24. SAVE FIGURE
################################################################################

tag <- if (
  HEATMAP_SCALE_METHOD ==
  "Zscore"
) {
  "Zscore"
} else {
  "MeanCentered"
}

png_file <- file.path(
  OUT_DIR,
  paste0(
    FIGURE_ID,
    "_RNAS_GEP1q_PolyIG_Heatmap_",
    tag,
    "_PATID_DiseaseType_PairedRNAS.png"
  )
)

pdf_file <- file.path(
  OUT_DIR,
  paste0(
    FIGURE_ID,
    "_RNAS_GEP1q_PolyIG_Heatmap_",
    tag,
    "_PATID_DiseaseType_PairedRNAS.pdf"
  )
)

save_ht_png(
  figure_2a_heatmap,
  png_file
)

save_ht_pdf(
  figure_2a_heatmap,
  pdf_file
)

################################################################################
# 25. SOURCE DATA
################################################################################

write_source(
  gep1q_set$retained,
  "Fig2A_Lasso_BestProbe_GEP1q_Retained_Selected_Probes.csv"
)

write_source(
  polyig_set$retained,
  "Fig2A_Lasso_BestProbe_PolyIG_Retained_Selected_Probes.csv"
)

write_source(
  tibble(
    Internal_Row_ID =
      rownames(
        scaled_heatmap_expression
      ),
    
    GENE =
      heatmap_row_labels,
    
    Signature =
      as.character(
        row_gene_set
      )
  ),
  "Fig2A_Lasso_BestProbe_Heatmap_RowAnnotations.csv"
)

################################################################################
# 26. SUMMARY
################################################################################

message(
  "\nStep 24C complete",
  "\nRNAS/RNBX pairing: PATID + disease type",
  "\nMGUS limma samples: ", n_mgus,
  "\nGEP70-high limma samples: ", n_gep70,
  "\nMGUS_ME genes: ", nrow(mgus_me),
  "\nGEP70_ME genes: ", nrow(gep70_me),
  "\nPaired RNAS heatmap samples: ", nrow(heatmap_metadata),
  "\nPolyIG heatmap rows: ", nrow(polyig_set$matrix),
  "\nGEP1q heatmap rows: ", nrow(gep1q_set$matrix),
  "\nHeatmap row order: PolyIG followed by GEP1q",
  "\nTop annotations: GEP Groups, PolyIG, Supervised ME score",
  "\nPNG: ", png_file,
  "\nPDF: ", pdf_file,
  "\n"
)