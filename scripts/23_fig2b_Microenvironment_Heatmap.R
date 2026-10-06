################################################################################
# STEP 23C — FIGURE 2B
# BEST PROBE PER GENE MICROENVIRONMENT HEATMAP
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

FIGURE_ID <- "Fig2B_Lasso_BestProbe"

ESET_PATH <- file.path(
  "data", "processed", "uams",
  "ESET_uams_all_gep1q_fish_polyig_pcme_binaries_lasso_groups_rnbx-matched.rds"
)

FIGURE_DIR <- file.path("results", "figures", "Fig2", FIGURE_ID)
SOURCE_DATA_DIR <- file.path(FIGURE_DIR, "source-data")
OUT_DIR <- file.path(FIGURE_DIR, paste0(FIGURE_ID, "_Microenvironment_Heatmap"))

MATCH_STATUS_COLUMN <- "RNAS_to_RNBX_PATID_match_MGUS_Lasso"
MATCH_TYPE_COLUMN <- "RNAS_to_RNBX_match_type_Lasso"
MATCHED_STATUS_VALUE <- "matched_one_to_one"

RNAS_NDMM_SAMPLE_GROUP <- "RNAS_CD138_NDMM"
RNBX_NDMM_SAMPLE_GROUP <- "RNBX_WB_NDMM"
RNAS_MGUS_SAMPLE_GROUP <- "RNAS_CD138_MGUS"
RNBX_MGUS_SAMPLE_GROUP <- "RNBX_WB_MGUS"

SAMPLE_ID_COL <- "chipid"
PATID_COL <- "patid"
SAMPLE_GROUP_COL <- "sample_group"

RNBX_GROUP_COLUMN <- "GEP_RNBX_Groups_MGUS_Lasso"
RNAS_GROUP_COLUMN <- "GEP_RNAS_Groups_MGUS_Lasso"

PCBM_BX_COL <- "PCBmBx_All"
POLYIG_COL <- "PolyIG_Score"

GENE_COL <- "GENE"
PC_ME_OTHER_COLUMN <- "PC_ME_Other"

MGUS_GROUP_VALUE <- "MGUS"
GEP70_HIGH_GROUP_VALUE <- "GEP70>0.66"

MGUS_PCBMBX_MAX <- 100
GEP70_PCBMBX_MIN <- 0
GEP70_PCBMBX_MAX_EXCLUSIVE <- 70

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

dir.create(SOURCE_DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

clean_text <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x %in% c("", "NA", "N/A", "<NA>", "NaN")] <- NA_character_
  x
}

clean_key <- function(x) toupper(clean_text(x))

clean_gene <- function(x) {
  x <- clean_text(x)
  x[tolower(x) == "no-gene"] <- NA_character_
  toupper(x)
}

write_source <- function(x, filename) {
  utils::write.csv(x, file.path(SOURCE_DATA_DIR, filename), row.names = FALSE)
}

format_group <- function(x) {
  x <- sub("^GEP1q-", "", as.character(x))
  gsub(">=", "≥", x, fixed = TRUE)
}

group_order <- function(x) {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & x != ""]
  if (!length(x)) stop("No nonmissing RNBX groups available.")
  
  c(
    HEATMAP_GROUP_ORDER[HEATMAP_GROUP_ORDER %in% x],
    setdiff(sort(x), HEATMAP_GROUP_ORDER)
  )
}

group_colors <- function(levels, custom) {
  missing <- setdiff(levels, names(custom))
  if (length(missing)) {
    stop("Missing custom colors for: ", paste(missing, collapse = ", "))
  }
  unname(custom[levels]) |> stats::setNames(levels)
}

quantile_scale <- function(x, probs, low, mid, high, digits = 3L) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) stop("No finite annotation values available.")
  
  lo <- unname(quantile(x, probs[1], na.rm = TRUE, names = FALSE))
  md <- median(x, na.rm = TRUE)
  hi <- unname(quantile(x, probs[2], na.rm = TRUE, names = FALSE))
  
  if (lo >= hi) {
    d <- max(abs(lo) * 0.01, 1e-6)
    lo <- lo - d
    hi <- hi + d
  }
  
  if (md <= lo || md >= hi) md <- (lo + hi) / 2
  at <- c(lo, md, hi)
  
  list(
    color_function = circlize::colorRamp2(at, c(low, mid, high)),
    legend_at = at,
    legend_labels = formatC(at, format = "fg", digits = digits)
  )
}

transform_rows <- function(x, method) {
  means <- rowMeans(x)
  centered <- sweep(x, 1, means, "-")
  
  if (method == "Center") {
    return(list(
      matrix = centered,
      row_sds = rep(NA_real_, nrow(x)),
      zero_variance = rep(FALSE, nrow(x))
    ))
  }
  
  sds <- apply(x, 1, sd)
  zero <- !is.finite(sds) | sds <= .Machine$double.eps
  scaled <- centered
  
  if (any(!zero)) {
    scaled[!zero, ] <- sweep(centered[!zero, , drop = FALSE], 1, sds[!zero], "/")
  }
  
  if (any(zero)) scaled[zero, ] <- 0
  
  list(
    matrix = scaled,
    row_sds = sds,
    zero_variance = zero
  )
}

save_ht_png <- function(ht, file) {
  if (file.exists(file)) unlink(file)
  png(file, width = FIGURE_WIDTH_IN, height = FIGURE_HEIGHT_IN,
      units = "in", res = PNG_RESOLUTION)
  on.exit(dev.off(), add = TRUE)
  ComplexHeatmap::draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    merge_legends = FALSE
  )
}

save_ht_pdf <- function(ht, file) {
  if (file.exists(file)) unlink(file)
  pdf(file, width = FIGURE_WIDTH_IN, height = FIGURE_HEIGHT_IN, useDingbats = FALSE)
  on.exit(dev.off(), add = TRUE)
  ComplexHeatmap::draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    merge_legends = FALSE
  )
}

################################################################################
# 3. LOAD AND ALIGN EXPRESSIONSET
################################################################################

if (!file.exists(ESET_PATH)) stop("Cannot find ESET:\n", ESET_PATH)

eset <- readRDS(ESET_PATH)
if (!methods::is(eset, "ExpressionSet")) stop("Input is not an ExpressionSet.")

expr <- as.matrix(Biobase::exprs(eset))
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
  RNBX_GROUP_COLUMN,
  RNAS_GROUP_COLUMN,
  MATCH_STATUS_COLUMN,
  MATCH_TYPE_COLUMN,
  PCBM_BX_COL,
  POLYIG_COL
)

missing_pdat <- setdiff(required_pdat, colnames(pdat))
if (length(missing_pdat)) {
  stop("Missing required pData columns: ", paste(missing_pdat, collapse = ", "))
}

required_fdat <- c(GENE_COL, PC_ME_OTHER_COLUMN)
missing_fdat <- setdiff(required_fdat, colnames(fdat))

if (length(missing_fdat)) {
  stop("Missing required fData columns: ", paste(missing_fdat, collapse = ", "))
}

expr_key <- clean_key(colnames(expr))
pdat_key <- clean_key(pdat[[SAMPLE_ID_COL]])

if (anyNA(expr_key) || anyNA(pdat_key) ||
    anyDuplicated(expr_key) || anyDuplicated(pdat_key)) {
  stop("Expression and pData sample identifiers must be nonmissing and unique.")
}

keep <- expr_key %in% pdat_key
if (!any(keep)) stop("No expression samples matched pData.")

expr <- expr[, keep, drop = FALSE]

idx <- match(clean_key(colnames(expr)), pdat_key)
pdat <- pdat[idx, , drop = FALSE]
rownames(pdat) <- colnames(expr)

fdat <- fdat[match(rownames(expr), rownames(fdat)), , drop = FALSE]
rownames(fdat) <- rownames(expr)

if (!identical(colnames(expr), rownames(pdat))) {
  stop("Expression/pData alignment failed.")
}

################################################################################
# 4. BUILD RNAS -> RNBX PAIRING BY PATID + DISEASE TYPE
################################################################################

meta <- tibble::tibble(
  Sample_ID = colnames(expr),
  PATID = pdat[[PATID_COL]],
  sample_group = pdat[[SAMPLE_GROUP_COL]],
  RNBX_Group = as.character(pdat[[RNBX_GROUP_COLUMN]]),
  RNAS_Group = as.character(pdat[[RNAS_GROUP_COLUMN]]),
  Match_Status = as.character(pdat[[MATCH_STATUS_COLUMN]]),
  Match_Type = as.character(pdat[[MATCH_TYPE_COLUMN]]),
  PCBmBx_All = pdat[[PCBM_BX_COL]],
  Row_PolyIG = pdat[[POLYIG_COL]]
)

meta$Is_RNAS_NDMM <- meta$sample_group == RNAS_NDMM_SAMPLE_GROUP
meta$Is_RNAS_MGUS <- meta$sample_group == RNAS_MGUS_SAMPLE_GROUP
meta$Is_RNBX_NDMM <- meta$sample_group == RNBX_NDMM_SAMPLE_GROUP
meta$Is_RNBX_MGUS <- meta$sample_group == RNBX_MGUS_SAMPLE_GROUP
meta$Is_RNBX_Target <- meta$Is_RNBX_NDMM | meta$Is_RNBX_MGUS

meta$Is_Strict_Match <- meta$Is_RNBX_Target &
  meta$Match_Status == MATCHED_STATUS_VALUE

meta$Disease_Type_Consistent <- dplyr::case_when(
  meta$Is_RNBX_NDMM & meta$Match_Type == "NDMM" ~ TRUE,
  meta$Is_RNBX_MGUS & meta$Match_Type == "MGUS" ~ TRUE,
  meta$Is_RNBX_Target ~ FALSE,
  TRUE ~ NA
)

bad_type <- meta |>
  dplyr::filter(Is_Strict_Match, !Disease_Type_Consistent)

if (nrow(bad_type)) {
  write_source(
    bad_type,
    "Fig2B_Lasso_BestProbe_Matched_RNBX_DiseaseType_Inconsistencies.csv"
  )
  stop("Some matched RNBX rows have inconsistent sample_group and match type.")
}

source_lookup <- meta |>
  dplyr::filter(Is_RNAS_NDMM | Is_RNAS_MGUS) |>
  dplyr::transmute(
    PATID,
    Match_Type = dplyr::if_else(Is_RNAS_NDMM, "NDMM", "MGUS"),
    Source_RNAS_Group = RNAS_Group,
    Paired_RNAS_PolyIG = Row_PolyIG
  )

duplicate_sources <- source_lookup |>
  dplyr::filter(!is.na(PATID)) |>
  dplyr::group_by(PATID, Match_Type) |>
  dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
  dplyr::filter(n > 1L)

if (nrow(duplicate_sources)) {
  write_source(
    duplicate_sources,
    "Fig2B_Lasso_BestProbe_Duplicate_RNAS_PATID_Type.csv"
  )
  stop("Multiple RNAS rows were found for the same PATID + disease type.")
}

rnbx <- meta |>
  dplyr::filter(Is_RNBX_Target) |>
  dplyr::left_join(
    source_lookup,
    by = c("PATID", "Match_Type")
  ) |>
  dplyr::mutate(
    Source_RNAS_Found = !is.na(Source_RNAS_Group),
    Source_Group_Matches_Target =
      Source_RNAS_Found &
      !is.na(RNBX_Group) &
      RNBX_Group == Source_RNAS_Group,
    
    Pairing_Integrity = dplyr::case_when(
      Match_Status != MATCHED_STATUS_VALUE ~
        "Not a Step 21B one-to-one match",
      !Source_RNAS_Found ~
        "Matched RNAS source not found",
      !Source_Group_Matches_Target ~
        "RNAS and RNBX group-label mismatch",
      TRUE ~
        "Complete"
    )
  )

strict_rnbx <- rnbx |>
  dplyr::filter(Pairing_Integrity == "Complete")

if (!nrow(strict_rnbx)) {
  stop("No RNBX rows passed Step 21B pairing-integrity checks.")
}

################################################################################
# 5. LIMMA SAMPLE SELECTION
################################################################################

strict_rnbx <- strict_rnbx |>
  dplyr::mutate(
    MGUS_Selected =
      RNBX_Group == MGUS_GROUP_VALUE &
      is.finite(PCBmBx_All) &
      PCBmBx_All <= MGUS_PCBMBX_MAX,
    
    GEP70_Selected =
      RNBX_Group == GEP70_HIGH_GROUP_VALUE &
      is.finite(PCBmBx_All) &
      PCBmBx_All >= GEP70_PCBMBX_MIN &
      PCBmBx_All < GEP70_PCBMBX_MAX_EXCLUSIVE,
    
    Analysis_Group = dplyr::case_when(
      MGUS_Selected ~ "MGUS",
      GEP70_Selected ~ "GEP70>0.66",
      TRUE ~ NA_character_
    )
  )

selected <- strict_rnbx |>
  dplyr::filter(!is.na(Analysis_Group)) |>
  dplyr::mutate(
    Analysis_Group = factor(
      Analysis_Group,
      levels = c("MGUS", "GEP70>0.66")
    )
  ) |>
  dplyr::arrange(Analysis_Group, Sample_ID)

n_mgus <- sum(selected$Analysis_Group == "MGUS")
n_gep70 <- sum(selected$Analysis_Group == "GEP70>0.66")

if (n_mgus < 2 || n_gep70 < 2) {
  stop("Each limma group must contain at least two selected samples.")
}

################################################################################
# 6. HEATMAP SAMPLE SELECTION
################################################################################

heatmap_metadata <- strict_rnbx |>
  dplyr::filter(!is.na(RNBX_Group))

if (!HEATMAP_INCLUDE_ALL_NONMISSING_RNBX_GROUPS) {
  heatmap_metadata <- selected
}

if (nrow(heatmap_metadata) < 2) {
  stop("Fewer than two matched RNBX samples are available.")
}

if (!any(is.finite(heatmap_metadata$Paired_RNAS_PolyIG))) {
  stop("Paired_RNAS_PolyIG contains no finite values after RNAS-to-RNBX pairing.")
}

################################################################################
# 7. PROBE-LEVEL LIMMA
################################################################################

gene <- clean_gene(fdat[[GENE_COL]])

pc_me <- suppressWarnings(
  as.integer(as.character(fdat[[PC_ME_OTHER_COLUMN]]))
)

valid_probe <- !is.na(gene) &
  !is.na(pc_me) &
  pc_me %in% c(0L, 1L, 2L) &
  rowSums(!is.finite(expr)) == 0

probe_expr <- expr[valid_probe, , drop = FALSE]
probe_gene <- gene[valid_probe]
probe_pcme <- pc_me[valid_probe]

selected_ids <- selected$Sample_ID
probe_expr_limma <- probe_expr[, selected_ids, drop = FALSE]

design <- model.matrix(~ 0 + selected$Analysis_Group)
colnames(design) <- c("MGUS", "GEP70_High")
rownames(design) <- selected_ids

if (!identical(colnames(probe_expr_limma), rownames(design))) {
  stop("Expression columns and limma design rows are not aligned.")
}

contrast <- limma::makeContrasts(
  GEP70_High_minus_MGUS = GEP70_High - MGUS,
  levels = design
)

fit <- limma::lmFit(probe_expr_limma, design)
fit <- limma::contrasts.fit(fit, contrast)
fit <- limma::eBayes(fit)

tt <- limma::topTable(
  fit,
  coef = "GEP70_High_minus_MGUS",
  number = Inf,
  sort.by = "none",
  p.value = 1,
  lfc = 0
)

mean_mgus <- rowMeans(
  probe_expr_limma[, selected$Analysis_Group == "MGUS", drop = FALSE]
)

mean_gep70 <- rowMeans(
  probe_expr_limma[, selected$Analysis_Group == "GEP70>0.66", drop = FALSE]
)

probe_results <- tibble::tibble(
  Probe_ID = rownames(probe_expr_limma),
  Gene = probe_gene,
  PC_ME_Other = probe_pcme,
  PC_ME_Other_Label = c(
    "0" = "Other",
    "1" = "PlasmaCells",
    "2" = "Microenvironment"
  )[as.character(probe_pcme)],
  N_MGUS = n_mgus,
  N_GEP70_High = n_gep70,
  Mean_MGUS = mean_mgus,
  Mean_GEP70_High = mean_gep70,
  Log2FC_GEP70_High_vs_MGUS = tt$logFC,
  FoldChange_GEP70_High_over_MGUS = 2^tt$logFC,
  T_Statistic = tt$t,
  B_Statistic = tt$B,
  PValue = tt$P.Value,
  FDR = tt$adj.P.Val
)

################################################################################
# 8. ONE REPRESENTATIVE PROBE PER GENE
################################################################################

best_probe <- probe_results |>
  dplyr::filter(
    !is.na(Gene),
    is.finite(FDR),
    is.finite(PValue),
    is.finite(Log2FC_GEP70_High_vs_MGUS)
  ) |>
  dplyr::group_by(Gene) |>
  dplyr::arrange(
    FDR,
    PValue,
    dplyr::desc(abs(Log2FC_GEP70_High_vs_MGUS)),
    Probe_ID,
    .by_group = TRUE
  ) |>
  dplyr::slice_head(n = 1L) |>
  dplyr::ungroup()

if (!nrow(best_probe)) stop("No representative probes remained.")

best_idx <- match(best_probe$Probe_ID, rownames(probe_expr))

gene_expr <- probe_expr[
  best_idx,
  ,
  drop = FALSE
]

rownames(gene_expr) <- best_probe$Gene

################################################################################
# 9. MICROENVIRONMENT GENE SETS
################################################################################

de <- best_probe |>
  dplyr::mutate(
    ME_GeneSet = dplyr::case_when(
      PC_ME_Other == 2L &
        Log2FC_GEP70_High_vs_MGUS < 0 &
        FDR < MGUS_ME_FDR_CUTOFF ~ "MGUS_ME",
      
      PC_ME_Other == 2L &
        Log2FC_GEP70_High_vs_MGUS > 0 &
        FDR < GEP70_ME_FDR_CUTOFF ~ "GEP70_ME",
      
      TRUE ~ "Not_selected"
    )
  )

mgus_me <- de |>
  dplyr::filter(ME_GeneSet == "MGUS_ME") |>
  dplyr::arrange(FDR, Log2FC_GEP70_High_vs_MGUS, Gene)

gep70_me <- de |>
  dplyr::filter(ME_GeneSet == "GEP70_ME") |>
  dplyr::arrange(FDR, dplyr::desc(Log2FC_GEP70_High_vs_MGUS), Gene)

selected_me <- dplyr::bind_rows(mgus_me, gep70_me)

if (!nrow(selected_me)) {
  stop("No microenvironment genes met the selection rules.")
}

################################################################################
# 10. BUILD HEATMAP MATRIX
################################################################################

heatmap_metadata <- heatmap_metadata |>
  dplyr::mutate(
    RNBX_Group = as.character(RNBX_Group)
  )

heatmap_group_levels <- group_order(
  heatmap_metadata$RNBX_Group
)

heatmap_metadata <- heatmap_metadata |>
  dplyr::mutate(
    RNBX_Group = factor(
      RNBX_Group,
      levels = heatmap_group_levels
    )
  ) |>
  dplyr::arrange(RNBX_Group, Sample_ID)

heatmap_expression <- gene_expr[
  selected_me$Gene,
  heatmap_metadata$Sample_ID,
  drop = FALSE
]

################################################################################
# 11. SUPERVISED ME SCORE AND COLUMN ORDER
################################################################################

if (USE_SUPERVISED_ME_HIGH_MINUS_LOW_SCORE) {
  
  high_genes <- intersect(gep70_me$Gene, rownames(gene_expr))
  low_genes <- intersect(mgus_me$Gene, rownames(gene_expr))
  
  if (!length(high_genes) || !length(low_genes)) {
    stop("Both GEP70_ME and MGUS_ME genes are required for the supervised score.")
  }
  
  high_mean <- colMeans(
    gene_expr[high_genes, heatmap_metadata$Sample_ID, drop = FALSE]
  )
  
  low_mean <- colMeans(
    gene_expr[low_genes, heatmap_metadata$Sample_ID, drop = FALSE]
  )
  
  heatmap_metadata$Mean_GEP70_ME_High_Genes <-
    unname(high_mean[heatmap_metadata$Sample_ID])
  
  heatmap_metadata$Mean_MGUS_ME_Low_Genes <-
    unname(low_mean[heatmap_metadata$Sample_ID])
  
  heatmap_metadata$Supervised_ME_HighMinusLow <-
    heatmap_metadata$Mean_GEP70_ME_High_Genes -
    heatmap_metadata$Mean_MGUS_ME_Low_Genes
  
  if (ORDER_COLUMNS_BY_SUPERVISED_ME_SCORE) {
    
    score_order <- heatmap_metadata$Supervised_ME_HighMinusLow
    
    if (SUPERVISED_ME_SCORE_LEFT_TO_RIGHT == "HighToLow") {
      score_order <- -score_order
    }
    
    idx <- if (SPLIT_SUPERVISED_ORDER_BY_RNBX_GROUP) {
      order(
        heatmap_metadata$RNBX_Group,
        score_order,
        heatmap_metadata$Sample_ID
      )
    } else {
      order(score_order, heatmap_metadata$Sample_ID)
    }
    
    heatmap_metadata <- heatmap_metadata[idx, , drop = FALSE]
    
    heatmap_expression <- heatmap_expression[
      ,
      heatmap_metadata$Sample_ID,
      drop = FALSE
    ]
  }
}

################################################################################
# 12. GROUP SPLITS
################################################################################

heatmap_group_levels <- group_order(
  heatmap_metadata$RNBX_Group
)

heatmap_metadata$RNBX_Group <- factor(
  heatmap_metadata$RNBX_Group,
  levels = heatmap_group_levels
)

group_counts <- table(
  heatmap_metadata$RNBX_Group
)

split_labels <- stats::setNames(
  paste0(
    format_group(names(group_counts)),
    "\n(n = ",
    as.integer(group_counts),
    ")"
  ),
  names(group_counts)
)

heatmap_metadata$Column_Split <- factor(
  split_labels[
    as.character(heatmap_metadata$RNBX_Group)
  ],
  levels = unname(
    split_labels[
      heatmap_group_levels
    ]
  )
)

################################################################################
# 13. ROW TRANSFORMATION
################################################################################

heatmap_transform <- transform_rows(
  heatmap_expression,
  HEATMAP_SCALE_METHOD
)

scaled_heatmap_expression <- heatmap_transform$matrix

row_gene_set <- factor(
  selected_me$ME_GeneSet,
  levels = c("MGUS_ME", "GEP70_ME")
)

################################################################################
# 14. COLORS
################################################################################

heatmap_color_limit <- if (HEATMAP_SCALE_METHOD == "Zscore") {
  Z_SCORE_COLOR_LIMIT
} else {
  max(abs(scaled_heatmap_expression), na.rm = TRUE)
}

if (!is.finite(heatmap_color_limit) || heatmap_color_limit <= 0) {
  heatmap_color_limit <- 1
}

heatmap_color_function <- circlize::colorRamp2(
  c(-heatmap_color_limit, 0, heatmap_color_limit),
  c("#2166AC", "white", "#B2182B")
)

heatmap_legend_at <- seq(
  -heatmap_color_limit,
  heatmap_color_limit,
  length.out = 5
)

row_gene_set_colors <- c(
  "MGUS_ME" = "#0072B2",
  "GEP70_ME" = "#D73027"
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

################################################################################
# 15. SUPERVISED SCORE COLOR SCALE
################################################################################

supervised_color_function <- NULL
supervised_legend_at <- NULL
supervised_legend_labels <- NULL

if (USE_SUPERVISED_ME_HIGH_MINUS_LOW_SCORE) {
  
  limit <- max(
    abs(
      heatmap_metadata$Supervised_ME_HighMinusLow
    ),
    na.rm = TRUE
  )
  
  if (!is.finite(limit) || limit <= 0) {
    limit <- 1
  }
  
  supervised_color_function <- circlize::colorRamp2(
    c(-limit, 0, limit),
    c("#2166AC", "white", "#B2182B")
  )
  
  supervised_legend_at <- seq(
    -limit,
    limit,
    length.out = 5
  )
  
  supervised_legend_labels <- formatC(
    supervised_legend_at,
    format = "fg",
    digits = 3
  )
}

################################################################################
# 16. COLUMN ORDER CONTROLS
################################################################################

use_supervised_column_order <-
  USE_SUPERVISED_ME_HIGH_MINUS_LOW_SCORE &&
  ORDER_COLUMNS_BY_SUPERVISED_ME_SCORE

if (use_supervised_column_order) {
  
  heatmap_column_split <- if (SPLIT_SUPERVISED_ORDER_BY_RNBX_GROUP) {
    heatmap_metadata$Column_Split
  } else {
    NULL
  }
  
  heatmap_cluster_columns <- FALSE
  heatmap_show_column_dendrogram <- FALSE
  
} else {
  
  heatmap_column_split <- heatmap_metadata$Column_Split
  heatmap_cluster_columns <- TRUE
  heatmap_show_column_dendrogram <- TRUE
}

################################################################################
# 17. TOP ANNOTATION
################################################################################

if (USE_SUPERVISED_ME_HIGH_MINUS_LOW_SCORE) {
  
  top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    `GEP Groups` = heatmap_metadata$RNBX_Group,
    PolyIG = heatmap_metadata$Paired_RNAS_PolyIG,
    ME_HighMinusLow = heatmap_metadata$Supervised_ME_HighMinusLow,
    
    col = list(
      `GEP Groups` = group_cols,
      PolyIG = polyig_scale$color_function,
      ME_HighMinusLow = supervised_color_function
    ),
    
    annotation_name_gp = grid::gpar(
      fontsize = 10,
      fontface = "bold"
    ),
    
    simple_anno_size = grid::unit(4, "mm"),
    na_col = "grey85",
    
    annotation_legend_param = list(
      `GEP Groups` = list(
        title = "GEP Groups",
        at = names(group_cols),
        labels = format_group(names(group_cols))
      ),
      
      PolyIG = list(
        title = "Paired RNAS\nPolyIG",
        at = polyig_scale$legend_at,
        labels = polyig_scale$legend_labels
      ),
      
      ME_HighMinusLow = list(
        title = "Supervised ME score\nmean high - mean low",
        at = supervised_legend_at,
        labels = supervised_legend_labels
      )
    )
  )
  
} else {
  
  top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    `GEP Groups` = heatmap_metadata$RNBX_Group,
    PolyIG = heatmap_metadata$Paired_RNAS_PolyIG,
    
    col = list(
      `GEP Groups` = group_cols,
      PolyIG = polyig_scale$color_function
    ),
    
    annotation_name_gp = grid::gpar(
      fontsize = 10,
      fontface = "bold"
    ),
    
    simple_anno_size = grid::unit(4, "mm"),
    na_col = "grey85",
    
    annotation_legend_param = list(
      `GEP Groups` = list(
        title = "GEP Groups",
        at = names(group_cols),
        labels = format_group(names(group_cols))
      ),
      
      PolyIG = list(
        title = "Paired RNAS\nPolyIG",
        at = polyig_scale$legend_at,
        labels = polyig_scale$legend_labels
      )
    )
  )
}

################################################################################
# 18. ROW ANNOTATION
################################################################################

row_annotation <- ComplexHeatmap::rowAnnotation(
  ME_GeneSet = row_gene_set,
  
  col = list(
    ME_GeneSet = row_gene_set_colors
  ),
  
  annotation_name_gp = grid::gpar(
    fontsize = 8,
    fontface = "bold"
  ),
  
  simple_anno_size = grid::unit(4, "mm"),
  
  annotation_legend_param = list(
    ME_GeneSet = list(
      title = "Microenvironment gene set"
    )
  )
)

################################################################################
# 19. CREATE FIGURE 2B HEATMAP
################################################################################

row_name_font_size <- if (nrow(scaled_heatmap_expression) <= 50) {
  6
} else if (nrow(scaled_heatmap_expression) <= 150) {
  4
} else {
  3.5
}

heatmap_title <-
  "Bone-Marrow Microenvironment Gene Signature Across MGUS and NDMM Risk Groups"

figure_2b_heatmap <- ComplexHeatmap::Heatmap(
  scaled_heatmap_expression,
  
  name = if (HEATMAP_SCALE_METHOD == "Zscore") {
    "Row z-score\nof log2 expression"
  } else {
    "Row mean-centered\nlog2 expression"
  },
  
  col = heatmap_color_function,
  top_annotation = top_annotation,
  left_annotation = row_annotation,
  
  row_split = row_gene_set,
  column_split = heatmap_column_split,
  
  cluster_rows = TRUE,
  cluster_row_slices = FALSE,
  cluster_columns = heatmap_cluster_columns,
  cluster_column_slices = FALSE,
  
  show_row_dend = TRUE,
  show_column_dend = heatmap_show_column_dendrogram,
  show_row_names = SHOW_ROW_NAMES,
  show_column_names = SHOW_COLUMN_NAMES,
  
  row_names_gp = grid::gpar(
    fontsize = row_name_font_size
  ),
  
  column_names_gp = grid::gpar(
    fontsize = 5
  ),
  
  row_title_gp = grid::gpar(
    fontsize = 12,
    fontface = "bold"
  ),
  
  column_title = heatmap_title,
  
  column_title_gp = grid::gpar(
    fontsize = 20,
    fontface = "bold"
  ),
  
  column_gap = grid::unit(2.5, "mm"),
  row_gap = grid::unit(2.5, "mm"),
  
  use_raster = ncol(scaled_heatmap_expression) > 200,
  raster_device = "png",
  
  heatmap_legend_param = list(
    title = if (HEATMAP_SCALE_METHOD == "Zscore") {
      "Row z-score\nof log2 expression"
    } else {
      "Row mean-centered\nlog2 expression"
    },
    
    at = heatmap_legend_at,
    
    labels = formatC(
      heatmap_legend_at,
      format = "fg",
      digits = 3
    ),
    
    title_gp = grid::gpar(
      fontface = "bold"
    )
  )
)

ComplexHeatmap::draw(
  figure_2b_heatmap,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legends = FALSE
)

################################################################################
# 20. SAVE FIGURE 2B
################################################################################

tag <- if (HEATMAP_SCALE_METHOD == "Zscore") {
  "Zscore"
} else {
  "MeanCentered"
}

png_file <- file.path(
  OUT_DIR,
  paste0(
    FIGURE_ID,
    "_Microenvironment_Heatmap_",
    tag,
    "_PATID_DiseaseType_PairedRNAS.png"
  )
)

pdf_file <- file.path(
  OUT_DIR,
  paste0(
    FIGURE_ID,
    "_Microenvironment_Heatmap_",
    tag,
    "_PATID_DiseaseType_PairedRNAS.pdf"
  )
)

save_ht_png(
  figure_2b_heatmap,
  png_file
)

save_ht_pdf(
  figure_2b_heatmap,
  pdf_file
)

################################################################################
# 21. KEY SOURCE DATA
################################################################################

write_source(
  probe_results,
  "Fig2B_Lasso_BestProbe_MGUS_vs_GEP70_LIMMA_AllProbeResults.csv"
)

write_source(
  best_probe,
  "Fig2B_Lasso_BestProbe_MGUS_vs_GEP70_LIMMA_BestProbePerGene.csv"
)

write_source(
  mgus_me,
  "Fig2B_Lasso_BestProbe_MGUS_ME_Selected.csv"
)

write_source(
  gep70_me,
  "Fig2B_Lasso_BestProbe_GEP70_ME_Selected.csv"
)

write_source(
  selected_me,
  "Fig2B_Lasso_BestProbe_Selected_Microenvironment_Genes_For_Heatmap.csv"
)

write_source(
  heatmap_metadata,
  "Fig2B_Lasso_BestProbe_Heatmap_Retained_RNBX_SampleMetadata_ColumnOrder.csv"
)

write_source(
  strict_rnbx,
  "Fig2B_Lasso_BestProbe_RNAS_RNBX_PATID_DiseaseType_Pairing_Audit.csv"
)

################################################################################
# 22. SUMMARY
################################################################################

message(
  "\nStep 23C complete",
  "\nRNAS/RNBX pairing: PATID + disease type",
  "\nMGUS LIMMA samples: ", n_mgus,
  "\nGEP70-high LIMMA samples: ", n_gep70,
  "\nProbe-level results: ", nrow(probe_results),
  "\nRepresentative probes/genes: ", nrow(best_probe),
  "\nMGUS_ME genes: ", nrow(mgus_me),
  "\nGEP70_ME genes: ", nrow(gep70_me),
  "\nHeatmap samples: ", nrow(heatmap_metadata),
  "\nHeatmap annotations: GEP Groups, PolyIG, Supervised ME score",
  "\nHeatmap scale: ", HEATMAP_SCALE_METHOD,
  "\nPNG: ", png_file,
  "\nPDF: ", pdf_file,
  "\n"
)