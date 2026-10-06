################################################################################
# STEP 22
# RNBX whole-bone-marrow limma: MGUS vs GEP70-high disease
#
# Purpose:
#   1. Load the Step 21B RNBX-matched UAMS ExpressionSet.
#   2. Select exact one-to-one matched RNBX whole-bone-marrow samples.
#   3. Select MGUS and GEP70-high groups using PCBmBx_All thresholds.
#   4. Collapse probe expression to gene-level mean expression.
#   5. Run limma comparing GEP70-high vs MGUS.
#   6. Generate volcano plots and source-data tables.
#
# Required Step 21 columns:
#   GEP_RNBX_Groups_MGUS_Lasso
#   RNAS_to_RNBX_PATID_match_MGUS_Lasso
#   RNAS_to_RNBX_match_type_Lasso
#
# RNAS_to_RNBX_source_CHIPID_Lasso is NOT required.
################################################################################

suppressPackageStartupMessages({
  library(Biobase)
  library(limma)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tibble)
  library(readr)
})


################################################################################
# 1. CONFIGURATION
################################################################################

FIGURE_ID <- "MGUS_RNBX_vs_GEP70_RNBX"

ESET_PATH <- file.path(
  "data",
  "processed",
  "uams",
  "ESET_uams_all_gep1q_fish_polyig_pcme_binaries_lasso_groups_rnbx-matched.rds"
)

FIGURE_DIR <- file.path(
  "results",
  "figures",
  "Extended",
  FIGURE_ID
)

SOURCE_DATA_DIR <- file.path(
  FIGURE_DIR,
  "source-data"
)

OUT_DIR <- file.path(
  FIGURE_DIR,
  paste0(
    FIGURE_ID,
    "_RNBX_LIMMA"
  )
)

dir.create(
  SOURCE_DATA_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


################################################################################
# 2. EXACT COLUMN NAMES
################################################################################

SAMPLE_ID_COL <- "chipid"
PATID_COL <- "patid"
SAMPLE_GROUP_COL <- "sample_group"

RNBX_GROUP_COL <- "GEP_RNBX_Groups_MGUS_Lasso"

MATCH_STATUS_COL <-
  "RNAS_to_RNBX_PATID_match_MGUS_Lasso"

MATCH_TYPE_COL <-
  "RNAS_to_RNBX_match_type_Lasso"

MATCHED_STATUS_VALUE <-
  "matched_one_to_one"

PCBM_BX_COL <- "PCBmBx_All"

GENE_COL <- "GENE"
PC_ME_OTHER_COL <- "PC_ME_Other"
CHR_ARM_COL <- "CHR_ARM"
LOCUS_TYPE_COL <- "locus_type_n54675"


################################################################################
# 3. ANALYSIS SETTINGS
################################################################################

MGUS_PCBmBx_MAX <- 10

GEP70_PCBmBx_MIN <- 0

GEP70_PCBmBx_MAX_EXCLUSIVE <- 60

VOLCANO_VIEW <- "AllPlots"

FDR_CUTOFF <- 0.05

ABS_LOG2FC_CUTOFF <- 0.001

LABEL_TOP_N_PER_SIDE <- 25

FOCUSED_LABELS_REQUIRE_FDR <- FALSE

SAVE_EXTRA_PLASMA_1Q_IG_PLOT <- TRUE

RNBX_WB_SAMPLE_GROUPS <- c(
  "RNBX_WB_NDMM",
  "RNBX_WB_MGUS"
)


################################################################################
# 4. HELPERS
################################################################################

clean_chr <- function(x) {
  
  x <- trimws(
    as.character(x)
  )
  
  x[
    x %in% c(
      "",
      "NA",
      "N/A",
      "<NA>",
      "NaN"
    )
  ] <- NA_character_
  
  x
}


clean_gene <- function(x) {
  
  x <- clean_chr(x)
  
  x[
    tolower(x) == "no-gene"
  ] <- NA_character_
  
  toupper(x)
}


write_source <- function(
    x,
    file
) {
  
  readr::write_csv(
    x,
    file.path(
      SOURCE_DATA_DIR,
      file
    )
  )
}


collapse_annotation <- function(
    x,
    gene_key,
    valid
) {
  
  values <- split(
    x[valid],
    gene_key[valid]
  )
  
  vapply(
    values,
    function(z) {
      
      z <- unique(
        clean_chr(z)
      )
      
      z <- z[
        !is.na(z)
      ]
      
      if (!length(z)) {
        NA_character_
      } else {
        z[1]
      }
    },
    character(1)
  )
}


select_top_labels <- function(
    dat,
    log2fc_col = "Log2FC_GEP70_High_vs_MGUS",
    fdr_col = "FDR",
    top_n = LABEL_TOP_N_PER_SIDE
) {
  
  neg <- dat[
    is.finite(
      dat[[log2fc_col]]
    ) &
      dat[[log2fc_col]] < 0,
    ,
    drop = FALSE
  ]
  
  pos <- dat[
    is.finite(
      dat[[log2fc_col]]
    ) &
      dat[[log2fc_col]] > 0,
    ,
    drop = FALSE
  ]
  
  neg <- head(
    neg[
      order(
        neg[[fdr_col]],
        neg[[log2fc_col]],
        na.last = TRUE
      ),
      ,
      drop = FALSE
    ],
    top_n
  )
  
  pos <- head(
    pos[
      order(
        pos[[fdr_col]],
        -pos[[log2fc_col]],
        na.last = TRUE
      ),
      ,
      drop = FALSE
    ],
    top_n
  )
  
  rbind(
    neg,
    pos
  )
}


focus_definition <- function(view) {
  
  switch(
    view,
    
    PlasmaCells = list(
      code = 1L,
      name = "Plasma cell genes",
      color = "#0072B2",
      legend = "Plasma cell genes"
    ),
    
    Microenvironment = list(
      code = 2L,
      name = "Microenvironment genes",
      color = "#D73027",
      legend = "Microenvironment genes"
    ),
    
    Other = list(
      code = 0L,
      name = "Other genes",
      color = "#009E73",
      legend = "Other genes"
    ),
    
    NULL
  )
}


################################################################################
# 5. LOAD EXPRESSIONSET
################################################################################

if (!file.exists(ESET_PATH)) {
  stop(
    "Input ExpressionSet not found:\n",
    ESET_PATH
  )
}

eset <- readRDS(
  ESET_PATH
)

if (!methods::is(
  eset,
  "ExpressionSet"
)) {
  stop(
    "Input file is not an ExpressionSet."
  )
}

expr_mat <- Biobase::exprs(
  eset
)

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


################################################################################
# 6. CHECK REQUIRED pData COLUMNS
################################################################################

REQUIRED_PDATA <- c(
  SAMPLE_ID_COL,
  PATID_COL,
  SAMPLE_GROUP_COL,
  RNBX_GROUP_COL,
  MATCH_STATUS_COL,
  MATCH_TYPE_COL,
  PCBM_BX_COL
)

missing_pdata <- setdiff(
  REQUIRED_PDATA,
  colnames(pdat)
)

if (length(missing_pdata) > 0) {
  stop(
    "Missing required pData column(s): ",
    paste(
      missing_pdata,
      collapse = ", "
    )
  )
}


################################################################################
# 7. CHECK REQUIRED fData COLUMNS
################################################################################

REQUIRED_FDATA <- c(
  GENE_COL,
  PC_ME_OTHER_COL,
  CHR_ARM_COL,
  LOCUS_TYPE_COL
)

missing_fdata <- setdiff(
  REQUIRED_FDATA,
  colnames(fdat)
)

if (length(missing_fdata) > 0) {
  stop(
    "Missing required fData column(s): ",
    paste(
      missing_fdata,
      collapse = ", "
    )
  )
}


################################################################################
# 8. ALIGN EXPRESSION MATRIX AND pData
################################################################################

expr_mat <- as.matrix(
  expr_mat
)

storage.mode(
  expr_mat
) <- "numeric"

expr_key <- toupper(
  trimws(
    colnames(expr_mat)
  )
)

pdat_key <- toupper(
  trimws(
    as.character(
      pdat$chipid
    )
  )
)

if (
  anyNA(pdat_key) ||
  any(!nzchar(pdat_key))
) {
  stop(
    "pData chipid contains missing or blank values."
  )
}

if (
  anyDuplicated(expr_key)
) {
  stop(
    "Expression matrix has duplicate sample IDs."
  )
}

if (
  anyDuplicated(pdat_key)
) {
  stop(
    "pData chipid values are duplicated."
  )
}

keep_expr <- expr_key %in% pdat_key

if (!any(keep_expr)) {
  stop(
    "No expression sample IDs matched pData chipid."
  )
}

expr_mat <- expr_mat[
  ,
  keep_expr,
  drop = FALSE
]

pdat <- pdat[
  match(
    toupper(
      trimws(
        colnames(expr_mat)
      )
    ),
    pdat_key
  ),
  ,
  drop = FALSE
]

rownames(pdat) <-
  colnames(expr_mat)

if (
  !identical(
    colnames(expr_mat),
    rownames(pdat)
  )
) {
  stop(
    "Expression columns and pData rows are not aligned."
  )
}


################################################################################
# 9. ALIGN fData
################################################################################

fdat <- fdat[
  match(
    rownames(expr_mat),
    rownames(fdat)
  ),
  ,
  drop = FALSE
]

rownames(fdat) <-
  rownames(expr_mat)


################################################################################
# 10. BUILD SAMPLE METADATA
################################################################################

meta <- data.frame(
  Sample_ID =
    colnames(expr_mat),
  
  PATID =
    pdat$patid,
  
  SAMPLE_GROUP =
    pdat$sample_group,
  
  RNBX_Group =
    as.character(
      pdat[[RNBX_GROUP_COL]]
    ),
  
  Match_Status =
    as.character(
      pdat[[MATCH_STATUS_COL]]
    ),
  
  Match_Type =
    as.character(
      pdat[[MATCH_TYPE_COL]]
    ),
  
  PCBmBx_All =
    pdat$PCBmBx_All,
  
  stringsAsFactors = FALSE
)


################################################################################
# 11. SELECT MATCHED RNBX WHOLE-BONE-MARROW COHORT
################################################################################

meta$RNBX_WB_Eligible <-
  meta$Match_Status ==
  MATCHED_STATUS_VALUE &
  meta$SAMPLE_GROUP %in%
  RNBX_WB_SAMPLE_GROUPS


################################################################################
# 12. VALIDATE DISEASE TYPE
################################################################################

meta$Disease_Type_Consistent <-
  (
    meta$SAMPLE_GROUP ==
      "RNBX_WB_NDMM" &
      meta$Match_Type ==
      "NDMM"
  ) |
  (
    meta$SAMPLE_GROUP ==
      "RNBX_WB_MGUS" &
      meta$Match_Type ==
      "MGUS"
  )

meta$Disease_Type_Consistent[
  is.na(
    meta$Disease_Type_Consistent
  )
] <- FALSE

meta$RNBX_WB_Eligible <-
  meta$RNBX_WB_Eligible &
  meta$Disease_Type_Consistent

matched_idx <-
  !is.na(
    meta$Match_Status
  ) &
  meta$Match_Status ==
  MATCHED_STATUS_VALUE

bad <- meta[
  matched_idx &
    !meta$Disease_Type_Consistent,
  ,
  drop = FALSE
]

if (nrow(bad) > 0) {
  
  write_source(
    bad,
    "FigS4_Matched_RNBX_DiseaseType_Inconsistencies.csv"
  )
  
  stop(
    "Matched RNBX rows have inconsistent disease type."
  )
}


################################################################################
# 13. SELECT MGUS COHORT
################################################################################

meta$MGUS_Selected <-
  meta$RNBX_WB_Eligible &
  meta$RNBX_Group ==
  "MGUS" &
  is.finite(
    meta$PCBmBx_All
  ) &
  meta$PCBmBx_All <=
  MGUS_PCBmBx_MAX


################################################################################
# 14. SELECT GEP70-HIGH COHORT
################################################################################

meta$GEP70_Selected <-
  meta$RNBX_WB_Eligible &
  meta$RNBX_Group ==
  "GEP70>0.66" &
  is.finite(
    meta$PCBmBx_All
  ) &
  meta$PCBmBx_All >=
  GEP70_PCBmBx_MIN &
  meta$PCBmBx_All <
  GEP70_PCBmBx_MAX_EXCLUSIVE

meta$MGUS_Selected[
  is.na(
    meta$MGUS_Selected
  )
] <- FALSE

meta$GEP70_Selected[
  is.na(
    meta$GEP70_Selected
  )
] <- FALSE


################################################################################
# 15. CREATE ANALYSIS GROUP
################################################################################

meta$Analysis_Group <-
  NA_character_

meta$Analysis_Group[
  meta$MGUS_Selected
] <- "MGUS"

meta$Analysis_Group[
  meta$GEP70_Selected
] <- "GEP70>0.66"

if (
  any(
    meta$MGUS_Selected &
    meta$GEP70_Selected
  )
) {
  stop(
    "A sample was assigned to both analysis groups."
  )
}


################################################################################
# 16. SELECT FINAL ANALYSIS SAMPLES
################################################################################

selected <- meta[
  !is.na(
    meta$Analysis_Group
  ),
  ,
  drop = FALSE
]

selected$Analysis_Group <- factor(
  selected$Analysis_Group,
  levels = c(
    "MGUS",
    "GEP70>0.66"
  )
)

selected <- selected[
  order(
    selected$Analysis_Group,
    selected$Sample_ID
  ),
  ,
  drop = FALSE
]

rownames(selected) <- NULL

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
    "Each group needs at least two samples.\n",
    "MGUS: ", n_mgus, "\n",
    "GEP70>0.66: ", n_gep70
  )
}

write_source(
  meta,
  "FigS4_AllSamples_RNBX_Group_PCBmBx_SelectionAudit.csv"
)

write_source(
  selected,
  "FigS4_Selected_MGUS_vs_GEP70_PCBmBx_Samples.csv"
)


################################################################################
# 17. COHORT SIZE PLOT
################################################################################

cohort_count_df <- data.frame(
  Cohort = factor(
    c(
      paste0(
        "MGUS\nPCBmBx_All <= ",
        MGUS_PCBmBx_MAX
      ),
      
      paste0(
        "GEP70>0.66\n",
        GEP70_PCBmBx_MIN,
        " <= PCBmBx_All < ",
        GEP70_PCBmBx_MAX_EXCLUSIVE
      )
    ),
    
    levels = c(
      paste0(
        "MGUS\nPCBmBx_All <= ",
        MGUS_PCBmBx_MAX
      ),
      
      paste0(
        "GEP70>0.66\n",
        GEP70_PCBmBx_MIN,
        " <= PCBmBx_All < ",
        GEP70_PCBmBx_MAX_EXCLUSIVE
      )
    )
  ),
  
  N = c(
    n_mgus,
    n_gep70
  )
)

write_source(
  cohort_count_df,
  "FigS4_Selected_Cohort_Counts.csv"
)

p_FigS4_Cohort_Counts <- ggplot(
  cohort_count_df,
  aes(
    Cohort,
    N,
    fill = Cohort
  )
) +
  geom_col(
    width = 0.70
  ) +
  geom_text(
    aes(
      label = N
    ),
    vjust = -0.35,
    size = 6,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(
      "#0072B2",
      "#D73027"
    )
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(
        0,
        0.15
      )
    )
  ) +
  labs(
    title =
      "Selected RNBX Whole-Bone-Marrow Cohort Sizes",
    x = NULL,
    y = "Number of samples",
    fill = NULL
  ) +
  theme_classic(
    base_size = 18
  ) +
  theme(
    plot.title =
      element_text(
        face = "bold",
        size = 22
      ),
    
    axis.text.x =
      element_text(
        size = 14
      ),
    
    axis.text.y =
      element_text(
        size = 14
      ),
    
    legend.position =
      "none"
  )

print(
  p_FigS4_Cohort_Counts
)

ggsave(
  file.path(
    OUT_DIR,
    "FigS4_Selected_Cohort_Counts.png"
  ),
  p_FigS4_Cohort_Counts,
  width = 8,
  height = 6,
  dpi = 600
)

ggsave(
  file.path(
    OUT_DIR,
    "FigS4_Selected_Cohort_Counts.pdf"
  ),
  p_FigS4_Cohort_Counts,
  width = 8,
  height = 6
)


################################################################################
# 18. COLLAPSE PROBES TO GENE-LEVEL MEAN EXPRESSION
################################################################################

gene_key <- clean_gene(
  fdat$GENE
)

valid_gene <-
  !is.na(
    gene_key
  )

if (!any(valid_gene)) {
  stop(
    "No valid gene annotations remain."
  )
}

selected_expr <- expr_mat[
  ,
  selected$Sample_ID,
  drop = FALSE
]

if (
  any(
    !is.finite(
      selected_expr
    )
  )
) {
  stop(
    "Selected expression matrix contains non-finite values."
  )
}

gene_expr_sum <- rowsum(
  selected_expr[
    valid_gene,
    ,
    drop = FALSE
  ],
  group =
    gene_key[
      valid_gene
    ],
  reorder = FALSE
)

probe_n <- table(
  gene_key[
    valid_gene
  ]
)

gene_expr <- sweep(
  gene_expr_sum,
  1,
  as.numeric(
    probe_n[
      rownames(
        gene_expr_sum
      )
    ]
  ),
  "/"
)


################################################################################
# 19. GENE-LEVEL ANNOTATION
################################################################################

probe_pc <- as.integer(
  fdat$PC_ME_Other
)

if (
  any(
    !is.na(probe_pc) &
    !probe_pc %in%
    c(
      0L,
      1L,
      2L
    )
  )
) {
  stop(
    "PC_ME_Other contains values other than 0, 1, 2, or NA."
  )
}

if (
  anyNA(
    probe_pc[
      valid_gene
    ]
  )
) {
  stop(
    "PC_ME_Other contains missing values for annotated probes."
  )
}

gene_pc <- tapply(
  probe_pc[
    valid_gene
  ],
  gene_key[
    valid_gene
  ],
  function(x) {
    
    values <- unique(x)
    
    if (
      length(values) != 1L
    ) {
      NA_integer_
    } else {
      as.integer(
        values
      )
    }
  }
)

bad_class <- names(
  gene_pc
)[
  is.na(
    gene_pc
  )
]

if (length(bad_class) > 0) {
  
  write_source(
    data.frame(
      Gene = bad_class
    ),
    "FigS4_Genes_With_Inconsistent_PC_ME_Other_Across_Probes.csv"
  )
  
  stop(
    "Some genes have inconsistent PC_ME_Other values across probes."
  )
}

gene_pc <- as.integer(
  gene_pc[
    rownames(
      gene_expr
    )
  ]
)

gene_class <- c(
  "0" = "Other",
  "1" = "Plasma cells",
  "2" = "Microenvironment"
)[
  as.character(
    gene_pc
  )
]

gene_chr <- collapse_annotation(
  fdat$CHR_ARM,
  gene_key,
  valid_gene
)[
  rownames(
    gene_expr
  )
]

gene_locus <- collapse_annotation(
  fdat$locus_type_n54675,
  gene_key,
  valid_gene
)[
  rownames(
    gene_expr
  )
]


################################################################################
# 20. LIMMA
################################################################################

analysis_group <- factor(
  as.character(
    selected$Analysis_Group
  ),
  levels = c(
    "MGUS",
    "GEP70>0.66"
  )
)

design <- model.matrix(
  ~ 0 + analysis_group
)

colnames(design) <- c(
  "MGUS",
  "GEP70_High"
)

rownames(design) <-
  selected$Sample_ID

if (
  !identical(
    colnames(gene_expr),
    rownames(design)
  )
) {
  stop(
    "Expression columns and limma design rows are not aligned."
  )
}

fit <- lmFit(
  gene_expr,
  design
)

fit <- contrasts.fit(
  fit,
  makeContrasts(
    GEP70_High_minus_MGUS =
      GEP70_High - MGUS,
    levels = design
  )
)

fit <- eBayes(
  fit
)

tt <- topTable(
  fit,
  coef =
    "GEP70_High_minus_MGUS",
  number = Inf,
  sort.by = "none",
  p.value = 1,
  lfc = 0
)

tt <- tt[
  match(
    rownames(gene_expr),
    rownames(tt)
  ),
  ,
  drop = FALSE
]

mean_mgus <- rowMeans(
  gene_expr[
    ,
    selected$Analysis_Group ==
      "MGUS",
    drop = FALSE
  ]
)

mean_gep70 <- rowMeans(
  gene_expr[
    ,
    selected$Analysis_Group ==
      "GEP70>0.66",
    drop = FALSE
  ]
)


################################################################################
# 21. DIFFERENTIAL EXPRESSION RESULTS
################################################################################

DE_RESULTS <- data.frame(
  Gene =
    rownames(gene_expr),
  
  PC_ME_Other =
    gene_pc,
  
  PC_ME_Other_Label =
    unname(
      gene_class
    ),
  
  CHR_ARM =
    gene_chr,
  
  locus_type_n54675 =
    gene_locus,
  
  N_MGUS =
    n_mgus,
  
  N_GEP70_High =
    n_gep70,
  
  Mean_MGUS =
    mean_mgus,
  
  Mean_GEP70_High =
    mean_gep70,
  
  Log2FC_GEP70_High_vs_MGUS =
    tt$logFC,
  
  FoldChange_GEP70_High_over_MGUS =
    2^tt$logFC,
  
  T_Statistic =
    tt$t,
  
  B_Statistic =
    tt$B,
  
  PValue =
    tt$P.Value,
  
  FDR =
    tt$adj.P.Val,
  
  stringsAsFactors = FALSE
)

DE_RESULTS$Direction <-
  "Not significant"

DE_RESULTS$Direction[
  DE_RESULTS$FDR <=
    FDR_CUTOFF &
    DE_RESULTS$Log2FC_GEP70_High_vs_MGUS >=
    ABS_LOG2FC_CUTOFF
] <- "Up in GEP70>0.66"

DE_RESULTS$Direction[
  DE_RESULTS$FDR <=
    FDR_CUTOFF &
    DE_RESULTS$Log2FC_GEP70_High_vs_MGUS <=
    -ABS_LOG2FC_CUTOFF
] <- "Up in MGUS"

DE_RESULTS <- DE_RESULTS[
  order(
    DE_RESULTS$FDR,
    -abs(
      DE_RESULTS$Log2FC_GEP70_High_vs_MGUS
    ),
    na.last = TRUE
  ),
  ,
  drop = FALSE
]

write_source(
  DE_RESULTS,
  "FigS4_MGUS_PCBmBx_vs_GEP70_PCBmBx_LIMMA_GeneResults.csv"
)

write_source(
  DE_RESULTS[
    DE_RESULTS$Direction !=
      "Not significant",
    ,
    drop = FALSE
  ],
  "FigS4_MGUS_PCBmBx_vs_GEP70_PCBmBx_LIMMA_SignificantGenes.csv"
)

for (
  pc_code in 0:2
) {
  
  write_source(
    DE_RESULTS[
      DE_RESULTS$PC_ME_Other ==
        pc_code,
      ,
      drop = FALSE
    ],
    paste0(
      "FigS4_MGUS_PCBmBx_vs_GEP70_PCBmBx_LIMMA_PC_ME_Other_",
      pc_code,
      "_GeneResults.csv"
    )
  )
}

write_source(
  data.frame(
    Sample_ID =
      rownames(design),
    
    design,
    
    check.names = FALSE
  ),
  "FigS4_LIMMA_DesignMatrix.csv"
)


################################################################################
# 22. VOLCANO DATA
################################################################################

volcano_dat <- DE_RESULTS[
  is.finite(
    DE_RESULTS$Log2FC_GEP70_High_vs_MGUS
  ) &
    is.finite(
      DE_RESULTS$FDR
    ),
  ,
  drop = FALSE
]

volcano_dat$MinusLog10FDR <- -log10(
  pmax(
    volcano_dat$FDR,
    .Machine$double.xmin
  )
)


################################################################################
# 23. VOLCANO FUNCTION
################################################################################

make_volcano <- function(
    dat,
    view
) {
  
  if (
    view == "All"
  ) {
    
    plot_dat <- dat
    
    plot_dat$Direction <- factor(
      plot_dat$Direction,
      levels = c(
        "Up in MGUS",
        "Not significant",
        "Up in GEP70>0.66"
      )
    )
    
    labels <- select_top_labels(
      plot_dat[
        plot_dat$Direction !=
          "Not significant",
        ,
        drop = FALSE
      ]
    )
    
    return(
      ggplot(
        plot_dat,
        aes(
          Log2FC_GEP70_High_vs_MGUS,
          MinusLog10FDR,
          color = Direction
        )
      ) +
        geom_point(
          alpha = 0.75,
          size = 1.5
        ) +
        geom_vline(
          xintercept = c(
            -ABS_LOG2FC_CUTOFF,
            ABS_LOG2FC_CUTOFF
          ),
          linetype = "dashed",
          linewidth = 0.45
        ) +
        geom_hline(
          yintercept =
            -log10(
              FDR_CUTOFF
            ),
          linetype = "dashed",
          linewidth = 0.45
        ) +
        ggrepel::geom_text_repel(
          data = labels,
          aes(
            label = Gene
          ),
          size = 3.7,
          max.overlaps = Inf,
          box.padding = 0.35,
          point.padding = 0.20,
          show.legend = FALSE
        ) +
        scale_color_manual(
          values = c(
            "Up in MGUS" =
              "#0072B2",
            
            "Not significant" =
              "grey70",
            
            "Up in GEP70>0.66" =
              "#D73027"
          ),
          drop = FALSE
        ) +
        labs(
          title =
            "MGUS vs GEP70-High",
          
          x =
            "log2 fold change (GEP70>0.66 vs MGUS)",
          
          y =
            expression(
              -log[10](
                "BH-adjusted P value"
              )
            ),
          
          color = NULL
        ) +
        theme_classic(
          base_size = 17
        ) +
        theme(
          plot.title =
            element_text(
              face = "bold",
              size = 22
            ),
          
          legend.position =
            "top",
          
          legend.text =
            element_text(
              size = 13
            )
        )
    )
  }
  
  
  fd <- focus_definition(
    view
  )
  
  focus <- dat[
    dat$PC_ME_Other ==
      fd$code,
    ,
    drop = FALSE
  ]
  
  if (!nrow(focus)) {
    return(NULL)
  }
  
  background <- dat
  
  background$Volcano_Legend <-
    "All genes"
  
  focus$Volcano_Legend <-
    fd$legend
  
  labels <- focus
  
  if (
    FOCUSED_LABELS_REQUIRE_FDR
  ) {
    
    labels <- labels[
      labels$FDR <=
        FDR_CUTOFF,
      ,
      drop = FALSE
    ]
  }
  
  labels <- select_top_labels(
    labels
  )
  
  colors <- c(
    "All genes" =
      "grey70",
    
    setNames(
      fd$color,
      fd$legend
    )
  )
  
  ggplot(
    dat,
    aes(
      Log2FC_GEP70_High_vs_MGUS,
      MinusLog10FDR
    )
  ) +
    geom_point(
      data = background,
      aes(
        color =
          Volcano_Legend
      ),
      alpha = 0.30,
      size = 1.35
    ) +
    geom_point(
      data = focus,
      aes(
        color =
          Volcano_Legend
      ),
      alpha = 0.90,
      size = 1.65
    ) +
    geom_vline(
      xintercept = c(
        -ABS_LOG2FC_CUTOFF,
        ABS_LOG2FC_CUTOFF
      ),
      linetype = "dashed",
      linewidth = 0.45
    ) +
    geom_hline(
      yintercept =
        -log10(
          FDR_CUTOFF
        ),
      linetype = "dashed",
      linewidth = 0.45
    ) +
    ggrepel::geom_text_repel(
      data = labels,
      aes(
        label = Gene
      ),
      color = fd$color,
      size = 3.7,
      max.overlaps = Inf,
      box.padding = 0.35,
      point.padding = 0.20,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = colors
    ) +
    guides(
      color = guide_legend(
        title = NULL,
        override.aes = list(
          alpha = c(
            0.30,
            0.90
          ),
          size = c(
            2,
            2
          )
        )
      )
    ) +
    labs(
      title = paste0(
        "MGUS vs GEP70-High: ",
        fd$name
      ),
      
      x =
        "log2 fold change (GEP70-High vs MGUS)",
      
      y =
        expression(
          -log[10](
            "BH-adjusted P value"
          )
        )
    ) +
    theme_classic(
      base_size = 17
    ) +
    theme(
      plot.title =
        element_text(
          face = "bold",
          size = 22
        ),
      
      legend.position =
        "top",
      
      legend.text =
        element_text(
          size = 13
        )
    )
}


################################################################################
# 24. EXTRA PLASMA-CELL 1q / IMMUNOGLOBULIN VOLCANO
################################################################################

make_plasma_1q_ig <- function(dat) {
  
  plasma <- dat[
    dat$PC_ME_Other ==
      1,
    ,
    drop = FALSE
  ]
  
  if (!nrow(plasma)) {
    return(NULL)
  }
  
  is_1q <-
    !is.na(
      plasma$CHR_ARM
    ) &
    tolower(
      trimws(
        plasma$CHR_ARM
      )
    ) == "1q"
  
  is_ig <-
    !is.na(
      plasma$locus_type_n54675
    ) &
    tolower(
      trimws(
        plasma$locus_type_n54675
      )
    ) == "immunoglobulin"
  
  oneq <- plasma[
    is_1q,
    ,
    drop = FALSE
  ]
  
  ig <- plasma[
    is_ig,
    ,
    drop = FALSE
  ]
  
  labels <- plasma[
    (
      is_1q |
        is_ig
    ) &
      plasma$FDR <=
      FDR_CUTOFF,
    ,
    drop = FALSE
  ]
  
  labels <- select_top_labels(
    labels
  )
  
  labels$Label_Group <- ifelse(
    !is.na(
      labels$CHR_ARM
    ) &
      tolower(
        trimws(
          labels$CHR_ARM
        )
      ) == "1q",
    
    "1q genes",
    
    "Immunoglobulin genes"
  )
  
  ggplot(
    plasma,
    aes(
      Log2FC_GEP70_High_vs_MGUS,
      MinusLog10FDR
    )
  ) +
    geom_point(
      aes(
        color =
          "All plasma cell genes"
      ),
      alpha = 0.30,
      size = 1.35
    ) +
    geom_point(
      data = oneq,
      aes(
        color =
          "1q genes"
      ),
      alpha = 0.95,
      size = 1.75
    ) +
    geom_point(
      data = ig,
      aes(
        color =
          "Immunoglobulin genes"
      ),
      alpha = 0.95,
      size = 1.75
    ) +
    geom_vline(
      xintercept = c(
        -ABS_LOG2FC_CUTOFF,
        ABS_LOG2FC_CUTOFF
      ),
      linetype = "dashed",
      linewidth = 0.45
    ) +
    geom_hline(
      yintercept =
        -log10(
          FDR_CUTOFF
        ),
      linetype = "dashed",
      linewidth = 0.45
    ) +
    ggrepel::geom_text_repel(
      data = labels,
      aes(
        label = Gene,
        color = Label_Group
      ),
      size = 3.7,
      max.overlaps = Inf,
      box.padding = 0.35,
      point.padding = 0.20,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = c(
        "All plasma cell genes" =
          "grey70",
        
        "1q genes" =
          "#D73027",
        
        "Immunoglobulin genes" =
          "#0072B2"
      )
    ) +
    guides(
      color = guide_legend(
        title = NULL,
        override.aes = list(
          alpha = c(
            0.30,
            0.95,
            0.95
          ),
          size = c(
            2,
            2,
            2
          )
        )
      )
    ) +
    labs(
      title = paste0(
        "MGUS vs GEP70-High\n",
        "Chromosome 1q and Immunoglobulin Genes Highlighted"
      ),
      
      x =
        "log2 fold change (GEP70-High vs MGUS)",
      
      y =
        expression(
          -log[10](
            "BH-adjusted P value"
          )
        )
    ) +
    theme_classic(
      base_size = 17
    ) +
    theme(
      plot.title =
        element_text(
          face = "bold",
          size = 22
        ),
      
      legend.position =
        "top",
      
      legend.text =
        element_text(
          size = 13
        )
    )
}


################################################################################
# 25. SAVE VOLCANO PLOTS
################################################################################

views <- if (
  VOLCANO_VIEW ==
  "AllPlots"
) {
  
  c(
    "All",
    "PlasmaCells",
    "Microenvironment",
    "Other"
  )
  
} else {
  
  VOLCANO_VIEW
}

saved <- list()

for (
  view in views
) {
  
  plot <- make_volcano(
    volcano_dat,
    view
  )
  
  if (
    is.null(plot)
  ) {
    
    saved[[length(saved) + 1]] <- data.frame(
      Volcano_View = view,
      PNG_File = NA_character_,
      PDF_File = NA_character_,
      Status = "Skipped",
      stringsAsFactors = FALSE
    )
    
    next
  }
  
  print(
    plot
  )
  
  stub <- paste0(
    "FigS4_MGUS_PCBmBx_vs_GEP70_PCBmBx_",
    view,
    "_Volcano"
  )
  
  ggsave(
    file.path(
      OUT_DIR,
      paste0(
        stub,
        ".png"
      )
    ),
    plot,
    width = 11,
    height = 8.5,
    dpi = 600
  )
  
  ggsave(
    file.path(
      OUT_DIR,
      paste0(
        stub,
        ".pdf"
      )
    ),
    plot,
    width = 11,
    height = 8.5
  )
  
  saved[[length(saved) + 1]] <- data.frame(
    Volcano_View = view,
    PNG_File = paste0(
      stub,
      ".png"
    ),
    PDF_File = paste0(
      stub,
      ".pdf"
    ),
    Status = "Saved",
    stringsAsFactors = FALSE
  )
}


################################################################################
# 26. SAVE EXTRA PLASMA-CELL 1q / IG VOLCANO
################################################################################

if (
  isTRUE(
    SAVE_EXTRA_PLASMA_1Q_IG_PLOT
  )
) {
  
  plot <- make_plasma_1q_ig(
    volcano_dat
  )
  
  if (
    !is.null(plot)
  ) {
    
    print(
      plot
    )
    
    png <-
      "FigS4_MGUS_PCBmBx_vs_GEP70_PCBmBx_PlasmaCells_1q_Immunoglobulin_Volcano.png"
    
    pdf <-
      "FigS4_MGUS_PCBmBx_vs_GEP70_PCBmBx_PlasmaCells_1q_Immunoglobulin_Volcano.pdf"
    
    ggsave(
      file.path(
        OUT_DIR,
        png
      ),
      plot,
      width = 11,
      height = 8.5,
      dpi = 600
    )
    
    ggsave(
      file.path(
        OUT_DIR,
        pdf
      ),
      plot,
      width = 11,
      height = 8.5
    )
    
    saved[[length(saved) + 1]] <- data.frame(
      Volcano_View =
        "PlasmaCells_1q_Immunoglobulin",
      
      PNG_File =
        png,
      
      PDF_File =
        pdf,
      
      Status =
        "Saved",
      
      stringsAsFactors =
        FALSE
    )
  }
}


################################################################################
# 27. SAVE VOLCANO MANIFEST
################################################################################

write_source(
  do.call(
    rbind,
    saved
  ),
  "FigS4_Volcano_Plots_Saved.csv"
)


################################################################################
# 28. RUN SUMMARY
################################################################################

write_source(
  tibble(
    Step = 22L,
    
    Figure =
      FIGURE_ID,
    
    Input_ESET =
      ESET_PATH,
    
    RNBX_Group_Column =
      RNBX_GROUP_COL,
    
    Match_Status_Column =
      MATCH_STATUS_COL,
    
    Match_Type_Column =
      MATCH_TYPE_COL,
    
    Match_status_required =
      MATCHED_STATUS_VALUE,
    
    NDMM_source_rule =
      "RNAS_CD138_NDMM",
    
    NDMM_target_rule =
      "RNBX_WB_NDMM",
    
    MGUS_source_rule =
      "RNAS_CD138_MGUS",
    
    MGUS_target_rule =
      "RNBX_WB_MGUS",
    
    Match_key =
      "PATID + disease type",
    
    n_MGUS =
      n_mgus,
    
    n_GEP70_high =
      n_gep70,
    
    MGUS_PCBmBx_max =
      MGUS_PCBmBx_MAX,
    
    GEP70_PCBmBx_min =
      GEP70_PCBmBx_MIN,
    
    GEP70_PCBmBx_max_exclusive =
      GEP70_PCBmBx_MAX_EXCLUSIVE,
    
    FDR_cutoff =
      FDR_CUTOFF,
    
    abs_log2FC_cutoff =
      ABS_LOG2FC_CUTOFF,
    
    volcano_view =
      VOLCANO_VIEW
  ),
  
  "FigS4_Run_Summary.csv"
)


################################################################################
# 29. COMPLETE
################################################################################

cat(
  "\nStep 22 complete.\n",
  "MGUS: ",
  n_mgus,
  "\n",
  "GEP70>0.66: ",
  n_gep70,
  "\n",
  "Output: ",
  normalizePath(
    OUT_DIR,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n",
  sep = ""
)