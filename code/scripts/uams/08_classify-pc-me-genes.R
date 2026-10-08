################################################################################
# STEP 08 — CLASSIFY PLASMA-CELL AND MICROENVIRONMENT GENES
#
# Purpose:
#   1. Match NDMM RNAS CD138 and RNBX whole-bone-marrow samples by PATID.
#   2. Collapse multiple probes to gene-level mean expression.
#   3. Perform paired limma differential expression:
#        RNBX_WB - RNAS_CD138
#   4. Classify genes as:
#        1 = Plasma-cell
#        2 = Microenvironment
#        0 = Other
#   5. Add classifications to fData.
#   6. Save paired results, volcano plot, and updated ExpressionSet.
#
# Required first:
#   Run through code/run_pipeline.R so DATA_ROOT and RESULTS_ROOT are defined.
################################################################################

library(limma)
library(ggplot2)


################################################################################
# 1. VERIFY PIPELINE ROOTS
################################################################################

if (!exists("DATA_ROOT")) {
  
  stop(
    "DATA_ROOT is not defined. ",
    "Run this script through code/run_pipeline.R."
  )
}

if (!exists("RESULTS_ROOT")) {
  
  stop(
    "RESULTS_ROOT is not defined. ",
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
  RESULTS_ROOT,
  recursive = TRUE,
  showWarnings = FALSE
)


################################################################################
# 2. FILES
################################################################################

ESET_FILE <- file.path(
  DATA_ROOT,
  "processed",
  "ESET_uams_all_gep1q_fish_polyig.rds"
)

OUTPUT_FILE <- file.path(
  DATA_ROOT,
  "processed",
  "ESET_uams_all_gep1q_fish_polyig_pcme.rds"
)

OUT_DIR <- file.path(
  RESULTS_ROOT,
  "figures",
  "Extended",
  "PlasmaCell-Microenvironment-genes"
)

if (!file.exists(ESET_FILE)) {
  
  stop(
    "PolyIG-scored ExpressionSet not found: ",
    ESET_FILE
  )
}

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


################################################################################
# 3. SETTINGS
################################################################################

P_VALUE_CUTOFF <- 0.05
FDR_CUTOFF <- 0.05
LOG2FC_CUTOFF <- 0


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
# 5. COLLAPSE PROBES TO GENE-LEVEL MEAN EXPRESSION
################################################################################

gene <- trimws(
  as.character(
    features$GENE
  )
)

keep_gene <- !is.na(gene) &
  gene != "" &
  tolower(gene) != "no-gene"

gene_key <- toupper(
  gene[
    keep_gene
  ]
)

gene_sum <- rowsum(
  expression[
    keep_gene,
    ,
    drop = FALSE
  ],
  group = gene_key
)

gene_count <- table(
  gene_key
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
# 6. SELECT NDMM RNAS AND RNBX SAMPLES
################################################################################

rnas <- phenotype[
  phenotype$sample_group == "RNAS_CD138_NDMM",
  ,
  drop = FALSE
]

rnbx <- phenotype[
  phenotype$sample_group == "RNBX_WB_NDMM",
  ,
  drop = FALSE
]

rnas$Sample_ID <- rownames(
  rnas
)

rnbx$Sample_ID <- rownames(
  rnbx
)

rnas$PATID <- as.character(
  rnas$patid
)

rnbx$PATID <- as.character(
  rnbx$patid
)

rnas <- rnas[
  !is.na(rnas$PATID) &
    rnas$PATID != "",
  ,
  drop = FALSE
]

rnbx <- rnbx[
  !is.na(rnbx$PATID) &
    rnbx$PATID != "",
  ,
  drop = FALSE
]


################################################################################
# 7. REQUIRE ONE RNAS AND ONE RNBX SAMPLE PER PATID
################################################################################

if (anyDuplicated(rnas$PATID)) {
  
  stop(
    "More than one RNAS CD138 NDMM sample exists for at least one PATID."
  )
}

if (anyDuplicated(rnbx$PATID)) {
  
  stop(
    "More than one RNBX NDMM sample exists for at least one PATID."
  )
}


################################################################################
# 8. MATCH RNAS AND RNBX SAMPLES BY PATID
################################################################################

paired <- merge(
  rnas[
    ,
    c(
      "PATID",
      "Sample_ID"
    )
  ],
  rnbx[
    ,
    c(
      "PATID",
      "Sample_ID"
    )
  ],
  by = "PATID",
  suffixes = c(
    "_RNAS",
    "_RNBX"
  )
)

paired <- paired[
  order(
    paired$PATID
  ),
  ,
  drop = FALSE
]

cat(
  "Paired NDMM samples:",
  nrow(paired),
  "\n"
)

write.csv(
  paired,
  file.path(
    OUT_DIR,
    "FigS4_RNAS_RNBX_Pairs.csv"
  ),
  row.names = FALSE
)


################################################################################
# 9. EXTRACT PAIRED GENE-EXPRESSION MATRICES
################################################################################

rnas_matrix <- gene_expression[
  ,
  paired$Sample_ID_RNAS,
  drop = FALSE
]

rnbx_matrix <- gene_expression[
  ,
  paired$Sample_ID_RNBX,
  drop = FALSE
]


################################################################################
# 10. BUILD PAIRED LIMMA MODEL
#
# Coefficient:
#   RNBX_WB - RNAS_CD138
################################################################################

pair_id <- rep(
  paired$PATID,
  2
)

compartment <- factor(
  c(
    rep(
      "RNBX_WB",
      nrow(paired)
    ),
    rep(
      "RNAS_CD138",
      nrow(paired)
    )
  ),
  levels = c(
    "RNAS_CD138",
    "RNBX_WB"
  )
)

limma_expression <- cbind(
  rnbx_matrix,
  rnas_matrix
)

design <- model.matrix(
  ~ 0 + factor(pair_id) + compartment
)


################################################################################
# 11. RUN PAIRED LIMMA ANALYSIS
################################################################################

fit <- lmFit(
  limma_expression,
  design
)

fit <- eBayes(
  fit
)

limma_results <- topTable(
  fit,
  coef = "compartmentRNBX_WB",
  number = Inf,
  sort.by = "none"
)

results <- data.frame(
  Gene = rownames(
    limma_results
  ),
  Log2FC_RNBX_vs_RNAS = limma_results$logFC,
  PValue = limma_results$P.Value,
  FDR = limma_results$adj.P.Val,
  stringsAsFactors = FALSE
)


################################################################################
# 12. CLASSIFY PLASMA-CELL AND MICROENVIRONMENT GENES
#
# 1 = Plasma-cell:
#       log2FC < 0 and nominal P < 0.05
#
# 2 = Microenvironment:
#       log2FC > 0 and nominal P < 0.05
#
# 0 = Other
################################################################################

results$PC_ME_Other <- 0L

results$PC_ME_Other[
  results$Log2FC_RNBX_vs_RNAS < LOG2FC_CUTOFF &
    results$PValue < P_VALUE_CUTOFF
] <- 1L

results$PC_ME_Other[
  results$Log2FC_RNBX_vs_RNAS > LOG2FC_CUTOFF &
    results$PValue < P_VALUE_CUTOFF
] <- 2L

results$PlasmaCellGenes <- as.integer(
  results$PC_ME_Other == 1
)

results$MicroenvironmentGenes <- as.integer(
  results$PC_ME_Other == 2
)

results$PC_ME_Label <- c(
  "Other",
  "Plasma cells",
  "Microenvironment"
)[
  results$PC_ME_Other + 1
]


################################################################################
# 13. SAVE GENE-LEVEL RESULTS
################################################################################

write.csv(
  results,
  file.path(
    OUT_DIR,
    "FigS4_RNBX_vs_RNAS_GeneResults.csv"
  ),
  row.names = FALSE
)


################################################################################
# 14. DEFINE VOLCANO-PLOT DIRECTION
#
# Up in RNBX:
#   log2FC > 0 and FDR < 0.05
#
# Up in RNAS:
#   log2FC < 0 and FDR < 0.05
################################################################################

results$Direction <- "Not significant"

results$Direction[
  results$FDR < FDR_CUTOFF &
    results$Log2FC_RNBX_vs_RNAS > LOG2FC_CUTOFF
] <- "Up in RNBX"

results$Direction[
  results$FDR < FDR_CUTOFF &
    results$Log2FC_RNBX_vs_RNAS < LOG2FC_CUTOFF
] <- "Up in RNAS"

results$MinusLog10FDR <- -log10(
  pmax(
    results$FDR,
    .Machine$double.xmin
  )
)


################################################################################
# 15. CREATE VOLCANO PLOT
################################################################################

p <- ggplot(
  results,
  aes(
    x = Log2FC_RNBX_vs_RNAS,
    y = MinusLog10FDR,
    color = Direction
  )
) +
  geom_point(
    alpha = 0.75,
    size = 1.5
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(
      FDR_CUTOFF
    ),
    linetype = "dashed"
  ) +
  scale_color_manual(
    values = c(
      "Up in RNAS" = "#0072B2",
      "Not significant" = "grey70",
      "Up in RNBX" = "#D73027"
    )
  ) +
  labs(
    title = "Paired NDMM RNBX vs RNAS",
    subtitle = paste0(
      "PATID-matched pairs: n = ",
      nrow(paired)
    ),
    x = "log2 fold change: RNBX / RNAS",
    y = expression(
      -log[10]("BH-adjusted P value")
    ),
    color = NULL
  ) +
  theme_classic(
    base_size = 16
  )

print(
  p
)


################################################################################
# 16. SAVE VOLCANO PLOT
################################################################################

ggsave(
  file.path(
    OUT_DIR,
    "FigS4_RNBX_vs_RNAS_Volcano.png"
  ),
  p,
  width = 9,
  height = 7,
  dpi = 600
)

ggsave(
  file.path(
    OUT_DIR,
    "FigS4_RNBX_vs_RNAS_Volcano.pdf"
  ),
  p,
  width = 9,
  height = 7
)


################################################################################
# 17. MAP GENE CLASSIFICATIONS BACK TO fData
################################################################################

feature_gene_key <- toupper(
  trimws(
    as.character(
      features$GENE
    )
  )
)

result_match <- match(
  feature_gene_key,
  results$Gene
)

features$PC_ME_Other <- 0L

matched <- !is.na(
  result_match
)

features$PC_ME_Other[
  matched
] <- results$PC_ME_Other[
  result_match[
    matched
  ]
]

features$PlasmaCellGenes <- as.integer(
  features$PC_ME_Other == 1
)

features$MicroenvironmentGenes <- as.integer(
  features$PC_ME_Other == 2
)


################################################################################
# 18. UPDATE EXPRESSIONSET
################################################################################

fData(
  eset
) <- features


################################################################################
# 19. SAVE UPDATED EXPRESSIONSET
################################################################################

dir.create(
  dirname(OUTPUT_FILE),
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  eset,
  OUTPUT_FILE
)


################################################################################
# 20. COMPLETE
################################################################################

cat(
  "Paired NDMM samples:",
  nrow(paired),
  "\n"
)

cat(
  "Plasma-cell genes:",
  sum(
    results$PC_ME_Other == 1
  ),
  "\n"
)

cat(
  "Microenvironment genes:",
  sum(
    results$PC_ME_Other == 2
  ),
  "\n"
)

cat(
  "Updated ExpressionSet saved:",
  OUTPUT_FILE,
  "\n"
)

cat(
  "Results saved:",
  OUT_DIR,
  "\n"
)