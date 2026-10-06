################################################################################
# STEP 11 — FIGURE 1C
#
# Purpose:
#   1. Recreate the Figure 1B TTP cohort.
#   2. Select plasma-cell probes with Cox FDR < 0.01.
#   3. Keep ALL qualifying probes.
#   4. Row Z-score the expression matrix.
#   5. Run consensus clustering.
#   6. Use K = 4 for Figure 1C.
#   7. Create the Figure 1C heatmap.
#
# Important:
#   Multiple probes mapping to the same gene are NOT collapsed.
################################################################################


################################################################################
# 1. Packages
################################################################################

library(Biobase)
library(dplyr)
library(readr)
library(matrixStats)
library(ConsensusClusterPlus)
library(ComplexHeatmap)
library(circlize)


################################################################################
# 2. Files and settings
################################################################################

ESET_FILE <- file.path(
  "data",
  "processed",
  "uams",
  "ESET_uams_all_gep1q_fish_polyig_pcme.rds"
)

COX_FILE <- file.path(
  "results",
  "figures",
  "Fig1",
  "Fig1B",
  "Fig1B_Cox_TTP_all-probe-results.csv"
)

OUT_DIR <- file.path(
  "results",
  "figures",
  "Fig1",
  "Fig1C_CC"
)

CCP_DIR <- file.path(
  OUT_DIR,
  "CCP"
)

HM_DIR <- file.path(
  OUT_DIR,
  "HM"
)

FDR_CUTOFF <- 0.01
PC_ME_OTHER_KEEP <- 1
CLIP_Z <- 2
MAX_K <- 8
REPS <- 1000
PITEM <- 0.80
PFEATURE <- 1.00
K_USE <- 4
DISTANCE <- "pearson"
CLUSTER_ALG <- "hc"
LINKAGE <- "ward.D2"

HEATMAP_CLUSTER_ORDER <- c(
  "1",
  "2",
  "3",
  "4"
)

MANUAL_CLUSTER_LABELS <- c(
  "1" = "Cluster 1",
  "2" = "Cluster 2",
  "3" = "Cluster 3",
  "4" = "Cluster 4"
)


################################################################################
# 3. Create output folders
################################################################################

dir.create(
  CCP_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  HM_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


################################################################################
# 4. Load ExpressionSet
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

fdata <- fData(
  eset
)


################################################################################
# 5. Recreate the Figure 1B TTP cohort
################################################################################

keep <- phenotype$sample_group == "RNAS_CD138_NDMM" &
  phenotype$tt_on_for_tx1 %in% c(
    "TT2_NoThal",
    "TT2_Thal",
    "TT3a",
    "TT3b",
    "TT4_S-TT3",
    "TT5"
  ) &
  as.numeric(
    phenotype$GEP70
  ) <= 0.66

keep[is.na(keep)] <- FALSE

phenotype_model <- phenotype[
  keep,
  ,
  drop = FALSE
]

phenotype_model$time <- as.numeric(
  phenotype_model$YearsTTP
)

phenotype_model$event <- as.numeric(
  phenotype_model$CensTTP
)

valid <- is.finite(
  phenotype_model$time
) &
  phenotype_model$time > 0 &
  phenotype_model$event %in% c(
    0,
    1
  )

phenotype_model <- phenotype_model[
  valid,
  ,
  drop = FALSE
]

sample_ids <- rownames(
  phenotype_model
)

expression <- expression[
  ,
  sample_ids,
  drop = FALSE
]

cat(
  "Samples:",
  ncol(expression),
  "\n"
)

cat(
  "Events:",
  sum(
    phenotype_model$event == 1
  ),
  "\n"
)


################################################################################
# 6. Load Cox results
################################################################################

cox <- read_csv(
  COX_FILE,
  show_col_types = FALSE
)


################################################################################
# 7. Select plasma-cell probes
#
#   PC_ME_Other == 1
#   Cox FDR < 0.01
#
#   ALL qualifying probes are retained.
################################################################################

selected <- cox |>
  filter(
    FDR < FDR_CUTOFF,
    PC_ME_Other == PC_ME_OTHER_KEEP
  ) |>
  arrange(
    FDR,
    p_value
  )

cat(
  "Selected probes:",
  nrow(selected),
  "\n"
)


################################################################################
# 8. Match selected probes to expression matrix
################################################################################

probe_index <- match(
  selected$Probe_ID,
  rownames(expression)
)

keep_probe <- !is.na(
  probe_index
)

selected <- selected[
  keep_probe,
  ,
  drop = FALSE
]

probe_index <- probe_index[
  keep_probe
]

expression_selected <- expression[
  probe_index,
  ,
  drop = FALSE
]

rownames(
  expression_selected
) <- selected$Probe_ID


################################################################################
# 9. Remove probes with non-finite expression
################################################################################

keep_probe <- apply(
  expression_selected,
  1,
  function(x) {
    all(
      is.finite(x)
    )
  }
)

expression_selected <- expression_selected[
  keep_probe,
  ,
  drop = FALSE
]

selected <- selected[
  keep_probe,
  ,
  drop = FALSE
]

cat(
  "Probes used for clustering:",
  nrow(expression_selected),
  "\n"
)


################################################################################
# 10. Row Z-score
################################################################################

row_mean <- rowMeans(
  expression_selected
)

row_sd <- matrixStats::rowSds(
  expression_selected
)

row_sd[
  !is.finite(row_sd) |
    row_sd == 0
] <- 1

zmat <- sweep(
  expression_selected,
  1,
  row_mean,
  "-"
)

zmat <- sweep(
  zmat,
  1,
  row_sd,
  "/"
)

zmat[
  zmat > CLIP_Z
] <- CLIP_Z

zmat[
  zmat < -CLIP_Z
] <- -CLIP_Z


################################################################################
# 11. Consensus clustering
################################################################################

set.seed(
  123
)

ccp <- ConsensusClusterPlus(
  d = zmat,
  maxK = MAX_K,
  reps = REPS,
  pItem = PITEM,
  pFeature = PFEATURE,
  clusterAlg = CLUSTER_ALG,
  distance = DISTANCE,
  innerLinkage = LINKAGE,
  finalLinkage = LINKAGE,
  plot = "png",
  title = CCP_DIR,
  seed = 123
)

################################################################################
# 12. Extract K = 4 clusters
################################################################################

clusters <- ccp[[K_USE]]$consensusClass

clusters <- as.integer(
  clusters[
    colnames(zmat)
  ]
)

cluster_df <- data.frame(
  CHIPID = colnames(zmat),
  Cluster = clusters
)


################################################################################
# 13. Save cluster assignments
################################################################################

write_csv(
  cluster_df,
  file.path(
    OUT_DIR,
    "Fig1C_ConsensusClusters_TTP_K4_Assignments.csv"
  )
)


################################################################################
# 14. Order samples
################################################################################

cluster_ordered <- factor(
  as.character(
    clusters
  ),
  levels = HEATMAP_CLUSTER_ORDER
)

column_order <- order(
  cluster_ordered
)

zmat_ordered <- zmat[
  ,
  column_order,
  drop = FALSE
]

cluster_ordered <- cluster_ordered[
  column_order
]


################################################################################
# 15. Cluster annotation
################################################################################

cluster_colors <- c(
  "Cluster 1" = "#E41A1C",
  "Cluster 2" = "#377EB8",
  "Cluster 3" = "#4DAF4A",
  "Cluster 4" = "#984EA3"
)

cluster_labels <- MANUAL_CLUSTER_LABELS[
  as.character(
    cluster_ordered
  )
]

ha_top <- HeatmapAnnotation(
  Cluster = factor(
    cluster_labels,
    levels = c(
      "Cluster 1",
      "Cluster 4",
      "Cluster 3",
      "Cluster 2"
    )
  ),
  annotation_name_side = "left",
  annotation_height = unit(
    4,
    "mm"
  ),
  col = list(
    Cluster = cluster_colors
  ),
  simple_anno_size = unit(
    4,
    "mm"
  ),
  border = TRUE
)


################################################################################
# 16. Chromosome-arm annotation
################################################################################

chr_arm <- as.character(
  fdata[
    selected$Probe_ID,
    "CHR_ARM"
  ]
)

chr_arm[
  is.na(chr_arm) |
    chr_arm == ""
] <- "Unplaced"

chr_arm_colors <- c(
  "1q" = "#E41A1C",
  "1p" = "#32BE40",
  "2q" = "#2EB5CD",
  "2p" = "#FF7F00",
  "3p" = "#645A1E",
  "3q" = "#009E73",
  "4q" = "#D391C0",
  "4p" = "#168169",
  "5q" = "#C0A421",
  "5p" = "#91A4E3",
  "6q" = "#CF3174",
  "6p" = "#F9837E",
  "7q" = "#8B470E",
  "7p" = "#AF6AF5",
  "8q" = "#348211",
  "8p" = "#57B69C",
  "9q" = "#1879AB",
  "10q" = "#AF8351",
  "10p" = "#8E4244",
  "11q" = "#889151",
  "12q" = "#7F67A1",
  "13q" = "#EE43BD",
  "14q" = "#3F51B5",
  "15q" = "#306437",
  "16q" = "#D0724C",
  "16p" = "#86436A",
  "17q" = "#C67384",
  "18q" = "#359E60",
  "18p" = "#92B231",
  "19q" = "#B33DC3",
  "19p" = "#8F6D14",
  "20q" = "#60ADE2",
  "20p" = "#86B179",
  "22q" = "#7B3294",
  "Xq" = "#0B9C97",
  "Xp" = "#DB992B",
  "Unplaced" = "#BDBDBD"
)

ig <- as.character(
  fdata[
    selected$Probe_ID,
    "LOCUS_TYPE"
  ]
)

ig <- ifelse(
  tolower(
    ig
  ) == "immunoglobulin",
  "IG",
  "Other"
)

row_annotation <- rowAnnotation(
  CHR_ARM = factor(
    chr_arm,
    levels = names(
      chr_arm_colors
    )
  ),
  IG = factor(
    ig,
    levels = c(
      "IG",
      "Other"
    )
  ),
  col = list(
    CHR_ARM = chr_arm_colors,
    IG = c(
      "IG" = "darkorange",
      "Other" = "grey85"
    )
  ),
  na_col = "grey90",
  annotation_name_side = "top",
  width = unit(
    7,
    "mm"
  ),
  border = TRUE
)


################################################################################
# 17. Heatmap colors
################################################################################

color_function <- circlize::colorRamp2(
  c(
    -CLIP_Z,
    0,
    CLIP_Z
  ),
  c(
    "blue",
    "white",
    "red"
  )
)


################################################################################
# 18. Create Figure 1C
################################################################################

full_heatmap <- Heatmap(
  zmat_ordered,
  name = "Row Z",
  col = color_function,
  top_annotation = ha_top,
  right_annotation = row_annotation,
  show_column_names = FALSE,
  show_row_names = FALSE,
  column_split = cluster_ordered,
  cluster_columns = TRUE,
  cluster_column_slices = FALSE,
  cluster_rows = TRUE,
  column_title = paste0(
    "Figure 1C: NDMM consensus clusters (K = ",
    K_USE,
    ") — plasma-cell Cox-selected TTP genes\n",
    "PC_ME_Other == ",
    PC_ME_OTHER_KEEP,
    "; Cox FDR < ",
    FDR_CUTOFF,
    "; all qualifying plasma-cell probes"
  ),
  heatmap_legend_param = list(
    at = c(
      -CLIP_Z,
      0,
      CLIP_Z
    )
  )
)


################################################################################
# 19. Save Figure 1C
################################################################################

pdf_file <- file.path(
  HM_DIR,
  "Fig1C_ConsensusClustering_TTP_K4.pdf"
)

png_file <- file.path(
  HM_DIR,
  "Fig1C_ConsensusClustering_TTP_K4.png"
)

pdf(
  pdf_file,
  width = 14,
  height = 10
)

draw(
  full_heatmap,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

dev.off()

png(
  png_file,
  width = 2600,
  height = 1800,
  res = 200
)

draw(
  full_heatmap,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

dev.off()


################################################################################
# 20. Save probe list
################################################################################

write_csv(
  selected,
  file.path(
    OUT_DIR,
    "Fig1C_AllQualifyingPlasmaCellProbes.csv"
  )
)


################################################################################
# 21. Complete
################################################################################

cat(
  "\nStep 11 complete.\n"
)

cat(
  "Samples:",
  ncol(zmat),
  "\n"
)

cat(
  "Probes clustered:",
  nrow(zmat),
  "\n"
)

cat(
  "Clusters:",
  K_USE,
  "\n"
)

cat(
  "Figure:",
  png_file,
  "\n"
)