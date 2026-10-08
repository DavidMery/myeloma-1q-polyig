################################################################################
# STEP 12 — FIGURE 1D: FOUR-QUADRANT CLUSTER ENRICHMENT
# Stripped-down version of 12b.
# Keeps the 12b figure appearance and analytical definitions.
# IMPORTANT:
#   Step 11 clusters ALL qualifying plasma-cell probes.
#   Step 12 then keeps ONE representative probe per gene for the Figure 1D
#   gene-level display, as in the original 12b script.
################################################################################

suppressPackageStartupMessages({
  library(Biobase)
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(limma)
  library(ggrepel)
})


################################################################################
# 0. VERIFY PIPELINE ROOTS
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
# 1. SETTINGS
################################################################################

ESET_PATH <- file.path(
  DATA_ROOT,
  "processed",
  "ESET_uams_all_gep1q_fish_polyig_pcme.rds"
)

FIGURE_DIR <- file.path(
  RESULTS_ROOT,
  "figures",
  "Fig1"
)

STEP11_DIR <- file.path(
  FIGURE_DIR,
  "Fig1C_CC"
)

OUT_DIR <- file.path(
  FIGURE_DIR,
  "Fig1D_Quadrant"
)

K_USE <- 4L
FDR_CUTOFF <- 0.01
PC_ME_OTHER_KEEP <- 1

MIN_BEST_LOG2FC <- 0
REQUIRE_BEST_FDR <- TRUE
BEST_FDR_CUTOFF <- 0.05

FILTER_CHR_ARMS_BY_MIN_GENES <- TRUE
MIN_GENES_PER_CHR_ARM <- 4L

# 12b quadrant placement
UPPER_LEFT_CLUSTER_ID  <- "1"
UPPER_RIGHT_CLUSTER_ID <- "2"
LOWER_LEFT_CLUSTER_ID  <- "3"
LOWER_RIGHT_CLUSTER_ID <- "4"

# 12b display labels
CLUSTER_LABELS <- c(
  "1" = "Cluster 1",
  "2" = "Cluster 2",
  "3" = "Cluster 3",
  "4" = "Cluster 4"
)


# 12b displayed cluster colors
CLUSTER_LABEL_COLORS <- c(
  "Cluster 1" = "firebrick2",
  "Cluster 2" = "steelblue",
  "Cluster 3" = "green4",
  "Cluster 4" = "purple4"
)

# 12b displayed cluster colors
#CLUSTER_LABEL_COLORS <- c(
#  "Cluster 1" = "#33A02C",
#  "Cluster 2" = "#A6CEE3",
#  "Cluster 3" = "#B2DF8A",
#  "Cluster 4" = "#1F78B4"
#)

# 12b plot appearance
POINT_SIZE <- 1.8
X_AXIS_LIMIT_LOG2FC <- 3.0
Y_AXIS_LIMIT_LOG2FC <- 3.0
AXIS_TICK_STEP_LOG2FC <- 0.5

PLOT_BASE_FONT_SIZE <- 18
PLOT_TITLE_SIZE <- 22
PLOT_SUBTITLE_SIZE <- 16
AXIS_TITLE_SIZE <- 18
AXIS_TEXT_SIZE <- 15
LEGEND_TITLE_SIZE <- 17
LEGEND_TEXT_SIZE <- 14
QUADRANT_LABEL_SIZE <- 7
QUADRANT_LABEL_FACE <- "bold"
LEGEND_SYMBOL_SIZE <- 5

UPPER_LEFT_LABEL_X_MULT <- -0.5
UPPER_LEFT_LABEL_Y_MULT <-  0.98
UPPER_RIGHT_LABEL_X_MULT <-  0.5
UPPER_RIGHT_LABEL_Y_MULT <-  0.98
LOWER_LEFT_LABEL_X_MULT <- -0.5
LOWER_LEFT_LABEL_Y_MULT <- -0.84
LOWER_RIGHT_LABEL_X_MULT <-  0.5
LOWER_RIGHT_LABEL_Y_MULT <- -0.84

PLOT_WIDTH <- 10
PLOT_HEIGHT <- 8
PLOT_DPI <- 600

# 12b chromosome-arm colors
CHR_ARM_COLORS <- c(
  "1q"="#E41A1C", "1p"="#32BE40", "2q"="#2EB5CD", "2p"="#FF7F00",
  "3p"="#645A1E", "3q"="#009E73", "4q"="#D391C0", "4p"="#168169",
  "5q"="#C0A421", "5p"="#91A4E3", "6q"="#CF3174", "6p"="#F9837E",
  "7q"="#8B470E", "7p"="#AF6AF5", "8q"="#348211", "8p"="#57B69C",
  "9q"="#1879AB", "10q"="#AF8351", "10p"="#8E4244", "11q"="#889151",
  "12q"="#7F67A1", "13q"="#EE43BD", "14q"="#3F51B5", "15q"="#306437",
  "16q"="#D0724C", "16p"="#86436A", "17q"="#C67384", "18q"="#359E60",
  "18p"="#92B231", "19q"="#B33DC3", "19p"="#8F6D14", "20q"="#60ADE2",
  "20p"="#86B179", "22q"="#7B3294", "Xq"="#0B9C97", "Xp"="#DB992B",
  "Unplaced"="#BDBDBD"
)

CHR_ARM_ORDER <- c(
  "1q","1p","2q","2p","3q","3p","4q","4p","5q","5p",
  "6q","6p","7q","7p","8q","8p","9q","9p","10q","10p",
  "11q","11p","12q","12p","13q","13p","14q","14p","15q","15p",
  "16q","16p","17q","17p","18q","18p","19q","19p","20q","20p",
  "21q","21p","22q","22p","Xq","Xp","Yq","Yp","Unplaced"
)

HIGHLIGHT_ONLY_SELECTED_CHR_ARMS <- FALSE
HIGHLIGHT_CHR_ARMS <- c("1q", "2p", "14q", "22q")
OTHER_CHR_ARM_COLOR <- "grey85"
CHR_ARM_LEGEND_NCOL <- 3

# Centroid appearance from 12b
MAKE_CENTROID_PLOT <- TRUE
CONGLOMERATE_INCLUDE_UNPLACED <- FALSE
CONGLOMERATE_MAX_POINT_SIZE <- 20
CONGLOMERATE_LABEL_MODE <- "highlighted_only"
CONGLOMERATE_LABEL_SIZE <- 4.8
EMPHASIZED_CENTROID_LABEL_SIZE <- 5.5
EMPHASIZED_CENTROID_LABEL_SEGMENT_SIZE <- 0.7

FORCED_CENTROID_LABEL_POSITIONS <- tibble::tribble(
  ~CHR_ARM, ~LabelX, ~LabelY, ~PointClearance, ~TextClearance, ~HJust, ~VJust,
  "1q",    1.05,  -0.55,  0.1, 0.1, 0.0, 0.5,
  "2p",   -1.05,   0.45,  0.1, 0.1, 1.0, 0.5,
  "14q",  -1.65,   1.75,  0.1, 0.1, 0.5, 0.5,
  "22q",  -1.65,   2.25,  0.1, 0.1, 0.0, 0.5
)

EMPHASIZED_CENTROID_LABEL_ARMS <- c(
  "1q",
  "2p",
  "14q",
  "22q"
)


################################################################################
# 2. FILES
################################################################################

CLUSTER_FILE <- file.path(
  STEP11_DIR,
  "Fig1C_ConsensusClusters_TTP_K4_Assignments.csv"
)

PROBE_FILE <- file.path(
  STEP11_DIR,
  "Fig1C_AllQualifyingPlasmaCellProbes.csv"
)

if (!file.exists(ESET_PATH)) {
  stop(
    "Processed ExpressionSet not found: ",
    ESET_PATH
  )
}

if (!file.exists(CLUSTER_FILE)) {
  stop(
    "Step 11 cluster file not found: ",
    CLUSTER_FILE
  )
}

if (!file.exists(PROBE_FILE)) {
  stop(
    "Step 11 probe file not found: ",
    PROBE_FILE
  )
}

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


################################################################################
# 3. LOAD ESET AND STEP-11 OUTPUTS
################################################################################

eset <- readRDS(
  ESET_PATH
)

expr_all <- Biobase::exprs(
  eset
)

fdata <- as.data.frame(
  Biobase::fData(
    eset
  ),
  check.names = FALSE
)

clusters <- read_csv(
  CLUSTER_FILE,
  show_col_types = FALSE
) |>
  mutate(
    CHIPID = as.character(CHIPID),
    Cluster = as.character(Cluster)
  )

selected <- read_csv(
  PROBE_FILE,
  show_col_types = FALSE
)


################################################################################
# 4. PREPARE FEATURE ANNOTATION
################################################################################

selected <- selected |>
  mutate(
    FeatureKey = as.character(Probe_ID),
    Probe_ID = as.character(Probe_ID)
  )

# Use the ExpressionSet fData as the authoritative source for annotation.
selected <- selected |>
  select(
    -any_of(
      c(
        "GENE",
        "CHR_ARM",
        "PC_ME_Other"
      )
    )
  ) |>
  left_join(
    fdata |>
      rownames_to_column(
        "FeatureKey"
      ) |>
      transmute(
        FeatureKey,
        GENE = as.character(GENE),
        CHR_ARM = as.character(CHR_ARM),
        PC_ME_Other = as.numeric(PC_ME_Other)
      ),
    by = "FeatureKey"
  ) |>
  filter(
    PC_ME_Other == PC_ME_OTHER_KEEP
  )

selected$CHR_ARM <- gsub(
  "^chr",
  "",
  selected$CHR_ARM,
  ignore.case = TRUE
)

selected$CHR_ARM[
  is.na(selected$CHR_ARM) |
    selected$CHR_ARM == ""
] <- "Unplaced"


################################################################################
# 5. ALIGN EXPRESSION TO CLUSTERED SAMPLES
################################################################################

sample_index <- match(
  clusters$CHIPID,
  colnames(expr_all)
)

if (anyNA(sample_index)) {
  stop(
    "Some Step-11 clustered samples are missing from the ExpressionSet."
  )
}

expr_selected <- expr_all[
  selected$FeatureKey,
  sample_index,
  drop = FALSE
]

colnames(
  expr_selected
) <- clusters$CHIPID

cluster_factor <- factor(
  clusters$Cluster,
  levels = c(
    "1",
    "2",
    "3",
    "4"
  )
)


################################################################################
# 6. LIMMA: EACH CLUSTER VS THE OTHER THREE
################################################################################

design <- model.matrix(
  ~ 0 + cluster_factor
)

colnames(
  design
) <- paste0(
  "C",
  1:4
)

rownames(
  design
) <- colnames(
  expr_selected
)

fit <- limma::lmFit(
  expr_selected,
  design
)

for (k in 1:4) {
  
  others <- setdiff(
    1:4,
    k
  )
  
  contrast <- limma::makeContrasts(
    contrasts = paste0(
      "C",
      k,
      " - (",
      paste0(
        "C",
        others,
        collapse = " + "
      ),
      ")/3"
    ),
    levels = design
  )
  
  fit_k <- limma::eBayes(
    limma::contrasts.fit(
      fit,
      contrast
    )
  )
  
  tmp <- limma::topTable(
    fit_k,
    number = Inf,
    sort.by = "none"
  )
  
  selected[[
    paste0(
      "log2FC_C",
      k
    )
  ]] <- tmp$logFC
  
  selected[[
    paste0(
      "PValue_C",
      k
    )
  ]] <- tmp$P.Value
  
  selected[[
    paste0(
      "FDR_C",
      k
    )
  ]] <- tmp$adj.P.Val
}


################################################################################
# 7. CLUSTER MEANS AND BEST CLUSTER
################################################################################

for (k in 1:4) {
  
  selected[[
    paste0(
      "MeanLog2Expr_C",
      k
    )
  ]] <- rowMeans(
    expr_selected[
      ,
      cluster_factor == k,
      drop = FALSE
    ]
  )
}

mean_mat <- as.matrix(
  selected[
    ,
    paste0(
      "MeanLog2Expr_C",
      1:4
    )
  ]
)

best <- max.col(
  mean_mat,
  ties.method = "first"
)

selected$BestCluster <- as.character(
  best
)

selected$BestCluster_Label <- CLUSTER_LABELS[
  selected$BestCluster
]

fdr_mat <- as.matrix(
  selected[
    ,
    paste0(
      "FDR_C",
      1:4
    )
  ]
)

logfc_mat <- as.matrix(
  selected[
    ,
    paste0(
      "log2FC_C",
      1:4
    )
  ]
)

selected$Best_log2FC_vsRest <- logfc_mat[
  cbind(
    seq_len(
      nrow(selected)
    ),
    best
  )
]

selected$Best_FDR <- fdr_mat[
  cbind(
    seq_len(
      nrow(selected)
    ),
    best
  )
]

selected$Best_FoldChange_vsRest <- 2^selected$Best_log2FC_vsRest


################################################################################
# 8. ONE PROBE PER GENE FOR FIGURE 1D
################################################################################

quadrant_tbl <- selected |>
  mutate(
    GENE = as.character(GENE),
    GENE_KEY = toupper(GENE)
  ) |>
  filter(
    !is.na(GENE),
    GENE != "",
    tolower(GENE) != "no-gene",
    is.finite(Best_log2FC_vsRest),
    Best_log2FC_vsRest > MIN_BEST_LOG2FC
  )

if (REQUIRE_BEST_FDR) {
  
  quadrant_tbl <- quadrant_tbl |>
    filter(
      is.finite(Best_FDR),
      Best_FDR < BEST_FDR_CUTOFF
    )
}

quadrant_tbl <- quadrant_tbl |>
  arrange(
    desc(Best_log2FC_vsRest),
    Best_FDR,
    Probe_ID
  ) |>
  group_by(
    GENE_KEY
  ) |>
  slice_head(
    n = 1
  ) |>
  ungroup()


################################################################################
# 9. REMOVE CHR_ARM GROUPS WITH <5 GENES
################################################################################

quadrant_tbl$CHR_ARM <- sapply(
  quadrant_tbl$CHR_ARM,
  function(x) as.character(x[[1]])
)

if (FILTER_CHR_ARMS_BY_MIN_GENES) {
  
  chr_counts <- table(
    quadrant_tbl$CHR_ARM
  )
  
  keep_chr <- names(
    chr_counts[
      chr_counts >= MIN_GENES_PER_CHR_ARM
    ]
  )
  
  quadrant_tbl <- quadrant_tbl[
    quadrant_tbl$CHR_ARM %in% keep_chr,
    ,
    drop = FALSE
  ]
}


################################################################################
# 10. FOUR-QUADRANT COORDINATES
################################################################################

quadrant_tbl <- quadrant_tbl |>
  mutate(
    Horizontal_log2FC = case_when(
      
      BestCluster == "1" ~
        MeanLog2Expr_C1 -
        (
          (
            MeanLog2Expr_C2 +
              MeanLog2Expr_C4
          ) / 2
        ),
      
      BestCluster == "2" ~
        MeanLog2Expr_C2 -
        (
          (
            MeanLog2Expr_C1 +
              MeanLog2Expr_C3
          ) / 2
        ),
      
      BestCluster == "3" ~
        MeanLog2Expr_C3 -
        (
          (
            MeanLog2Expr_C2 +
              MeanLog2Expr_C4
          ) / 2
        ),
      
      BestCluster == "4" ~
        MeanLog2Expr_C4 -
        (
          (
            MeanLog2Expr_C1 +
              MeanLog2Expr_C3
          ) / 2
        )
    ),
    
    Vertical_log2FC = case_when(
      
      BestCluster == "1" ~
        MeanLog2Expr_C1 -
        (
          (
            MeanLog2Expr_C3 +
              MeanLog2Expr_C4
          ) / 2
        ),
      
      BestCluster == "2" ~
        MeanLog2Expr_C2 -
        (
          (
            MeanLog2Expr_C3 +
              MeanLog2Expr_C4
          ) / 2
        ),
      
      BestCluster == "3" ~
        MeanLog2Expr_C3 -
        (
          (
            MeanLog2Expr_C1 +
              MeanLog2Expr_C2
          ) / 2
        ),
      
      BestCluster == "4" ~
        MeanLog2Expr_C4 -
        (
          (
            MeanLog2Expr_C1 +
              MeanLog2Expr_C2
          ) / 2
        )
    ),
    
    Horizontal_log2FC = pmax(
      Horizontal_log2FC,
      0
    ),
    
    Vertical_log2FC = pmax(
      Vertical_log2FC,
      0
    ),
    
    PlotX = if_else(
      BestCluster %in% c(
        UPPER_LEFT_CLUSTER_ID,
        LOWER_LEFT_CLUSTER_ID
      ),
      -Horizontal_log2FC,
      Horizontal_log2FC
    ),
    
    PlotY = if_else(
      BestCluster %in% c(
        LOWER_LEFT_CLUSTER_ID,
        LOWER_RIGHT_CLUSTER_ID
      ),
      -Vertical_log2FC,
      Vertical_log2FC
    )
  )


################################################################################
# 11. CHR_ARM COLORS
################################################################################

observed_arms <- intersect(
  CHR_ARM_ORDER,
  unique(
    quadrant_tbl$CHR_ARM
  )
)

if (HIGHLIGHT_ONLY_SELECTED_CHR_ARMS) {
  
  quadrant_tbl$CHR_ARM_PLOT <- ifelse(
    quadrant_tbl$CHR_ARM %in%
      HIGHLIGHT_CHR_ARMS,
    quadrant_tbl$CHR_ARM,
    "Other"
  )
  
  chr_colors <- c(
    CHR_ARM_COLORS[
      HIGHLIGHT_CHR_ARMS
    ],
    Other = OTHER_CHR_ARM_COLOR
  )
  
} else {
  
  quadrant_tbl$CHR_ARM_PLOT <- quadrant_tbl$CHR_ARM
  
  chr_colors <- CHR_ARM_COLORS[
    observed_arms
  ]
}

quadrant_tbl$CHR_ARM_PLOT <- factor(
  quadrant_tbl$CHR_ARM_PLOT,
  levels = names(
    chr_colors
  )
)


################################################################################
# 12. AXES AND CLUSTER LABELS
################################################################################

axis_breaks <- seq(
  -X_AXIS_LIMIT_LOG2FC,
  X_AXIS_LIMIT_LOG2FC,
  by = AXIS_TICK_STEP_LOG2FC
)

quadrant_counts <- table(
  factor(
    quadrant_tbl$BestCluster,
    levels = c(
      "1",
      "2",
      "3",
      "4"
    )
  )
)


################################################################################
# 13. INDIVIDUAL-GENE FIGURE 1D
################################################################################

p <- ggplot(
  quadrant_tbl,
  aes(
    PlotX,
    PlotY
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.9,
    color = "grey"
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.9,
    color = "grey"
  ) +
  geom_point(
    aes(
      color = CHR_ARM_PLOT
    ),
    shape = 16,
    size = POINT_SIZE,
    alpha = 1
  ) +
  annotate(
    "text",
    x = -X_AXIS_LIMIT_LOG2FC * 0.5,
    y = Y_AXIS_LIMIT_LOG2FC * 0.98,
    label = paste0(
      CLUSTER_LABELS[[
        UPPER_LEFT_CLUSTER_ID
      ]],
      "(n = ",
      quadrant_counts[[
        UPPER_LEFT_CLUSTER_ID
      ]],
      ")"
    ),
    color = CLUSTER_LABEL_COLORS[[
      CLUSTER_LABELS[[
        UPPER_LEFT_CLUSTER_ID
      ]]
    ]],
    size = QUADRANT_LABEL_SIZE,
    fontface = QUADRANT_LABEL_FACE
  ) +
  annotate(
    "text",
    x = X_AXIS_LIMIT_LOG2FC * 0.5,
    y = Y_AXIS_LIMIT_LOG2FC * 0.98,
    label = paste0(
      CLUSTER_LABELS[[
        UPPER_RIGHT_CLUSTER_ID
      ]],
      "(n = ",
      quadrant_counts[[
        UPPER_RIGHT_CLUSTER_ID
      ]],
      ")"
    ),
    color = CLUSTER_LABEL_COLORS[[
      CLUSTER_LABELS[[
        UPPER_RIGHT_CLUSTER_ID
      ]]
    ]],
    size = QUADRANT_LABEL_SIZE,
    fontface = QUADRANT_LABEL_FACE
  ) +
  annotate(
    "text",
    x = -X_AXIS_LIMIT_LOG2FC * 0.5,
    y = -Y_AXIS_LIMIT_LOG2FC * 0.84,
    label = paste0(
      CLUSTER_LABELS[[
        LOWER_LEFT_CLUSTER_ID
      ]],
      "(n = ",
      quadrant_counts[[
        LOWER_LEFT_CLUSTER_ID
      ]],
      ")"
    ),
    color = CLUSTER_LABEL_COLORS[[
      CLUSTER_LABELS[[
        LOWER_LEFT_CLUSTER_ID
      ]]
    ]],
    size = QUADRANT_LABEL_SIZE,
    fontface = QUADRANT_LABEL_FACE
  ) +
  annotate(
    "text",
    x = X_AXIS_LIMIT_LOG2FC * 0.5,
    y = -Y_AXIS_LIMIT_LOG2FC * 0.84,
    label = paste0(
      CLUSTER_LABELS[[
        LOWER_RIGHT_CLUSTER_ID
      ]],
      "(n = ",
      quadrant_counts[[
        LOWER_RIGHT_CLUSTER_ID
      ]],
      ")"
    ),
    color = CLUSTER_LABEL_COLORS[[
      CLUSTER_LABELS[[
        LOWER_RIGHT_CLUSTER_ID
      ]]
    ]],
    size = QUADRANT_LABEL_SIZE,
    fontface = QUADRANT_LABEL_FACE
  ) +
  scale_color_manual(
    values = chr_colors,
    limits = names(
      chr_colors
    ),
    breaks = names(
      chr_colors
    ),
    drop = FALSE
  ) +
  scale_x_continuous(
    breaks = axis_breaks
  ) +
  scale_y_continuous(
    breaks = axis_breaks
  ) +
  coord_equal(
    xlim = c(
      -X_AXIS_LIMIT_LOG2FC,
      X_AXIS_LIMIT_LOG2FC
    ),
    ylim = c(
      -Y_AXIS_LIMIT_LOG2FC,
      Y_AXIS_LIMIT_LOG2FC
    ),
    clip = "off"
  ) +
  labs(
    title = paste0(
      "Cluster-specific enrichment\n",
      "of plasma-cell Cox-selected genes (manually relabeled clusters)"
    ),
    subtitle = paste0(
      "Each solid point is one plasma-cell gene (PC_ME_Other == ",
      PC_ME_OTHER_KEEP,
      "). Coordinates are directional log2 fold-change values."
    ),
    x = "Directional enrichment (log2FC)",
    y = "Directional enrichment (log2FC)",
    color = "CHR_ARM"
  ) +
  theme_classic(
    base_size = PLOT_BASE_FONT_SIZE
  ) +
  theme(
    plot.title = element_text(
      size = PLOT_TITLE_SIZE,
      face = "bold",
      lineheight = 1.05,
      margin = margin(
        b = 6
      )
    ),
    plot.subtitle = element_text(
      size = PLOT_SUBTITLE_SIZE,
      lineheight = 1.10,
      margin = margin(
        b = 10
      )
    ),
    axis.title = element_text(
      size = AXIS_TITLE_SIZE
    ),
    axis.text = element_text(
      size = AXIS_TEXT_SIZE
    ),
    legend.title = element_text(
      size = LEGEND_TITLE_SIZE,
      face = "bold"
    ),
    legend.text = element_text(
      size = LEGEND_TEXT_SIZE
    ),
    legend.key = element_blank(),
    legend.spacing.y = grid::unit(
      4,
      "pt"
    ),
    plot.margin = margin(
      12,
      38,
      12,
      12
    )
  ) +
  guides(
    color = guide_legend(
      ncol = CHR_ARM_LEGEND_NCOL,
      byrow = TRUE,
      override.aes = list(
        shape = 16,
        size = LEGEND_SYMBOL_SIZE,
        alpha = 1
      )
    )
  )


################################################################################
# 14. SAVE FIGURE 1D
################################################################################

png_path <- file.path(
  OUT_DIR,
  "Fig1D_ClusterEnrichment_ALL_COX_SELECTED_K4_PC1_QuadrantLabels.png"
)

pdf_path <- file.path(
  OUT_DIR,
  "Fig1D_ClusterEnrichment_ALL_COX_SELECTED_K4_PC1_QuadrantLabels.pdf"
)

ggsave(
  png_path,
  p,
  width = PLOT_WIDTH,
  height = PLOT_HEIGHT,
  units = "in",
  dpi = PLOT_DPI,
  bg = "white"
)

ggsave(
  pdf_path,
  p,
  width = PLOT_WIDTH,
  height = PLOT_HEIGHT,
  units = "in",
  bg = "white"
)


################################################################################
# 15. CHR_ARM CENTROID FIGURE
################################################################################

if (MAKE_CENTROID_PLOT) {
  
  centroid <- quadrant_tbl |>
    filter(
      is.finite(PlotX),
      is.finite(PlotY),
      is.finite(Best_log2FC_vsRest)
    ) |>
    filter(
      if (CONGLOMERATE_INCLUDE_UNPLACED) {
        TRUE
      } else {
        CHR_ARM != "Unplaced"
      }
    ) |>
    group_by(
      CHR_ARM,
      CHR_ARM_PLOT
    ) |>
    summarise(
      n_genes = n(),
      Mean_PlotX = mean(PlotX),
      Mean_PlotY = mean(PlotY),
      .groups = "drop"
    )
  
  label_centroid <- centroid |>
    filter(
      CHR_ARM %in% HIGHLIGHT_CHR_ARMS
    )
  
  centroid_plot <- ggplot(
    centroid,
    aes(
      Mean_PlotX,
      Mean_PlotY
    )
  ) +
    geom_hline(
      yintercept = 0,
      linewidth = 0.9,
      color = "black"
    ) +
    geom_vline(
      xintercept = 0,
      linewidth = 0.9,
      color = "black"
    ) +
    geom_point(
      aes(
        fill = CHR_ARM_PLOT,
        size = n_genes
      ),
      shape = 21,
      color = "black",
      stroke = 0.45,
      alpha = 1
    ) +
    annotate(
      "text",
      x = -X_AXIS_LIMIT_LOG2FC * 0.5,
      y = Y_AXIS_LIMIT_LOG2FC * 0.98,
      label = CLUSTER_LABELS[[
        UPPER_LEFT_CLUSTER_ID
      ]],
      color = CLUSTER_LABEL_COLORS[[
        CLUSTER_LABELS[[
          UPPER_LEFT_CLUSTER_ID
        ]]
      ]],
      size = QUADRANT_LABEL_SIZE,
      fontface = QUADRANT_LABEL_FACE
    ) +
    annotate(
      "text",
      x = X_AXIS_LIMIT_LOG2FC * 0.5,
      y = Y_AXIS_LIMIT_LOG2FC * 0.98,
      label = CLUSTER_LABELS[[
        UPPER_RIGHT_CLUSTER_ID
      ]],
      color = CLUSTER_LABEL_COLORS[[
        CLUSTER_LABELS[[
          UPPER_RIGHT_CLUSTER_ID
        ]]
      ]],
      size = QUADRANT_LABEL_SIZE,
      fontface = QUADRANT_LABEL_FACE
    ) +
    annotate(
      "text",
      x = -X_AXIS_LIMIT_LOG2FC * 0.5,
      y = -Y_AXIS_LIMIT_LOG2FC * 0.84,
      label = CLUSTER_LABELS[[
        LOWER_LEFT_CLUSTER_ID
      ]],
      color = CLUSTER_LABEL_COLORS[[
        CLUSTER_LABELS[[
          LOWER_LEFT_CLUSTER_ID
        ]]
      ]],
      size = QUADRANT_LABEL_SIZE,
      fontface = QUADRANT_LABEL_FACE
    ) +
    annotate(
      "text",
      x = X_AXIS_LIMIT_LOG2FC * 0.5,
      y = -Y_AXIS_LIMIT_LOG2FC * 0.84,
      label = CLUSTER_LABELS[[
        LOWER_RIGHT_CLUSTER_ID
      ]],
      color = CLUSTER_LABEL_COLORS[[
        CLUSTER_LABELS[[
          LOWER_RIGHT_CLUSTER_ID
        ]]
      ]],
      size = QUADRANT_LABEL_SIZE,
      fontface = QUADRANT_LABEL_FACE
    ) +
    scale_fill_manual(
      values = chr_colors,
      limits = names(
        chr_colors
      ),
      breaks = names(
        chr_colors
      ),
      drop = FALSE
    ) +
    scale_size_area(
      max_size = CONGLOMERATE_MAX_POINT_SIZE,
      name = "Number of genes"
    ) +
    scale_x_continuous(
      breaks = axis_breaks
    ) +
    scale_y_continuous(
      breaks = axis_breaks
    ) +
    coord_equal(
      xlim = c(
        -X_AXIS_LIMIT_LOG2FC,
        X_AXIS_LIMIT_LOG2FC
      ),
      ylim = c(
        -Y_AXIS_LIMIT_LOG2FC,
        Y_AXIS_LIMIT_LOG2FC
      ),
      clip = "off"
    ) +
    labs(
      title = paste0(
        "Chromosome-arm mean enrichment\n",
        "across plasma-cell Cox-selected genes"
      ),
      subtitle = paste0(
        "Each dot is the mean signed position of plasma-cell genes on that\n",
        "chromosome arm. Dot area represents the number of genes."
      ),
      x = "Mean directional enrichment (log2FC)",
      y = "Mean directional enrichment (log2FC)",
      fill = "CHR_ARM"
    ) +
    theme_classic(
      base_size = PLOT_BASE_FONT_SIZE
    ) +
    theme(
      plot.title = element_text(
        size = PLOT_TITLE_SIZE,
        face = "bold",
        lineheight = 1.05,
        margin = margin(
          b = 6
        )
      ),
      plot.subtitle = element_text(
        size = PLOT_SUBTITLE_SIZE,
        lineheight = 1.10,
        margin = margin(
          b = 10
        )
      ),
      axis.title = element_text(
        size = AXIS_TITLE_SIZE
      ),
      axis.text = element_text(
        size = AXIS_TEXT_SIZE
      ),
      legend.title = element_text(
        size = LEGEND_TITLE_SIZE,
        face = "bold"
      ),
      legend.text = element_text(
        size = LEGEND_TEXT_SIZE
      ),
      legend.key = element_blank(),
      legend.spacing.y = grid::unit(
        4,
        "pt"
      ),
      plot.margin = margin(
        12,
        42,
        12,
        12
      )
    ) +
    guides(
      fill = guide_legend(
        ncol = CHR_ARM_LEGEND_NCOL,
        byrow = TRUE,
        override.aes = list(
          shape = 21,
          size = LEGEND_SYMBOL_SIZE,
          color = "black"
        )
      ),
      size = guide_legend(
        override.aes = list(
          fill = "grey75",
          color = "black",
          shape = 21
        )
      )
    )
  
  
  ##############################################################################
  # 15A. EMPHASIZED CHR_ARM LABELS WITH FORCED POSITIONS
  ##############################################################################
  
  forced <- label_centroid |>
    inner_join(
      FORCED_CENTROID_LABEL_POSITIONS,
      by = "CHR_ARM"
    ) |>
    mutate(
      dx = LabelX - Mean_PlotX,
      dy = LabelY - Mean_PlotY,
      d = sqrt(
        dx^2 +
          dy^2
      ),
      d = ifelse(
        is.finite(d) &
          d > 0,
        d,
        1
      ),
      SegmentX =
        Mean_PlotX +
        PointClearance *
        dx / d,
      SegmentY =
        Mean_PlotY +
        PointClearance *
        dy / d,
      SegmentXend =
        LabelX -
        TextClearance *
        dx / d,
      SegmentYend =
        LabelY -
        TextClearance *
        dy / d
    )
  
  if (nrow(forced) > 0) {
    
    centroid_plot <- centroid_plot +
      geom_segment(
        data = forced,
        aes(
          x = SegmentX,
          y = SegmentY,
          xend = SegmentXend,
          yend = SegmentYend
        ),
        inherit.aes = FALSE,
        linewidth = EMPHASIZED_CENTROID_LABEL_SEGMENT_SIZE,
        lineend = "round",
        color = "black"
      ) +
      geom_label(
        data = forced,
        aes(
          x = LabelX,
          y = LabelY,
          label = paste0(
            CHR_ARM,
            " (n=",
            n_genes,
            ")"
          ),
          hjust = HJust,
          vjust = VJust
        ),
        inherit.aes = FALSE,
        color = "black",
        fill = "white",
        alpha = 0.97,
        size = EMPHASIZED_CENTROID_LABEL_SIZE,
        fontface = "bold",
        label.size = 0,
        label.padding = grid::unit(
          0.08,
          "lines"
        ),
        lineheight = 0.95,
        show.legend = FALSE
      )
  }
  
  ordinary <- label_centroid |>
    filter(
      !(
        CHR_ARM %in%
          forced$CHR_ARM
      )
    ) |>
    filter(
      !(
        CHR_ARM %in%
          EMPHASIZED_CENTROID_LABEL_ARMS
      )
    )
  
  if (nrow(ordinary) > 0) {
    
    centroid_plot <- centroid_plot +
      geom_text_repel(
        data = ordinary,
        aes(
          label = paste0(
            CHR_ARM,
            " (n=",
            n_genes,
            ")"
          )
        ),
        inherit.aes = FALSE,
        x = ordinary$Mean_PlotX,
        y = ordinary$Mean_PlotY,
        color = "black",
        size = CONGLOMERATE_LABEL_SIZE,
        seed = 123,
        max.overlaps = Inf,
        box.padding = 0.55,
        point.padding = 0.35,
        min.segment.length = 0,
        force = 1.2,
        segment.size = 0.45
      )
  }
  
  
  ##############################################################################
  # 15B. SAVE CHROMOSOME-ARM CENTROID FIGURE
  ##############################################################################
  
  centroid_png <- file.path(
    OUT_DIR,
    "Fig1D_Supporting_CHR_ARM_Centroids_ALL_COX_SELECTED_K4_PC1_QuadrantLabels.png"
  )
  
  centroid_pdf <- file.path(
    OUT_DIR,
    "Fig1D_Supporting_CHR_ARM_Centroids_ALL_COX_SELECTED_K4_PC1_QuadrantLabels.pdf"
  )
  
  ggsave(
    centroid_png,
    centroid_plot,
    width = PLOT_WIDTH,
    height = PLOT_HEIGHT,
    units = "in",
    dpi = PLOT_DPI,
    bg = "white"
  )
  
  ggsave(
    centroid_pdf,
    centroid_plot,
    width = PLOT_WIDTH,
    height = PLOT_HEIGHT,
    units = "in",
    bg = "white"
  )
}


################################################################################
# 16. SOURCE DATA
################################################################################

write_csv(
  quadrant_tbl,
  file.path(
    OUT_DIR,
    "Fig1D_ClusterEnrichment_ALL_COX_SELECTED_K4_PC1_source-data.csv"
  )
)

write_csv(
  tibble(
    CHR_ARM = names(
      chr_colors
    ),
    Color = unname(
      chr_colors
    )
  ),
  file.path(
    OUT_DIR,
    "Fig1D_ClusterEnrichment_ALL_COX_SELECTED_K4_CHR_ARM_colors.csv"
  )
)


################################################################################
# 17. COMPLETE
################################################################################

cat(
  "\nStep 12 complete.\n",
  "Samples: ",
  ncol(expr_selected),
  "\n",
  "Genes plotted: ",
  nrow(quadrant_tbl),
  "\n",
  "FDR cutoff: ",
  FDR_CUTOFF,
  "\n",
  "Minimum best log2FC: ",
  MIN_BEST_LOG2FC,
  "\n",
  "CHR_ARM minimum genes: ",
  MIN_GENES_PER_CHR_ARM,
  "\n",
  "Figure: ",
  png_path,
  "\n",
  sep = ""
)