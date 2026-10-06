################################################################################
# STEP 10 — FIGURE 1B
#
# Purpose:
#   1. Restrict to RNAS CD138+ NDMM treated on TT2-TT5.
#   2. Restrict to GEP70-standard-risk disease (GEP70 <= 0.66).
#   3. Fit one continuous Cox model per probe.
#   4. Report effects per 1 SD higher expression.
#   5. Restrict Figure 1B to plasma-cell probes (PC_ME_Other == 1).
#   6. Highlight immunoglobulin probes and chromosome 1q probes.
#   7. Use parallel processing for probe-wise Cox models.
#   8. Save the Cox results, source data, and Figure 1B.
################################################################################

library(survival)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(parallel)


################################################################################
# 1. Settings
################################################################################

ENDPOINT <- "TTP"

GEP70_CUTOFF <- 0.66

FDR_LABEL_THRESHOLD <- 0.01

N_CORES <- 12L

TT_KEEP <- c(
  "TT2_NoThal",
  "TT2_Thal",
  "TT3a",
  "TT3b",
  "TT4_S-TT3",
  "TT5"
)


################################################################################
# 2. Files
################################################################################

ESET_FILE <- file.path(
  "data",
  "processed",
  "uams",
  "ESET_uams_all_gep1q_fish_polyig_pcme.rds"
)

OUT_DIR <- file.path(
  "results",
  "figures",
  "Fig1",
  "Fig1B"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


################################################################################
# 3. Load ExpressionSet
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
# 4. Select endpoint
################################################################################

if (ENDPOINT == "TTP") {
  
  TIME_COL <- "YearsTTP"
  EVENT_COL <- "CensTTP"
  
} else if (ENDPOINT == "OS") {
  
  TIME_COL <- "yearsOS_post_tx1"
  EVENT_COL <- "censOS_post_tx1"
  
} else {
  
  stop(
    "ENDPOINT must be TTP or OS."
  )
}


################################################################################
# 5. Restrict to GEP70-standard-risk NDMM treated on TT2-TT5
################################################################################

keep <- phenotype$sample_group == "RNAS_CD138_NDMM" &
  phenotype$tt_on_for_tx1 %in% TT_KEEP &
  as.numeric(
    phenotype$GEP70
  ) <= GEP70_CUTOFF

keep[is.na(keep)] <- FALSE

phenotype_model <- phenotype[
  keep,
  ,
  drop = FALSE
]

expression_model <- expression[
  ,
  rownames(phenotype_model),
  drop = FALSE
]


################################################################################
# 6. Prepare survival data
################################################################################

phenotype_model$time <- as.numeric(
  phenotype_model[
    ,
    TIME_COL
  ]
)

phenotype_model$event <- as.numeric(
  phenotype_model[
    ,
    EVENT_COL
  ]
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

expression_model <- expression_model[
  ,
  rownames(phenotype_model),
  drop = FALSE
]

cat(
  "Samples:",
  ncol(expression_model),
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
# 7. Prepare feature annotations
################################################################################

gene <- trimws(
  as.character(
    features$GENE
  )
)

chr_arm <- trimws(
  as.character(
    features$CHR_ARM
  )
)

pc_me <- as.numeric(
  features$PC_ME_Other
)


################################################################################
# 8. Identify immunoglobulin probes
################################################################################

ig_column <- which(
  colnames(features) == "locus_type_n54675"
)[1]

if (!is.na(ig_column)) {
  
  locus_type <- trimws(
    as.character(
      features[
        ,
        ig_column
      ]
    )
  )
  
} else {
  
  locus_type <- rep(
    NA_character_,
    nrow(features)
  )
}


################################################################################
# 9. Keep probes with valid gene annotations
################################################################################

keep_probe <- !is.na(gene) &
  gene != "" &
  tolower(gene) != "no-gene"

expression_model <- expression_model[
  keep_probe,
  ,
  drop = FALSE
]

gene <- gene[
  keep_probe
]

chr_arm <- chr_arm[
  keep_probe
]

pc_me <- pc_me[
  keep_probe
]

locus_type <- locus_type[
  keep_probe
]


################################################################################
# 10. Fit one Cox model per probe
#
# Effect:
#   HR per 1 SD higher expression
################################################################################

fit_probe <- function(i) {
  
  x <- as.numeric(
    expression_model[
      i,
      ,
      drop = TRUE
    ]
  )
  
  keep <- is.finite(
    x
  )
  
  x <- x[
    keep
  ]
  
  time <- phenotype_model$time[
    keep
  ]
  
  event <- phenotype_model$event[
    keep
  ]
  
  sd_x <- sd(
    x
  )
  
  if (
    length(x) < 20 ||
    sum(event == 1) < 5 ||
    !is.finite(sd_x) ||
    sd_x == 0
  ) {
    
    return(
      c(
        beta = NA_real_,
        HR = NA_real_,
        HR_lower95 = NA_real_,
        HR_upper95 = NA_real_,
        p_value = NA_real_
      )
    )
  }
  
  expression_z <- (
    x - mean(x)
  ) / sd_x
  
  fit <- coxph(
    Surv(
      time,
      event
    ) ~ expression_z
  )
  
  fit_summary <- summary(
    fit
  )
  
  c(
    beta = unname(
      coef(fit)[1]
    ),
    
    HR = unname(
      exp(
        coef(fit)[1]
      )
    ),
    
    HR_lower95 = unname(
      exp(
        confint(fit)[1]
      )
    ),
    
    HR_upper95 = unname(
      exp(
        confint(fit)[2]
      )
    ),
    
    p_value = unname(
      fit_summary$coef[
        1,
        "Pr(>|z|)"
      ]
    )
  )
}


################################################################################
# 11. Run probe-wise Cox analysis in parallel
################################################################################

cl <- makeCluster(
  N_CORES
)

clusterExport(
  cl,
  c(
    "expression_model",
    "phenotype_model",
    "fit_probe"
  ),
  envir = environment()
)

clusterEvalQ(
  cl,
  library(survival)
)

cox_results <- parLapply(
  cl,
  seq_len(
    nrow(expression_model)
  ),
  fit_probe
)

stopCluster(
  cl
)

cox_results <- do.call(
  rbind,
  cox_results
)

cox_results <- as.data.frame(
  cox_results
)


################################################################################
# 12. Assemble results
################################################################################

results <- data.frame(
  Probe_ID = rownames(
    expression_model
  ),
  GENE = gene,
  CHR_ARM = chr_arm,
  LOCUS_TYPE = locus_type,
  PC_ME_Other = pc_me,
  cox_results,
  stringsAsFactors = FALSE
)

results$FDR <- p.adjust(
  results$p_value,
  method = "BH"
)

results$neg_log10_FDR <- -log10(
  pmax(
    results$FDR,
    .Machine$double.xmin
  )
)


################################################################################
# 13. Define immunoglobulin and chromosome 1q probes
################################################################################

results$is_IG <- !is.na(
  results$LOCUS_TYPE
) &
  tolower(
    results$LOCUS_TYPE
  ) == "immunoglobulin"

results$is_1q <- tolower(
  gsub(
    "[[:space:]_-]+",
    "",
    results$CHR_ARM
  )
) %in% c(
  "1q",
  "chr1q"
)

results$Probe_Class <- "Other"

results$Probe_Class[
  results$is_IG
] <- "Immunoglobulin"

results$Probe_Class[
  !results$is_IG &
    results$is_1q
] <- "Chromosome 1q"


################################################################################
# 14. Save complete Cox results
################################################################################

write.csv(
  results,
  file.path(
    OUT_DIR,
    paste0(
      "Fig1B_Cox_",
      ENDPOINT,
      "_all-probe-results.csv"
    )
  ),
  row.names = FALSE
)


################################################################################
# 15. Restrict Figure 1B to plasma-cell probes
################################################################################

plot_data <- results |>
  filter(
    PC_ME_Other == 1,
    is.finite(beta),
    is.finite(FDR)
  )

label_data <- plot_data |>
  filter(
    FDR < FDR_LABEL_THRESHOLD,
    is_IG | is_1q
  )


################################################################################
# 16. Save Figure 1B source data
################################################################################

write.csv(
  plot_data,
  file.path(
    OUT_DIR,
    paste0(
      "Fig1B_Cox_",
      ENDPOINT,
      "_source-data.csv"
    )
  ),
  row.names = FALSE
)


################################################################################
# 17. Create Figure 1B
################################################################################

p <- ggplot(
  plot_data,
  aes(
    x = beta,
    y = neg_log10_FDR
  )
) +
  geom_hline(
    yintercept = -log10(
      FDR_LABEL_THRESHOLD
    ),
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = 0
  ) +
  geom_point(
    data = plot_data |>
      filter(
        Probe_Class == "Other"
      ),
    color = "grey65",
    alpha = 0.4,
    size = 1.4
  ) +
  geom_point(
    data = plot_data |>
      filter(
        Probe_Class != "Other"
      ),
    aes(
      color = Probe_Class
    ),
    alpha = 0.7,
    size = 1.9
  ) +
  scale_color_manual(
    values = c(
      "Immunoglobulin" = "steelblue",
      "Chromosome 1q" = "firebrick"
    )
  ) +
  labs(
    title = ifelse(
      ENDPOINT == "TTP",
      "Genes Associated with Time to Progression in GEP70 Standard-Risk NDMM",
      "Genes Associated with Overall Survival in GEP70 Standard-Risk NDMM"
    ),
    subtitle = paste0(
      "Baseline CD138+ NDMM; n = ",
      ncol(expression_model),
      "; events = ",
      sum(
        phenotype_model$event == 1
      ),
      "; effect per 1 SD higher expression"
    ),
    x = "Cox coefficient: log(HR) per 1 SD higher expression",
    y = expression(
      -log[10]("FDR")
    ),
    color = "Probe class"
  ) +
  theme_classic(
    base_size = 13
  )


################################################################################
# 18. Add IG and 1q labels at FDR < 0.01
################################################################################

if (
  nrow(
    label_data
  ) > 0
) {
  
  p <- p +
    geom_text_repel(
      data = label_data,
      aes(
        label = GENE,
        color = Probe_Class
      ),
      size = 3,
      max.overlaps = Inf,
      show.legend = FALSE
    )
}


################################################################################
# 19. Save Figure 1B
################################################################################

print(
  p
)

ggsave(
  file.path(
    OUT_DIR,
    paste0(
      "Fig1B_Cox_",
      ENDPOINT,
      ".png"
    )
  ),
  p,
  width = 9,
  height = 6,
  dpi = 600
)

ggsave(
  file.path(
    OUT_DIR,
    paste0(
      "Fig1B_Cox_",
      ENDPOINT,
      ".pdf"
    )
  ),
  p,
  width = 9,
  height = 6
)


################################################################################
# 20. Complete
################################################################################

cat(
  "Cox cohort:",
  ncol(expression_model),
  "samples\n"
)

cat(
  "Events:",
  sum(
    phenotype_model$event == 1
  ),
  "\n"
)

cat(
  "Probes analyzed:",
  nrow(results),
  "\n"
)

cat(
  "Plasma-cell probes plotted:",
  nrow(plot_data),
  "\n"
)

cat(
  "IG/1q probes labeled at FDR <",
  FDR_LABEL_THRESHOLD,
  ":",
  nrow(label_data),
  "\n"
)