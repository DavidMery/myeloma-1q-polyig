#!/usr/bin/env bash

set -euo pipefail

cd "$(dirname "$0")/.."

Rscript code/scripts/00_prepare-environment.R
Rscript code/scripts/01_check-missing-values.R
Rscript code/scripts/02_prepare-uams-data.R
Rscript code/scripts/03_train-gep1q-elastic-net-cv.R
Rscript code/scripts/04_evaluate-gep1q-penalized-heldout-confusion-matrix.R
Rscript code/scripts/05_apply-gep1q-lasso-score.R
Rscript code/scripts/06_add-fish-data.R
Rscript code/scripts/07_train-and-apply-polyig-score.R
Rscript code/scripts/08_classify-pc-me-genes.R
Rscript code/scripts/09_fig1a-gep70-km.R
Rscript code/scripts/10_fig1b-cox-volcano.R
Rscript code/scripts/11_fig1c-consensus-clustering-heatmap.R
Rscript code/scripts/12_fig1d_cluster_enrichment.R
Rscript code/scripts/13_fish1q-ttp-km.R
Rscript code/scripts/14_unigh-ttp-km.R
Rscript code/scripts/15_fish1q-unigh-ttp-km.R
Rscript code/scripts/16_gep1q-ttp-km.R
Rscript code/scripts/17_polyig-ttp-km.R
Rscript code/scripts/18_fig1e-gep1q-polyig-km.R
Rscript code/scripts/19_add-uams-binary-pdata-columns.R
Rscript code/scripts/20_add-uams-composite-rnas-groups.R
Rscript code/scripts/21_match-rnas-groups-to-rnbx-by-patid.R
Rscript code/scripts/22_rnbx-mgus-vs-gep70-limma.R
Rscript code/scripts/23_fig2b_Microenvironment_Heatmap.R
Rscript code/scripts/24_fig2a_Microenvironment_Heatmap.R
Rscript code/scripts/25_fig3a-gep1q-polyig-km-tt2.R
Rscript code/scripts/26_fig3b-gep1q-polyig-km-tt3plus.R
Rscript code/scripts/27_Fig4A-D_4_6_8yr_Full_followup.R
Rscript code/scripts/28_fig5a-gep1q-polyig-km.R
Rscript code/scripts/29_fig5b-fish1q-unigh-ttp-km.R
Rscript code/scripts/30_overlay-gep-and-clinical-ttr-curves.R
Rscript code/scripts/31_clinical-gep-confusion-matrix.R
Rscript code/scripts/32_univariate-cox.R
Rscript code/scripts/33_multivariate-cox.R