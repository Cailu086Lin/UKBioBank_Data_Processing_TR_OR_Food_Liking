# UKBioBank_Data_Processing_TR_OR_Food_Liking
Analysis provided in the manuscript titled "Genes promoting liking onions lower blood pressure and type 2 diabetes risk"

This repository is meant to provide the source code for the analysis provided in the above manuscript by Liang-Dar Hwang, Cailu Lin, David M Evans, Nicholas G Martin, Danielle R Reed, Paule V Joseph

# Abstract

Food preferences shape dietary habits and influence long-term health. While hundreds of loci have been associated with food intake, they often reflect an effect of disease on diet, rather than diet on risk of disease. Taste and olfactory perception influence food preferences and choices before the onset of disease. Here, we investigate the influence of nonsynonymous genetic variation within taste and olfactory receptor genes on food preferences in UK Biobank (mean age 57). We identify 700 associations, of which 84 are also associated with their corresponding food intake traits. We replicate 45 of these associations in the independent Avon Longitudinal Study of Parents and Children cohort, and show that the effect sizes are larger in this younger cohort (mean age 25). Using Mendelian randomization, we show that a genetic predisposition to liking onions lowers blood pressure and the risk of type 2 diabetes. These findings advance our understanding of the direct genetic influences on food preferences and their consequent impacts on eating behaviors and disease susceptibility. Furthermore, our study demonstrates how the genetics of chemosensory perception can be used to strengthen causal inference in nutritional epidemiology, informing public health strategies to prevent and manage conditions like hypertension and type 2 diabetes.

# sessionInfo()
R version 4.4.1 (2024-06-14) Platform: x86_64-pc-linux-gnu Running under: Ubuntu 22.04.5 LTS

Matrix products: default BLAS: /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so; LAPACK version 3.10.0

locale: [1] LC_CTYPE=en_US.UTF-8 LC_NUMERIC=C LC_TIME=en_US.UTF-8 LC_COLLATE=en_US.UTF-8
[5] LC_MONETARY=en_US.UTF-8 LC_MESSAGES=en_US.UTF-8 LC_PAPER=en_US.UTF-8 LC_NAME=C
[9] LC_ADDRESS=C LC_TELEPHONE=C LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

time zone: America/New_York tzcode source: system (glibc)

attached base packages: [1] stats graphics grDevices datasets utils methods base

other attached packages: [1] Rmpfr_1.0-0 gmp_0.7-5 readxl_1.4.3 PredictABEL_1.2-4 lm.beta_1.7-2 scales_1.3.0
[7] tidyr_1.3.1 dplyr_1.1.4 broom_1.0.7 vcfR_1.15.0

loaded via a namespace (and not attached): [1] gtable_0.3.6 xfun_0.51 ggplot2_3.5.1 htmlwidgets_1.6.4 lattice_0.22-6
[6] vctrs_0.6.5 tools_4.4.1 generics_0.1.3 parallel_4.4.1 tibble_3.2.1
[11] cluster_2.1.6 pacman_0.5.1 pkgconfig_2.0.3 Matrix_1.7-2 data.table_1.17.0
[16] checkmate_2.3.2 lifecycle_1.0.4 compiler_4.4.1 stringr_1.5.1 munsell_0.5.1
[21] permute_0.9-7 htmltools_0.5.8.1 brant_0.3-0 yaml_2.3.10 htmlTable_2.4.3
[26] Formula_1.2-5 pillar_1.10.1 MASS_7.3-61 vegan_2.6-10 sessioninfo_1.2.3
[31] Hmisc_5.2-2 rpart_4.1.23 nlme_3.1-167 tidyselect_1.2.1 digest_0.6.37
[36] stringi_1.8.4 purrr_1.0.4 splines_4.4.1 cowplot_1.1.3 fastmap_1.2.0
[41] grid_4.4.1 colorspace_2.1-1 cli_3.6.4 magrittr_2.0.3 base64enc_0.1-3
[46] XML_3.99-0.18 ape_5.8-1 ROCR_1.0-11 foreign_0.8-88 backports_1.5.0
[51] rmarkdown_2.29 nnet_7.3-19 gridExtra_2.3 cellranger_1.1.0 evaluate_1.0.3
[56] knitr_1.49 tcltk_4.4.1 viridisLite_0.4.2 mgcv_1.9-1 rlang_1.1.5
[61] Rcpp_1.0.14 pinfsc50_1.3.0 xtable_1.8-4 glue_1.8.0 renv_1.0.7
[66] rstudioapi_0.17.1 R6_2.6.1 PBSmodelling_2.69.3
