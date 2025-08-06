# Univariate Analysis Log

This file contains the console output from the analysis script run on: 2025-08-06 13:03:11 




<<<<< RUNNING ALL ANALYSES FOR NEGATIVE MODE >>>>>

==================================================
STARTING KRUSKAL-WALLIS ANALYSIS FOR MODE: NEG 
==================================================
--------------------------------------------------
Running Kruskal-Wallis Analysis: All Groups Comparison 
Successfully saved KW results for: All Groups Comparison 
Found 523 features from Kruskal-Wallis test. Generating outputs for all...
Finished generating individual Kruskal-Wallis feature outputs.
--------------------------------------------------
Running Kruskal-Wallis Analysis: Healthy vs TR vs PR at Timepoint 1 
Successfully saved KW results for: Healthy vs TR vs PR at Timepoint 1 
Found 523 features from Kruskal-Wallis test. Generating outputs for all...
Finished generating individual Kruskal-Wallis feature outputs.
--------------------------------------------------
Running Kruskal-Wallis Analysis: Healthy vs TR vs PR at Timepoint 2 
Successfully saved KW results for: Healthy vs TR vs PR at Timepoint 2 
Found 523 features from Kruskal-Wallis test. Generating outputs for all...
Finished generating individual Kruskal-Wallis feature outputs.

==================================================
STARTING WILCOXON ANALYSIS FOR MODE: NEG 
==================================================
--------------------------------------------------
Running Wilcoxon Analysis: Paired TR Timepoint 1 vs 2 
Successfully saved Wilcoxon results for: Paired TR Timepoint 1 vs 2 
Found 523 features from Wilcoxon test. Generating outputs for all...
Finished generating individual Wilcoxon feature outputs.
--------------------------------------------------
Running Wilcoxon Analysis: Paired PR Timepoint 1 vs 2 
Successfully saved Wilcoxon results for: Paired PR Timepoint 1 vs 2 
Found 523 features from Wilcoxon test. Generating outputs for all...
Finished generating individual Wilcoxon feature outputs.
--------------------------------------------------
Running Wilcoxon Analysis: Unpaired Healthy vs TR at Timepoint 1 
Successfully saved Wilcoxon results for: Unpaired Healthy vs TR at Timepoint 1 
Found 523 features from Wilcoxon test. Generating outputs for all...
Finished generating individual Wilcoxon feature outputs.

==================================================
STARTING INTERACTION ANALYSIS FOR MODE: NEG 
==================================================
--------------------------------------------------
Running Interaction Analysis: Interaction of Change (T2-T1) between TR and PR groups 
Successfully saved Interaction results for: Interaction of Change (T2-T1) between TR and PR groups 



<<<<< RUNNING ALL ANALYSES FOR POSITIVE MODE >>>>>

==================================================
STARTING KRUSKAL-WALLIS ANALYSIS FOR MODE: POS 
==================================================
--------------------------------------------------
Running Kruskal-Wallis Analysis: All Groups Comparison 
Successfully saved KW results for: All Groups Comparison 
Found 1005 features from Kruskal-Wallis test. Generating outputs for all...
Finished generating individual Kruskal-Wallis feature outputs.
--------------------------------------------------
Running Kruskal-Wallis Analysis: Healthy vs TR vs PR at Timepoint 1 
Successfully saved KW results for: Healthy vs TR vs PR at Timepoint 1 
Found 1005 features from Kruskal-Wallis test. Generating outputs for all...
Finished generating individual Kruskal-Wallis feature outputs.
--------------------------------------------------
Running Kruskal-Wallis Analysis: Healthy vs TR vs PR at Timepoint 2 
Successfully saved KW results for: Healthy vs TR vs PR at Timepoint 2 
Found 1005 features from Kruskal-Wallis test. Generating outputs for all...
Finished generating individual Kruskal-Wallis feature outputs.

==================================================
STARTING WILCOXON ANALYSIS FOR MODE: POS 
==================================================
--------------------------------------------------
Running Wilcoxon Analysis: Paired TR Timepoint 1 vs 2 
Successfully saved Wilcoxon results for: Paired TR Timepoint 1 vs 2 
Found 1005 features from Wilcoxon test. Generating outputs for all...
Finished generating individual Wilcoxon feature outputs.
--------------------------------------------------
Running Wilcoxon Analysis: Paired PR Timepoint 1 vs 2 
Successfully saved Wilcoxon results for: Paired PR Timepoint 1 vs 2 
Found 1005 features from Wilcoxon test. Generating outputs for all...
Finished generating individual Wilcoxon feature outputs.
--------------------------------------------------
Running Wilcoxon Analysis: Unpaired Healthy vs TR at Timepoint 1 
Successfully saved Wilcoxon results for: Unpaired Healthy vs TR at Timepoint 1 
Found 1005 features from Wilcoxon test. Generating outputs for all...
Finished generating individual Wilcoxon feature outputs.

==================================================
STARTING INTERACTION ANALYSIS FOR MODE: POS 
==================================================
--------------------------------------------------
Running Interaction Analysis: Interaction of Change (T2-T1) between TR and PR groups 
Successfully saved Interaction results for: Interaction of Change (T2-T1) between TR and PR groups 

All analyses completed successfully!
R version 4.4.2 (2024-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 26100)

Matrix products: default


locale:
[1] LC_COLLATE=English_Israel.utf8  LC_CTYPE=English_Israel.utf8    LC_MONETARY=English_Israel.utf8 LC_NUMERIC=C                    LC_TIME=English_Israel.utf8    

time zone: Asia/Jerusalem
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] rstatix_0.7.2        ggpubr_0.6.1         scales_1.4.0         vegan_2.7-1          permute_0.9-8        ggrepel_0.9.6        pals_1.10            randomForest_4.7-1.2
 [9] ropls_1.38.0         lubridate_1.9.4      forcats_1.0.0        stringr_1.5.1        dplyr_1.1.4          purrr_1.1.0          readr_2.1.5          tidyr_1.3.1         
[17] tibble_3.3.0         ggplot2_3.5.2        tidyverse_2.0.0     

loaded via a namespace (and not attached):
 [1] tidyselect_1.2.1            farver_2.1.2                timechange_0.3.0            lifecycle_1.0.4             cluster_2.1.6               statmod_1.5.0              
 [7] magrittr_2.0.3              compiler_4.4.2              rlang_1.1.5                 tools_4.4.2                 calibrate_1.7.7             ggsignif_0.6.4             
[13] S4Arrays_1.6.0              labeling_0.4.3              DelayedArray_0.32.0         mapproj_1.2.12              RColorBrewer_1.1-3          abind_1.4-8                
[19] withr_3.0.2                 BiocGenerics_0.52.0         grid_4.4.2                  stats4_4.4.2                colorspace_2.1-1            MASS_7.3-61                
[25] MultiAssayExperiment_1.32.0 dichromat_2.0-0.1           SummarizedExperiment_1.36.0 cli_3.6.4                   crayon_1.5.3                ragg_1.4.0                 
[31] generics_0.1.4              rstudioapi_0.17.1           httr_1.4.7                  tzdb_0.5.0                  zlibbioc_1.52.0             splines_4.4.2              
[37] maps_3.4.3                  parallel_4.4.2              BiocManager_1.30.26         XVector_0.46.0              matrixStats_1.5.0           vctrs_0.6.5                
[43] Matrix_1.7-1                carData_3.0-5               jsonlite_2.0.0              car_3.1-3                   IRanges_2.40.1              hms_1.1.3                  
[49] S4Vectors_0.44.0            Formula_1.2-5               qqman_0.1.9                 systemfonts_1.2.3           limma_3.62.2                MultiDataSet_1.34.0        
[55] glue_1.8.0                  stringi_1.8.7               gtable_0.3.6                GenomeInfoDb_1.42.3         GenomicRanges_1.58.0        UCSC.utils_1.2.0           
[61] pillar_1.11.0               GenomeInfoDbData_1.2.13     R6_2.6.1                    textshaping_1.0.1           lattice_0.22-6              Biobase_2.66.0             
[67] backports_1.5.0             broom_1.0.9                 Rcpp_1.1.0                  SparseArray_1.6.2           nlme_3.1-166                mgcv_1.9-1                 
[73] MatrixGenerics_1.18.1       pkgconfig_2.0.3            
pproj_1.2.12              RColorBrewer_1.1-3          abind_1.4-8                
[19] withr_3.0.2                 BiocGenerics_0.52.0         grid_4.4.2                  stats4_4.4.2                colorspace_2.1-1            MASS_7.3-61                
[25] MultiAssayExperiment_1.32.0 dichromat_2.0-0.1           SummarizedExperiment_1.36.0 cli_3.6.4                   crayon_1.5.3                ragg_1.4.0                 
[31] generics_0.1.4              rstudioapi_0.17.1           httr_1.4.7                  tzdb_0.5.0                  zlibbioc_1.52.0             splines_4.4.2              
[37] maps_3.4.3                  parallel_4.4.2              BiocManager_1.30.26         XVector_0.46.0              matrixStats_1.5.0           vctrs_0.6.5                
[43] Matrix_1.7-1                carData_3.0-5               jsonlite_2.0.0              car_3.1-3                   IRanges_2.40.1              hms_1.1.3                  
[49] S4Vectors_0.44.0            Formula_1.2-5               qqman_0.1.9                 systemfonts_1.2.3           limma_3.62.2                MultiDataSet_1.34.0        
[55] glue_1.8.0                  stringi_1.8.7               gtable_0.3.6                GenomeInfoDb_1.42.3         GenomicRanges_1.58.0        UCSC.utils_1.2.0           
[61] pillar_1.11.0               GenomeInfoDbData_1.2.13     R6_2.6.1                    textshaping_1.0.1           lattice_0.22-6              Biobase_2.66.0             
[67] backports_1.5.0             broom_1.0.9                 Rcpp_1.1.0                  SparseArray_1.6.2           nlme_3.1-166                mgcv_1.9-1                 
[73] MatrixGenerics_1.18.1       pkgconfig_2.0.3            
     crayon_1.5.3                ragg_1.4.0                 
[31] generics_0.1.4              rstudioapi_0.17.1           httr_1.4.7                  tzdb_0.5.0                  zlibbioc_1.52.0             splines_4.4.2              
[37] maps_3.4.3                  parallel_4.4.2              BiocManager_1.30.26         XVector_0.46.0              matrixStats_1.5.0           vctrs_0.6.5                
[43] Matrix_1.7-1                carData_3.0-5               jsonlite_2.0.0              car_3.1-3                   IRanges_2.40.1              hms_1.1.3                  
[49] S4Vectors_0.44.0            Formula_1.2-5               qqman_0.1.9                 systemfonts_1.2.3           limma_3.62.2                MultiDataSet_1.34.0        
[55] glue_1.8.0                  stringi_1.8.7               gtable_0.3.6                GenomeInfoDb_1.42.3         GenomicRanges_1.58.0        UCSC.utils_1.2.0           
[61] pillar_1.11.0               GenomeInfoDbData_1.2.13     R6_2.6.1                    textshaping_1.0.1           lattice_0.22-6              Biobase_2.66.0             
[67] backports_1.5.0             broom_1.0.9                 Rcpp_1.1.0                  SparseArray_1.6.2           nlme_3.1-166                mgcv_1.9-1                 
[73] MatrixGenerics_1.18.1       pkgconfig_2.0.3            
ts for all...
Finished generating individual Wilcoxon feature outputs.
--------------------------------------------------
Running Wilcoxon Analysis: Unpaired Healthy vs TR at Timepoint 1 
Successfully saved Wilcoxon results for: Unpaired Healthy vs TR at Timepoint 1 
Found 1005 features from Wilcoxon test. Generating outputs for all...
Finished generating individual Wilcoxon feature outputs.

==================================================
STARTING INTERACTION ANALYSIS FOR MODE: POS 
==================================================
--------------------------------------------------
Running Interaction Analysis: Interaction of Change (T2-T1) between TR and PR groups 
Successfully saved Interaction results for: Interaction of Change (T2-T1) between TR and PR groups 

All analyses completed successfully!
R version 4.4.2 (2024-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 26100)

Matrix products: default


locale:
[1] LC_COLLATE=English_Israel.utf8  LC_CTYPE=English_Israel.utf8    LC_MONETARY=English_Israel.utf8 LC_NUMERIC=C                    LC_TIME=English_Israel.utf8    

time zone: Asia/Jerusalem
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] rstatix_0.7.2        ggpubr_0.6.1         scales_1.4.0         vegan_2.7-1          permute_0.9-8        ggrepel_0.9.6        pals_1.10            randomForest_4.7-1.2
 [9] ropls_1.38.0         lubridate_1.9.4      forcats_1.0.0        stringr_1.5.1        dplyr_1.1.4          purrr_1.1.0          readr_2.1.5          tidyr_1.3.1         
[17] tibble_3.3.0         ggplot2_3.5.2        tidyverse_2.0.0     

loaded via a namespace (and not attached):
 [1] tidyselect_1.2.1            farver_2.1.2                timechange_0.3.0            lifecycle_1.0.4             cluster_2.1.6               statmod_1.5.0              
 [7] magrittr_2.0.3              compiler_4.4.2              rlang_1.1.5                 tools_4.4.2                 calibrate_1.7.7             ggsignif_0.6.4             
[13] S4Arrays_1.6.0              labeling_0.4.3              DelayedArray_0.32.0         mapproj_1.2.12              RColorBrewer_1.1-3          abind_1.4-8                
[19] withr_3.0.2                 BiocGenerics_0.52.0         grid_4.4.2                  stats4_4.4.2                colorspace_2.1-1            MASS_7.3-61                
[25] MultiAssayExperiment_1.32.0 dichromat_2.0-0.1           SummarizedExperiment_1.36.0 cli_3.6.4                   crayon_1.5.3                ragg_1.4.0                 
[31] generics_0.1.4              rstudioapi_0.17.1           httr_1.4.7                  tzdb_0.5.0                  zlibbioc_1.52.0             splines_4.4.2              
[37] maps_3.4.3                  parallel_4.4.2              BiocManager_1.30.26         XVector_0.46.0              matrixStats_1.5.0           vctrs_0.6.5                
[43] Matrix_1.7-1                carData_3.0-5               jsonlite_2.0.0              car_3.1-3                   IRanges_2.40.1              hms_1.1.3                  
[49] S4Vectors_0.44.0            Formula_1.2-5               qqman_0.1.9                 systemfonts_1.2.3           limma_3.62.2                MultiDataSet_1.34.0        
[55] glue_1.8.0                  stringi_1.8.7               gtable_0.3.6                GenomeInfoDb_1.42.3         GenomicRanges_1.58.0        UCSC.utils_1.2.0           
[61] pillar_1.11.0               GenomeInfoDbData_1.2.13     R6_2.6.1                    textshaping_1.0.1           lattice_0.22-6              Biobase_2.66.0             
[67] backports_1.5.0             broom_1.0.9                 Rcpp_1.1.0                  SparseArray_1.6.2           nlme_3.1-166                mgcv_1.9-1                 
[73] MatrixGenerics_1.18.1       pkgconfig_2.0.3            
