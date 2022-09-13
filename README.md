### Impact of elexacaftor-tezacaftor-ivacaftor therapy on the cystic fibrosis airway microbial metagenome 
<br>
Pallenberg, Sophia T.<sup>1,2</sup>, Pust, Marie-Madlen<sup>1,2</sup>, Rosenboom, Ilona<sup>1,2</sup>, Hansen, Gesine<sup>1,2</sup>, Wiehlmann, Lutz<sup>1,2</sup> Dittrich, Anna-Maria<sup>1,2</sup>, Tümmler, Burkhard<sup>1,2</sup> 
<br><br>
<sup>1</sup>Department for Pediatric Pneumology, Allergology and Neonatology, Hannover Medical School, Hannover, Germany <br>
<sup>2</sup>German Center for Lung Research, Biomedical Research in Endstage and Obstructive Lung Disease (BREATH), Hannover Medical School, Hannover, Germany <br>

<br><br>

# R session 
```
R version 4.1.1 (2021-08-10)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19042)

Matrix products: default

locale:
[1] LC_COLLATE=German_Germany.1252  LC_CTYPE=German_Germany.1252    LC_MONETARY=German_Germany.1252 LC_NUMERIC=C                    LC_TIME=German_Germany.1252    

attached base packages:
[1] tcltk     stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] forcats_0.5.1              tibble_3.1.3               tidyverse_1.3.1            ggtext_0.1.1               GGally_2.1.2               rcompanion_2.4.1          
 [7] PerformanceAnalytics_2.0.4 xts_0.12.1                 zoo_1.8-9                  psych_2.1.6                magrittr_2.0.1             ggsignif_0.6.2            
[13] taxa_0.4.0                 vctrs_0.3.8                metacoder_0.3.5            RColorBrewer_1.1-2         repr_1.1.4                 seriation_1.3.0           
[19] corrplot_0.90              asbio_1.7                  stringr_1.4.0              Boruta_7.0.0               randomForest_4.6-14        readxl_1.3.1              
[25] plyr_1.8.6                 hrbrthemes_0.8.0           ggthemes_4.2.4             ggpubr_0.4.0               matrixStats_0.60.0         igraph_1.2.6              
[31] dplyr_1.0.7                Hmisc_4.5-0                Formula_1.2-4              survival_3.2-11            tidyr_1.1.3                viridis_0.6.1             
[37] viridisLite_0.4.0          ggrepel_0.9.1              ggplot2_3.3.5              vegan_2.5-7                lattice_0.20-45            permute_0.9-5             
[43] purrr_0.3.4                readr_2.0.1               

loaded via a namespace (and not attached):
  [1] utf8_1.2.2           tidyselect_1.1.1     htmlwidgets_1.5.3    grid_4.1.1           TSP_1.1-10           combinat_0.0-8       ranger_0.13.1       
  [8] munsell_0.5.0        codetools_0.2-18     withr_2.4.2          colorspace_2.0-2     knitr_1.33           rstudioapi_0.13      stats4_4.1.1        
 [15] DescTools_0.99.42    Rttf2pt1_1.3.9       labeling_0.4.2       gWidgets2tcltk_1.0-8 mnormt_2.0.2         bit64_4.0.5          farver_2.1.0        
 [22] generics_0.1.0       TH.data_1.0-10       xfun_0.25            R6_2.5.1             cachem_1.0.6         reshape_0.8.8        assertthat_0.2.1    
 [29] scales_1.1.1         vroom_1.5.4          multcomp_1.4-17      nnet_7.3-16          rootSolve_1.8.2.2    gtable_0.3.0         multcompView_0.1-8  
 [36] lmom_2.8             sandwich_3.0-1       rlang_0.4.11         systemfonts_1.0.2    scatterplot3d_0.3-41 splines_4.1.1        rstatix_0.7.0       
 [43] extrafontdb_1.0      broom_0.7.9          checkmate_2.0.0      modelr_0.1.8         yaml_2.2.1           reshape2_1.4.4       abind_1.4-5         
 [50] backports_1.2.1      gridtext_0.1.4       extrafont_0.17       tools_4.1.1          ellipsis_0.3.2       proxy_0.4-26         Rcpp_1.0.7          
 [57] progress_1.2.2       base64enc_0.1-3      prettyunits_1.1.1    rpart_4.1-15         cowplot_1.1.1        deSolve_1.30         haven_2.4.3         
 [64] cluster_2.1.2        fs_1.5.0             data.table_1.14.0    pixmap_0.4-12        openxlsx_4.2.4       reprex_2.0.1         lmtest_0.9-38       
 [71] tmvnsim_1.0-2        mvtnorm_1.1-2        hms_1.1.0            evaluate_0.14        rio_0.5.27           jpeg_0.1-9           gridExtra_2.3       
 [78] compiler_4.1.1       crayon_1.4.1         htmltools_0.5.1.1    mgcv_1.8-36          tzdb_0.1.2           libcoin_1.0-8        expm_0.999-6        
 [85] Exact_2.1            lubridate_1.8.0      DBI_1.1.1            gWidgets2_1.0-9      dbplyr_2.1.1         MASS_7.3-54          boot_1.3-28         
 [92] Matrix_1.3-4         car_3.0-11           cli_3.0.1            quadprog_1.5-8       parallel_4.1.1       pkgconfig_2.0.3      registry_0.5-1      
 [99] coin_1.4-1           foreign_0.8-81       xml2_1.3.2           foreach_1.5.1        rvest_1.0.2          digest_0.6.27        rmarkdown_2.10      
[106] cellranger_1.1.0     htmlTable_2.2.1      nortest_1.0-4        gld_2.6.2            gdtools_0.2.4        curl_4.3.2           modeltools_0.2-23   
[113] lifecycle_1.0.0      nlme_3.1-152         jsonlite_1.7.2       carData_3.0-4        fansi_0.5.0          pillar_1.6.2         ggsci_2.9           
[120] fastmap_1.1.0        httr_1.4.2           plotrix_3.8-1        glue_1.4.2           zip_2.2.0            png_0.1-7            iterators_1.0.13    
[127] bit_4.0.4            class_7.3-19         stringi_1.7.3        latticeExtra_0.6-29  memoise_2.0.0        e1071_1.7-8     
```
