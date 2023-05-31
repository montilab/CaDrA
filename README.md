
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CaDrA

![build](https://github.com/montilab/cadra/workflows/rcmdcheck/badge.svg)
![Gitter](https://img.shields.io/gitter/room/montilab/cadra) ![GitHub
issues](https://img.shields.io/github/issues/montilab/cadra) ![GitHub
last commit](https://img.shields.io/github/last-commit/montilab/cadra)

**Ca**ndidate **Dr**ivers **A**nalysis: Multi-Omic Search for Candidate
Drivers of Functional Signatures

**CaDrA** is an R package that supports a heuristic search framework
aimed at identifying candidate drivers of a molecular phenotype of
interest.

The main function takes two inputs:

1)  A binary multi-omics dataset, which can be represented as a matrix
    of binary features or a **SummarizedExperiment** class object where
    the rows are 1/0 vectors indicating the presence/absence of ‘omics’
    features (e.g. somatic mutations, copy number alterations,
    epigenetic marks, etc.), and the columns are the samples.
2)  A molecular phenotype of interest, which can be represented as a
    vector of continuous scores (e.g. protein expression, pathway
    activity, etc.)

Based on these two inputs, **CaDrA** implements a forward/backward
search algorithm to find a set of features that together is maximally
associated with the observed input scores, based on one of several
scoring functions (*Kolmogorov-Smirnov*, *Conditional Mutual
Information*, *Wilcoxon*, *custom-defined scoring function*), making it
useful to find complementary omics features likely driving the input
molecular phenotype.

For more information, please see the associated manuscript [Kartha et
al. (2019)](https://www.frontiersin.org/articles/10.3389/fgene.2019.00121/full)

## (1) Installation

``` r
library(devtools)
devtools::install_github("montilab/CaDrA")
```

## (2) Quickstart

``` r
library(CaDrA)
library(SummarizedExperiment)
```

## (3) CaDrA Query of BRCA YAP/TAZ Activity

### (i) Load & Format Data Inputs

``` r

## Read in BRCA GISTIC+Mutation object
data(BRCA_GISTIC_MUT_SIG)
eset_mut_scna <- BRCA_GISTIC_MUT_SIG

## Read in input score
data(TAZYAP_BRCA_ACTIVITY)
input_score <- TAZYAP_BRCA_ACTIVITY

## Samples to keep based on the overlap between the two inputs
overlap <- intersect(names(input_score), colnames(eset_mut_scna))
eset_mut_scna <- eset_mut_scna[,overlap]
input_score <- input_score[overlap]

## Binarize FS to only have 0's and 1's
assay(eset_mut_scna)[assay(eset_mut_scna) > 1] <- 1.0

## Pre-filter FS based on occurrence frequency
eset_mut_scna_flt <- CaDrA::prefilter_data(
  FS = eset_mut_scna,
  max_cutoff = 0.6,  # max event frequency (60%)
  min_cutoff = 0.03  # min event frequency (3%)
)  
```

### (ii) Run CaDrA

Here, we repeat the candidate search starting from each of the top ‘N’
features and report the combined results as a heatmap (to summarize the
number of times each feature is selected across repeated runs).

IMPORTANT NOTE: The legacy function `topn_eval()` is equivalent to the
recommended `candidate_search()` function

``` r

topn_res <- CaDrA::candidate_search(
  FS = eset_mut_scna_flt,
  input_score = input_score,
  method = "ks_pval",          # Use Kolmogorow-Smirnow scoring function 
  weight = NULL,               # If weights is provided, perform a weighted-KS test
  alternative = "less",        # Use one-sided hypothesis testing
  search_method = "both",      # Apply both forward and backward search
  top_N = 7,                   # Evaluate top 7 starting points for each search
  max_size = 7,                # Maximum size a meta-feature matrix can extend to
  do_plot = FALSE,             # Plot after finding the best features
  best_score_only = FALSE      # Return meta-feature set, observed input scores and calculated best score
)
```

### (iii) Visualize Best Results

``` r

## Fetch the meta-feature set corresponding to its best scores over top N features searches
topn_best_meta <- CaDrA::topn_best(topn_res)

# Visualize the best results with the meta-feature plot
CaDrA::meta_plot(topn_best_list = topn_best_meta, input_score_label = "YAP/TAZ Activity")
```

![](./man/figures/visualize.best-1.png)<!-- -->

### (iv) Summarize Top N Results

``` r

# Evaluate results across top N features you started from
CaDrA::topn_plot(topn_res) 
```

![](./man/figures/summarize-1.png)<!-- -->

## (4) SessionInfo

``` r
sessionInfo()
R version 4.3.0 (2023-04-21)
Platform: x86_64-apple-darwin20 (64-bit)
Running under: macOS Ventura 13.3.1

Matrix products: default
BLAS:   /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRblas.0.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/New_York
tzcode source: internal

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] CaDrA_0.99.5                SummarizedExperiment_1.30.1
 [3] Biobase_2.60.0              GenomicRanges_1.52.0       
 [5] GenomeInfoDb_1.36.0         IRanges_2.34.0             
 [7] S4Vectors_0.38.1            BiocGenerics_0.46.0        
 [9] MatrixGenerics_1.12.0       matrixStats_0.63.0         

loaded via a namespace (and not attached):
 [1] gtable_0.3.3            xfun_0.39               ggplot2_3.4.2          
 [4] caTools_1.18.2          lattice_0.21-8          vctrs_0.6.2            
 [7] tools_4.3.0             bitops_1.0-7            generics_0.1.3         
[10] parallel_4.3.0          tibble_3.2.1            fansi_1.0.4            
[13] highr_0.10              pkgconfig_2.0.3         R.oo_1.25.0            
[16] Matrix_1.5-4            KernSmooth_2.23-21      lifecycle_1.0.3        
[19] GenomeInfoDbData_1.2.10 R.cache_0.16.0          farver_2.1.1           
[22] compiler_4.3.0          stringr_1.5.0           gplots_3.1.3           
[25] munsell_0.5.0           codetools_0.2-19        misc3d_0.9-1           
[28] htmltools_0.5.5         RCurl_1.98-1.12         yaml_2.3.7             
[31] pillar_1.9.0            crayon_1.5.2            MASS_7.3-60            
[34] R.utils_2.12.2          DelayedArray_0.26.3     iterators_1.0.14       
[37] foreach_1.5.2           gtools_3.9.4            tidyselect_1.2.0       
[40] digest_0.6.31           stringi_1.7.12          dplyr_1.1.2            
[43] reshape2_1.4.4          labeling_0.4.2          fastmap_1.1.1          
[46] grid_4.3.0              colorspace_2.1-0        cli_3.6.1              
[49] magrittr_2.0.3          S4Arrays_1.0.4          utf8_1.2.3             
[52] withr_2.5.0             scales_1.2.1            rmarkdown_2.21         
[55] XVector_0.40.0          R.methodsS3_1.8.2       evaluate_0.21          
[58] knitr_1.42              tcltk_4.3.0             doParallel_1.0.17      
[61] rlang_1.1.1             Rcpp_1.0.10             glue_1.6.2             
[64] ppcor_1.1               rstudioapi_0.14         R6_2.5.1               
[67] plyr_1.8.8              zlibbioc_1.46.0        
```

## (5) Acknowledgements

This project is funded in part by the
[NIH/NIDCR](https://www.nidcr.nih.gov/) (3R01DE030350-01A1S1,
R01DE031831), [Find the Cause Breast Cancer
Foundation](https://findthecausebcf.org), and
[NIH/NIA](https://www.nia.nih.gov/) (UH3 AG064704).
