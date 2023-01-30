
library(SummarizedExperiment)



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
eset_mut_scna_flt <- prefilter_data(
  FS = eset_mut_scna,
  max_cutoff = 0.6,  # max event frequency (60%)
  min_cutoff = 0.03  # min event frequency (3%)
)  



topn_res <- candidate_search(
  FS = eset_mut_scna_flt,
  input_score = input_score,
  method = "ks",               # Use Kolmogorow-Smirnow scoring function 
  weight = NULL,               # If weights is provided, perform a weighted-KS test
  alternative = "less",        # Use one-sided hypothesis testing
  metric = "pval",             # Use p-value to search for best features
  search_method = "both",      # Apply both forward and backward search
  top_N = 7,                   # Evaluate top 7 starting points for each search
  max_size = 7,                # Maximum size a meta-feature matrix can extend to
  do_plot = FALSE,             # Plot after finding the best features
  best_score_only = FALSE      # Return meta-feature set, observed input scores and calculated best score
)


## Fetch the meta-feature set corresponding to its best scores over top N features searches
topn_best_meta <- topn_best(topn_res)

# Visualize the best results with the meta-feature plot
meta_plot(topn_best_list = topn_best_meta, input_score_label = "YAP/TAZ Activity")
