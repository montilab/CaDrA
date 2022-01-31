
#' Top 'N' evaluate
#' 
#' Generates and evaluates stepwise search results for the top 'N' starting indices, checking for overlapping resulting features from each case. This function is mainly used to evaluate search results over the top 'N' best starting features for a given dataset.
#' @param ESet an ordered expression set object with the same sample ordering and features as processed by the stepwise.search() function when performing a step-wise heuristic search
#' @param input_score a vector of ranked or continuous values (required). 
#' @param method a character string specifying the method used to compute scores for features, must be one of "ks" or "wilcox" or "mi" (mutually exclusive method from REVEALER) or "custom" (a personal customization method). If input_score contains ranked scores, then 'ks' method is used by default. Otherwise, 'mi" is the default method
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided", "greater" or "less". Default is "less" for left-skewed significance testing.
#' @param metric a character string specifying which metric to use for candidate search. One of either 'pval' or 'stat' may be used, corresponding to the score p-value or statistic. Default is 'pval'
#' @param search_method a character string specifying which method to perform or filter out the best candidates. Default is 'forward'.
#' @param max_size an integer specifying the maximum size a meta-feature can extend do for a given search. Default is 7.
#' @param best_score_only a logical indicating whether or not to only return the best meta-feature score over the top 'N' evaluation. Default is TRUE
#' @param N an integer specifying the number of features to start the search over, starting from the top 'N' features in each case. Default is 1
#' @param do_plot a logical indicating whether you want to plot the resulting evaluation matrix. Default is TRUE
#' @param verbose a logical indicating whether or not to verbose diagnostic messages. Default is TRUE. 
#' 
#' @return Default is a list of lists, where each list entry is one that is returned by the stepwise search run for a given starting index (See stepwise.search()). If best_score_only is set to TRUE, only the best score over the top N space is returned (useful for permutation-based testing)
#' 1's and 0's represent whether a feature in any given row is present in a meta-feature along with a starting feature in the corresponding column.
#' @export
#' @import Biobase gplots
topn_eval <- function(
  ESet,
  input_score,
  method = "ks", 
  alternative = "less", 
  metric = "pval", 
  search_method = "both", 
  max_size = 7,
  best_score_only = TRUE,
  N = 1,
  do_plot = TRUE,
  verbose = TRUE
){
  
  # Set up verbose option
  options(verbose=verbose)
  
  if (N > nrow(ESet))
    stop("Please specify an N value that is less than the number of features in the ESet..\n")
  
  if (N > 10)
    warning("N value specified is greater than 10. This may result in longer search time..\n")
  
  verbose("Evaluating search over top features: ", 1:N, "\n\n")
  
  #Performs stepwise search over top N indices
  topn_l <- sapply(1:N, function(x){ 
    candidate_search(ES=ESet, input_score=input_score, method=method, alternative=alternative, metric=metric, search_method=search_method, max_size=max_size, search_start=x) 
  }, simplify = FALSE) 
  
  if(best_score_only==TRUE){
    scores_l <- lapply(topn_l, "[[", 2)
    
    # Working with scores for each top N run
    s <- unlist(scores_l)
    
    #Fetch the best score from the iterations
    # This ASSUMES you're using metric = "pval"
    # NEEDS UPDATING TO ACCOMODATE STATISTIC 
    best_score <- s[order(s)][1] #Based on the p-values, the lowest value will be the most sig 
    
    return(best_score)
  } #best_score_only
  
  
  if(do_plot){
    topn_plot(topN_list=topn_l)  
  } #do_plot
  
  return(topn_l) #Default is to return the top N stepwise search results as a list of lists
  
}

