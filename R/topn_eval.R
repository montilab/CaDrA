
#' Evaluate Top 'N' Best Features
#' 
#' Generates and evaluates candidate search results for the top 'N' starting indices, checking for overlapping resulting features from each case. This function is mainly used to evaluate search results over the top 'N' best starting features for a given dataset.
#' @param ES an expression set of binary features (required). It must be a BioBase expressionSet object. The rownames of the expression set must contain unique features which are used in the search.   
#' @param input_score a vector of continuous values for a target profile (required). The input_score must have names or labels that matches the colnames of the expression matrix.
#' @param method a character string specifies a method to compute the score for each feature (\code{"ks"} or \code{"wilcox"} or \code{"revealer"} (conditional mutual information from REVEALER) or \code{"custom"} (a customized method)). Default is \code{ks}.
#' @param custom_function if method is \code{"custom"}, specifies the customized function here. Default is \code{NULL}.
#' @param custom_parameters if method is \code{"custom"}, specifies a list of arguments to be passed to the custom_function(). Default is \code{NULL}.
#' @param alternative a character string specifies an alternative hypothesis testing (\code{"two.sided"} or \code{"greater"} or \code{"less"}). Default is \code{less} for left-skewed significance testing.
#' @param metric a character string specifies a metric to use for candidate search criteria. \code{"pval"} or \code{"stat"} may be used, corresponding to the score p-value or statistic. Default is \code{pval}.
#' @param weights a vector of weights use to perform a weighted-KS testing. Default is \code{NULL}.   
#' @param target_match a direction of target matching (\code{"negative"} or \code{"positive"}) from REVEALER. Use \code{"positive"} to match the higher values of the target, \code{"negative"} to match the lower values. Default is \code{positive}. 
#' @param top_N an integer specifies the number of features to start the search over, starting from the top 'N' features in each case. Default is \code{1}.
#' @param search_method a character string specifies a method to filter out the best candidates (\code{"forward"} or \code{"both"}). Default is \code{both} (backward and forward).
#' @param max_size an integer specifies a maximum size that a meta-feature can extend to do for a given search. Default is \code{7}.
#' @param best_score_only a logical value indicates whether or not the function should return only the score corresponding to the search results. Default is \code{FALSE}.
#' @param do_plot a logical value indicates whether or not to plot the resulting evaluation matrix. Default is \code{TRUE}.
#' @param verbose a logical value indicates whether or not to print the diagnostic messages. Default is \code{FALSE}. 
#' 
#' @return By default, this function will return a list of lists, where each list entry is one that is returned by the candidate search for a given starting index (See \code{candidate_search()}). If \code{best_score_only} is set to \code{TRUE}, only the best score over the top N space is returned (useful for permutation-based testing)
#' 1's and 0's represent whether a feature in any given row is present in a meta-feature along with a starting feature in the corresponding column.
#' @examples
#' 
#' # Load R library
#' library(Biobase)
#'
#' # Load pre-computed expression set
#' data(sim.ES)
#' 
#' # set seed
#' set.seed(123)
#' 
#' # Provide a vector of continuous scores for a target profile with names to each score value 
#' input_score = rnorm(n = ncol(sim.ES))
#' names(input_score) <- colnames(sim.ES)
#' 
#' # Define additional parameters and run the function
#' topn_eval <- topn_eval(
#'   ES = sim.ES, input_score = input_score, method = "ks",
#'   alternative = "less", metric = "pval", top_N = 3, search_method = "both", 
#'   max_size = 7, best_score_only = FALSE
#' )
#' 
#' @export
#' @import Biobase gplots
topn_eval <- function(
  ES,
  input_score,
  method = "ks", 
  custom_function = NULL,
  custom_parameters = NULL,
  alternative = "less", 
  metric = "pval", 
  weights = NULL,
  target_match = "positive",
  top_N = 1,
  search_method = "both", 
  max_size = 7,
  best_score_only = TRUE,
  do_plot = TRUE,
  verbose = FALSE
){
  
  # Set up verbose option
  options(verbose = FALSE)
  
  if(top_N > nrow(ES))
    stop("Please specify an top_N value that is less than the number of features in the ES.\n")
  
  if(top_N > 10)
    warning("top_N value specified is greater than 10. This may result in a longer search time.\n")
  
  verbose("Evaluating search over top features: ", 1:top_N, "\n\n")
  
  # Performs candidate search over top top_N indices
  topn_l <- sapply(1:top_N, function(x){ 
    
    candidate_search(
      ES = ES, 
      input_score = input_score, 
      method = method, 
      custom_function = custom_function,
      custom_parameters = custom_parameters,
      alternative = alternative, 
      metric = metric, 
      weights = weights,
      target_match = target_match,
      search_start = x,
      search_method = search_method, 
      max_size = max_size,
      best_score_only = FALSE
    ) 
    
  }, simplify = FALSE) 
  
  # do_plot
  if(do_plot){
    
    topn_plot(topn_list = topn_l)  
    
  }
  
  # best_score_only
  if(best_score_only == TRUE){
    
    scores_l <- lapply(1:length(topn_l), function(l){ topn_l[[l]][['Score']] })
    
    # Working with scores for each top N run
    s <- unlist(scores_l)
    
    #Fetch the best score from the iterations
    # This ASSUMES you're using metric = "pval"
    # NEEDS UPDATING TO ACCOMODATE STATISTIC 
    best_score <- s[order(s)][1] #Based on the p-values, the lowest value will be the most significant 
    
    return(best_score)
    
  }
  
  return(topn_l) #Default is to return the top N stepwise search results as a list of lists
  
}

