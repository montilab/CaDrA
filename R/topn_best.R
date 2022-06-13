
#' Top 'N' best meta-features
#' 
#' Takes the resulting list of meta-features returned from topn_eval() function and fetches the meta-features with best (local) score
#' @param topn_list A nested list object that is returned by topn_eval() function when best_score_only is set to FALSE. This contains both the ESets as well as the scores and its input scores for each resulting meta-feature in the top 'N' search mode.
#' @return A list containing the (local) best meta-feature ESet, as well as its corresponding search score and input scores
#' @examples
#' 
#' # Load pre-computed Top-N list generated for sim.ES dataset
#' data(topn.list)
#' 
#' # Get the best meta-features list
#' topn_best_meta <- topn_best(topn_list=topn.list) 
#' 
#' @export
#' @import Biobase
topn_best <- function(topn_list){
  
  # get best score list
  scores_l <- lapply(1:length(topn_list), function(l){ topn_list[[l]][['Score']] })
  
  # Working with scores for each top N run
  scores <- unlist(scores_l)
  
  #
  # NEEDS UPDATING TO ACCOMODATE STATISTIC 
  #
  
  # Fetch the index housing the best ESet (this wil be the one with the best score)
  n <- which.min(scores)
  
  # Also store the score
  top_score <- min(scores)
  
  # Corresponding ESet object
  best_meta <- topn_list[[n]]$ESet
  
  # Correspoding input_score
  best_input_score <- topn_list[[n]]$input_score
  
  return(list("ESet" = best_meta, "Score" = top_score, "input_score" = best_input_score))
  
}
