
#' Top 'N' best
#' 
#' Takes the resulting list of meta-features returned from topn_eval() function and fetches the meta-feature with best (local) score
#' @param topn_list The nested list object that is returned by topn_eval() function when best.score.only is set to FALSE. This contains both the ESets as well as the scores for each resulting meta-feature in the top 'N'  search mode.
#' @return A list containing the (local) best meta-feature ESet, as well as its corresponding search score
#' @export
#' @import Biobase
topn_best <- function(topn_list){
  
  # Fetch the index housing the best ESet (this wil be the one with the best score)
  n <- which.min(base::sapply(topn_list,"[[",2))
  
  # Also store the score
  top_score <- min(base::sapply(topn_list,"[[",2))
  
  # Corresponding ESet object
  best_meta <- topn_list[[n]]$ESet
  
  return(list("ESet" = best_meta, "Score" = top_score))
  
}
