#' Top 'N' Best Meta-Features
#'
#' Take the resulting list of meta-features returned from 
#' \code{candidate_search} over top N feature searches and 
#' fetch the meta-feature with the best score.
#' @param topn_list A nested list of objects that are returned 
#' from \code{candidate_search} using the following parameters:
#' \code{FS = sim_FS}, \code{input_score = sim_Scores}, 
#' \code{top_N = 7}, \code{method = "ks_pval"}, \code{alternative = "less"}, 
#' \code{search_method = "both"}, \code{max_size = 10}, and 
#' \code{best_score_only = FALSE}. 
#' 
#' @return A list of objects containing the best meta-feature matrix, its 
#' corresponding best score, its observed input scores, rank of best 
#' meta-features based on their scores, its marginal and cumulative best scores.
#' 
#' @examples
#'
#' # Load pre-computed Top-N list generated for sim_FS dataset
#' data(topn_list)
#'
#' # Get the best meta-features list
#' topn_best_meta <- topn_best(topn_list = topn_list)
#'
#' @export
topn_best <- function(topn_list){

  # get best score list
  scores_l <- lapply(seq_along(topn_list),
                     function(l){ topn_list[[l]][['score']] })
  
  # Working with scores for each top N run
  scores <- unlist(scores_l)
  
  # Fetch the index housing the best FS
  # Based on the statistics, the largest value will be the most significant
  n <- which.max(scores)
  
  # Obtain the best score
  top_score <- scores[n]
  
  # Corresponding FS object
  best_meta <- topn_list[[n]]$feature_set

  # Corresponding input_score
  best_input_score <- topn_list[[n]]$input_score
  
  # get indices of best features
  best_meta_indices <- topn_list[[n]]$best_indices
  
  # Get scores of best features
  marginal_best_scores <- topn_list[[n]]$marginal_best_scores

  # Get scores of best features
  cumulative_best_scores <- topn_list[[n]]$cumulative_best_scores  
  
  return(
    list(
      "feature_set" = best_meta,
      "input_score" = best_input_score,
      "score" = top_score,
      "best_indices" = best_meta_indices,
      "marginal_best_scores" = marginal_best_scores,
      "cumulative_best_scores" = cumulative_best_scores
    )
  )

}
