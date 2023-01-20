
#' Top 'N' Best Meta-Features
#'
#' Takes the resulting list of meta-features returned from \code{candidate_search()}
#' and fetches the meta-features with the best score
#' @param topn_list A nested list of objects that is returned from \code{candidate_search()}
#' by setting \code{best_score_only = FALSE}. The nested list contains the
#' meta-feature sets as well as the observed input scores and its corresponding best score 
#' for each top 'N' search.
#' @return A list of objects containing the best meta-feature, its corresponding best score 
#' and associated input scores
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
#' @import SummarizedExperiment
topn_best <- function(topn_list){

  # get best score list
  scores_l <- lapply(seq_along(topn_list),
                     function(l){ topn_list[[l]][['score']] })

  # Working with scores for each top N run
  scores <- unlist(scores_l)

  # Obtain the metric used in candidate search
  metric <- lapply(
    seq_along(topn_list),
    function(l){ topn_list[[l]][['metric']] }
  ) |>
    unlist() |>
    unique()

  # Fetch the best score from the iterations
  # NEEDS UPDATING TO ACCOMODATE STATISTIC
  if(metric == "pval"){
    # Fetch the index housing the best FS
    # Based on the p-values, the lowest value will be the most significant
    n <- which.min(scores)

    # Also obtain the best score
    top_score <- min(scores)
  }else{
    # Fetch the index housing the best FS
    # Based on the statistics, the largest value will be the most significant
    n <- which.max(scores)

    # Obtain the best score
    top_score <- max(scores)
  }

  # Corresponding FS object
  best_meta <- topn_list[[n]]$feature_set

  # Correspoding input_score
  best_input_score <- topn_list[[n]]$input_score

  return(
    list(
      "feature_set" = best_meta,
      "input_score" = best_input_score,
      "score" = top_score
    )
  )

}
