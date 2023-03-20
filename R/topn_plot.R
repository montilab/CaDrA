#' Top 'N' Plot
#'
#' Generate a heatmap representation of overlapping meta-features across
#' top N feature searches using \code{candidate_search()} function
#' @param topn_list a list of objects returned from \code{candidate_search()} 
#' using simulated dataset \code{FS = sim_FS}, \code{input_score = sim_Scores}, 
#' \code{top_N = 7}, \code{method = "ks_pval"}, \code{alternative = "less"}, 
#' \code{search_method = "both"}, \code{max_size = 10},
#' and \code{best_score_only = FALSE} as inputs to the function.
#' 
#' The resulting list contains a set of meta-features in form of SummarizedExperiment 
#' object, a vector of observed input scores, and its corresponding best score 
#' over top N feature searches
#'
#' @return a heatmap of overlapping meta-features for a given top N feature searches
#' @examples
#'
#' # Load pre-computed Top-N list generated for sim_FS dataset
#' data(topn_list)
#'
#' # Get the overlapping top N plot
#' topn_plot(topn_list = topn_list)
#'
#' @export
#' @import gplots
#' @importFrom graphics legend
topn_plot <- function(
    topn_list
){

  # Get feature_set and best scores for top n features
  feature_set_l <- lapply(seq_along(topn_list),
                   function(l){ topn_list[[l]][['feature_set']] })
  
  scores_l <- lapply(seq_along(topn_list),
                     function(l){ topn_list[[l]][['score']] })

  # Get the list of feature names from each FS object
  f_list <- lapply(feature_set_l, rownames)

  # Get the union of all features that were returned across all top N runs
  f_union <- Reduce(f = union, f_list)

  f_checklist <- lapply(f_list, function(x, ref = f_union){
    return(f_union %in% x)
  })

  # Working with scores for each top N run
  scores <- unlist(scores_l)
  
  # Make a matrix indicating which features are found across each top n run
  mat <- do.call(cbind, f_checklist)*1
  
  # Assign row and column names to mat
  colnames(mat) <- names(scores)
  rownames(mat) <- f_union  
  
  if(ncol(mat) >= 2){
    
    # If ncol(mat) > 1, we can order matrix columns in decreasing order
    mat <- mat[, order(scores, decreasing = TRUE)]
    
    # Add the index number and its best score to each of the columns names of mat
    colnames(mat) <- paste(colnames(mat), " [", seq(1, ncol(mat)), "] ",
                         round(scores,3), sep="")
    
    # Color all the overlapping features as red and others as white
    colcode <-
      if (all(mat == 1)) c("firebrick2", "white") else c("white", "firebrick2")
    
    verbose("Generating top N overlap heatmap...\n")
    
    # Create the overlapping heatmap
    heatmap.2(
      x = mat,
      col = colcode,
      Colv = FALSE,
      dendrogram = "none",
      margins = c(10,10),
      cexRow = 0.7,
      cexCol = 0.7,
      cex.main = 0.8,
      key = FALSE,
      trace = "none",
      sepwidth = c(0.1,0.1),
      sepcolor = "grey90",
      colsep = seq_len(ncol(mat)),
      rowsep = seq_len(nrow(mat))
    )
    
    legend(
      "topleft",
      legend=c("Present","Absent"),
      fill=c("firebrick2","white"),
      bty="n"
    )
    
  }else{
    
    warning("Cannot plot overlap matrix for N = 1. ",
            "Please use a larger value for evaluating top N visualization.")
    
  }
}


