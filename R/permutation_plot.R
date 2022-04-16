
#' Permutation Plot
#' 
#' Plots the Emperical Null distribution
#' @param permutation_results a list of lists, where each list entry is one that is returned by the candidate search run for a given starting index. This is computed within and can be returned by the topn.eval() function.
#' 
#' @return a density plot
#' @examples
#' # Load pre-computed Top-N list generated for sim.ES dataset
#' data(topn.list)
#' topn_plot(topn.list)
#' 
#' @export
#' @import gplots
#' @importFrom graphics legend
permutation_plot <- function(permutation_results){
  
  ## Extract values from the the permutation results
  top_N <- permutation_results[["top_N"]]
  search_start <- permutation_results[["search_start"]]
  perm_best_scores <- permutation_results[["perm_best_scores"]]
  perm_pval <- permutation_results[["perm_pval"]]
  obs_best_score <- permutation_results[["obs_best_score"]]
  
  plot_title <- paste("Emperical Null distribution (N = ", length(perm_best_scores), ")\n Permutation p-val <= ", round(perm_pval, 5), "\nBest observed score: ", round(obs_best_score, 5), sep="")
  
  if(!is.null(top_N)){
    plot_title <- paste(plot_title,"\n Top N: ", top_N, sep="")
  }else{
    plot_title <- paste(plot_title,"\n Seed: ", search_start, sep="")
  }
  
  #Here, let us plot the absolute values of the permutation p-values, for simplicity
  #You only consider absolute values when calculating the permutation p-values
  #Katia: adding ".data" to avoid a warning during check:
  # no visible binding for global variable 
  g <- ggplot(data = data.frame("x" = perm_best_scores), aes(x = .data$x)) +
    geom_histogram(fill = "black", color = "gray") +
    theme_classic() +
    theme(
      axis.line.x=element_line(color = "black"),
      axis.line.y=element_line(color = "black")
    )
  
  g <- g + geom_vline(xintercept = obs_best_score, linetype = "longdash", size = 1.5, colour = "red") +
    labs(
      title = plot_title,
      x = "Score",
      y = "Count"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0))
  
  print(g)
  
}


