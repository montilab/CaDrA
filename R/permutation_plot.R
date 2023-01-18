
#' Permutation Best Scores Plot
#'
#' Plot the Empirical Null Distribution of Permutation Best Scores Based on
#' a given number of permutations (\code{n_perm})
#'
#' @param perm_res a list of objects returned from CaDrA() function
#' using the simulated dataset (\code{sim.ES}) and random generated
#' \code{input_score = rnorm(n = ncol(sim.ES))} with set.seed(123).
#'
#' @return a density plot
#' @examples
#'
#' # Load pre-computed permutation results
#' data(perm_res)
#'
#' # Plot the permutation results
#' permutation_plot(perm_res)
#'
#' @export
#' @import gplots
#' @importFrom graphics legend
permutation_plot <- function(perm_res){

  ## Extract values from the the permutation results
  top_N <- perm_res[["key"]][["top_N"]]
  search_start <- perm_res[["key"]][["search_start"]]
  perm_best_scores <- perm_res[["perm_best_scores"]]
  perm_pval <- perm_res[["perm_pval"]]
  obs_best_score <- perm_res[["obs_best_score"]]

  plot_title <- paste("Emperical Null distribution (N = ",
                      length(perm_best_scores),
                      ")\n Permutation p-val <= ",
                      round(perm_pval, 5),
                      "\nBest observed score: ",
                      round(obs_best_score, 5), sep="")

  if(!is.null(top_N)){
    plot_title <- paste(plot_title,"\n Top N: ", top_N, sep="")
  }else{
    plot_title <- paste(plot_title,"\n Seed: ", search_start, sep="")
  }

  # Here, let us plot the absolute values of the permutation p-values,
  # for simplicity
  #You only consider absolute values when calculating the permutation p-values
  g <- ggplot(data = data.frame("x" = perm_best_scores), aes(x = .data$x)) +
    geom_histogram(fill = "black", color = "gray") +
    theme_classic() +
    theme(
      axis.line.x=element_line(color = "black"),
      axis.line.y=element_line(color = "black")
    )

  g <- g + geom_vline(xintercept = obs_best_score, linetype = "longdash",
                      size = 1.5, colour = "red") +
    labs(
      title = plot_title,
      x = "Score",
      y = "Count"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0))

  g

}


