#' Pre-computed permutation results from simulated data \code{sim_FS} 
#'
#' The permutation result returned from \code{CaDrA()} using simulated dataset
#' (\code{FS = sim_FS}) and pre-simulated input scores
#' \code{input_score = sim_Scores}
#'
#' @docType data
#'
#' @usage data(perm_res)
#'
#' @format A list of objects returned from \code{CaDrA()} search. The result
#' contains a list of key parameters that are used to run the
#' permutation-based testing, a vector of permuted best scores,
#' an observed best score, and a permuted p-value.
#'
#' To visualize the Empirical Null Distribution of the permuted best scores over 
#' n_perm iterations, just pass the resulting list to \code{permutation_plot()}.
#' 
#' See \code{permutation_plot()} for more details.
#'
#' @references Kartha VK, Kern JG, Sebastiani P, Zhang L,
#' Varelas X, Monti S (2017) CaDrA: A computational framework for performing
#' candidate driver analyses using binary genomic features.
#' (\href{https://www.frontiersin.org/articles/10.3389/fgene.2019.00121/full}{Frontiers in Genetics})
#'
#' @examples
#'
#' # Load the pre-computed permutation results
#' data(perm_res)
#'
#' # Plot the NULL permutation distribution against the observed best score
#' permutation_plot(perm_res = perm_res)
#'
"perm_res"
