#' Pre-computed permutation results
#'
#' The permutation result returned from \code{CaDrA()} using simulated dataset
#' (\code{FS = sim_FS}) and pre-simulated input scores
#' \code{input_score = sim_Scores}
#'
#' @docType data
#'
#' @usage data(perm_res)
#'
#' @format A list of objects returned from \code{CaDrA} search. The result
#' contains a list of key parameters that are used for the
#' permutation-based testing, a vector of permuted best scores,
#' an observed best score, and a permutation p-value.
#'
#' To visualize the Empirical Null Distribution of permuted best scores,
#' just pass the resulting list to \code{\link[CaDrA]{permutation_plot}} function.
#' See \code{\link[CaDrA]{CaDrA}} and \code{\link[CaDrA]{permutation_plot}}
#' for more details.
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
