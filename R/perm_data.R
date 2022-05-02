#' Permutation Results
#'
#' A list of objects returned from \code{cadra_search()} when it runs
#' on simulated dataset (\code{ES = sim.ES}) with \code{input_score = ncol(sim.ES):1} and \code{method = "ks"}. 
#'
#' @docType data
#'
#' @usage data(perm.res)
#'
#' @format A list of objects returned from \code{cadra_search()}
#' containing top_N value, search_start value, permutation best scores across n_perm, permutation p-value and observed best score. 
#' The resulting list can be used to visualize the NULL permutation distribution by passing it to \code{\link[CaDrA]{permutation_plot}}
#' See \code{\link[CaDrA]{cadra_search}} and \code{\link[CaDrA]{permutation_plot}} more details. 
#'
#' @references Kartha VK, Kern JG, Sebastiani P, Zhang L,
#' Varelas X, Monti S (2017) CaDrA: A computational framework for performing 
#' candidate driver analyses using binary genomic features. 
#' (\href{https://www.biorxiv.org/content/early/2017/11/23/221846}{bioRxiv})
#'
#' @examples
#' 
#' # Load the pre-computed permutation results
#' data(perm.res)
#' 
#' # Plot the NULL permutation distribution against the observed best score
#' permutation_plot(perm.res)
#' 
"perm.res"
