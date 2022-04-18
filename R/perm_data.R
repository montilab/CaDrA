#' Permutation Results
#'
#' A list of objects returned from cadra_search() when it runs
#' on the \code{sim.ES} simulated dataset with ks method. 
#'
#' @docType data
#'
#' @usage data(perm.res)
#'
#' @format A list of objects returned by \code{cadra_search()}
#' containing parameters such as top N features, search start, permutation best scores, permutation p-values and observed best score. 
#' This resulting object can be used to visualize the null permutation distribution plot.
#' See \code{\link[CaDrA]{cadra_search}} for more information. 
#'
#' @references Kartha VK, Kern JG, Sebastiani P, Zhang L,
#' Varelas X, Monti S (2017) CaDrA: A computational framework for performing 
#' candidate driver analyses using binary genomic features. 
#' (\href{https://www.biorxiv.org/content/early/2017/11/23/221846}{bioRxiv})
#'
#' @examples
#' data(perm.res)
"perm.res"
