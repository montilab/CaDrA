#' Top-N result list for simulated binary genomic data
#'
#' The resulting list object returned by topn.eval() when run
#' on the \code{sim.ES} simulated dataset. All default parameters were used
#' when running topn.eval() except for max.size which was set to 10 
#' to account for the presence of 10 left-skewed (i.e. true positive or TP) 
#' features in the sim.ES dataset.
#'
#' @docType data
#'
#' @usage data(topn.list)
#'
#' @format A list object returned by \code{topn.eval}
#' containing the ExpressionSet objects (and corresponding meta-feature scores)
#' for each search run when starting with each of the top 7 ranked features.
#' Returned ExpressionSet columns (i.e. samples) are ranked by the user-defined ranking variable.  
#' See \code{\link[CaDrA]{topn.eval}} for more information.
#'
#' @references Kartha VK, Kern JG, Sebastiani P, Zhang L,
#' Varelas X, Monti S (2017) CaDrA: A computational framework for performing 
#' candidate driver analyses using binary genomic features. 
#' (\href{https://www.biorxiv.org/content/early/2017/11/23/221846}{bioRxiv})
#'
#' @examples
#' data(topn.list)
#' 
#' # To fetch each search result (Expression Set and score of the corresponding meta-feature)
#' # For N=1 (result when the search is initiated with the top-scoring starting feature)
#' topn.list[[1]]
#' 
#' # To fetch just the ExpressionSet object
#' topn.list[[1]]$ESet
#' 
#' # ExpressionSet for the search result when starting with the second-best feature
#' topn.list[[2]]$ESet
#' 
#' # Or we can find the result that had the best score over the top N (N=7) runs
#' best.meta <- topn.best(topn.list)
#' 
#' # You can now visualize this result using the meta.plot() function 
#' meta.plot(best.meta$ESet)
#' 
#' # Or you can visualize the overlap of features across the top N (N=7) returned meta-features
#' # We do this by passing the topn.list object to the topn.plot() function
#' topn.plot(topn.list)
"topn.list"
