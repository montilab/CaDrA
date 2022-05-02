#' Top-N Results for Simulated Data (\code{sim.ES})
#'
#' A list of lists returned from topn_eval() using simulated dataset (\code{ES = sim.ES}),  
#' \code{input_score = ncol(sim.ES):1}, \code{top_N = 7}, \code{method = "ks"},
#' and \code{max_size = 10} as inputs to the function.
#' \code{Note:} \code{max_size} is set to 10 as we would like to account for the presence of 10 left-skewed 
#' (i.e. true positive or TP) features in the \code{sim.ES} dataset.
#' Over top_N = 7 feature searches, a list of ESet, along with its corresponding best score, 
#' and input_score are returned from each search. 
#'
#' @docType data
#'
#' @usage data(topn.list)
#'
#' @format A list of objects returned from \code{topn_eval()} function
#' containing ESet, corresponding best scores, and input_score of the top_N feature searches.
#' See \code{\link[CaDrA]{topn_eval}} for more information.
#'
#' @references Kartha VK, Kern JG, Sebastiani P, Zhang L,
#' Varelas X, Monti S (2017) CaDrA: A computational framework for performing 
#' candidate driver analyses using binary genomic features. 
#' (\href{https://www.biorxiv.org/content/early/2017/11/23/221846}{bioRxiv})
#'
#' @examples
#' 
#' # Load pre-computed Top-N list generated for sim.ES dataset
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
#' best_meta <- topn_best(topn.list)
#' 
#' # You can now visualize this result using the meta_plot() function 
#' meta_plot(best_meta)
#' 
#' # Or you can visualize the overlap of features across the top N (N=7) returned meta-features
#' # We do this by passing the topn.list object to the topn_plot() function
#' topn_plot(topn.list)
#' 
"topn.list"
