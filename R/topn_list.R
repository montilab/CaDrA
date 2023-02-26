#' Top-N Results for Simulated Data (\code{sim_FS})
#'
#' A list of objects returned from \code{candidate_search()} using simulated dataset
#' \code{FS = sim_FS}, \code{input_score = sim_Scores}, \code{top_N = 7},
#' \code{method = "ks_pval"}, \code{alternative = "less"}, 
#' \code{search_method = "both"}, \code{max_size = 10},
#' and \code{best_score_only = FALSE} as inputs to the function.
#' \code{NOTE:} \code{max_size} is set to 10 as we would like to account
#' for the presence of 10 left-skewed (i.e. true positive or TP) features 
#' in \code{sim_FS} dataset.
#' 
#' Over top_N = 7 feature searches, a set of meta-features in form of SummarizedExperiment 
#' object, along with a vector of observed input scores and its corresponding best score 
#' are returned from each search.
#'
#' @docType data
#'
#' @usage data(topn_list)
#'
#' @format A list of objects returned from \code{candidate_search()} including 
#' a set of meta-features in form of SummarizedExperiment objects, 
#' its observed input_score, and corresponding best score pertaining to each 
#' top N feature searches.
#' 
#' See \code{\link[CaDrA]{candidate_search}} for more information.
#'
#' @references Kartha VK, Kern JG, Sebastiani P, Zhang L,
#' Varelas X, Monti S (2019) CaDrA: A computational framework for performing
#' candidate driver analyses using binary genomic features.
#' (\href{https://www.frontiersin.org/articles/10.3389/fgene.2019.00121/full}{Frontiers in Genetics})
#'
#' @examples
#'
#' # Load pre-computed Top-N list generated for sim_FS and sim_Scores dataset
#' data(topn_list)
#'
#' # Fetch the first meta-feature
#' topn_list[[1]]$feature_set
#'
#' # Fetch the second meta-feature
#' topn_list[[2]]$feature_set
#'
#' # Retrieve the meta-feature with the best score among top_N = 7 runs
#' topn_best_meta <- topn_best(topn_list = topn_list)
#'
#' # Visualize the best meta-feature using meta_plot() function
#' meta_plot(topn_best_list = topn_best_meta)
#'
#' # Visualize overlap of meta-features across top_N = 7 using topn_plot() function
#' topn_plot(topn_list = topn_list)
#'
"topn_list"
