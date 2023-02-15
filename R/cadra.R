#' CaDrA Search
#'
#' Perform permutation-based testings on a sample of permuted input scores 
#' using \code{candidate_search} as the main iterative function for each run.
#'
#' @param FS a SummarizedExperiment object containing binary features where
#' rows represent features of interest (e.g. genes, transcripts, exons, etc...)
#' and columns represent the samples.
#' @param input_score a vector of continuous scores representing a phenotypic
#' readout of interest such as protein expression, pathway activity, etc.
#' 
#' NOTE: \code{input_score} object must have names or labels that match the column
#' names of \code{FS} object.
#' @param method a character string specifies a scoring method that is
#' used in the search. There are 6 options: (\code{"ks_pval"} or \code{ks_score}
#' or \code{"wilcox_pval"} or \code{wilcox_score} or 
#' \code{"revealer"} (conditional mutual information from REVEALER) or
#' \code{"custom"} (a customized scoring method)). 
#' Default is \code{ks_pval}.
#' @param custom_function if method is \code{"custom"}, specifies
#' the name of the customized function here. Default is \code{NULL}.
#' 
#' NOTE: custom_function() must take \code{FS} and \code{input_score} 
#' as its input arguments, and its final result must return a vector of row-wise 
#' scores ordered from most significant to least significant where its labels or 
#' names matched the row names of \code{FS} object.
#' @param custom_parameters if method is \code{"custom"}, specifies a list of
#' additional arguments (excluding \code{FS} and \code{input_score}) to be 
#' passed to \code{custom_function}. Default is \code{NULL}.
#' @param alternative a character string specifies an alternative hypothesis
#' testing (\code{"two.sided"} or \code{"greater"} or \code{"less"}).
#' Default is \code{less} for left-skewed significance testing.
#' NOTE: this argument only apply to KS and Wilcoxon method
#' @param weight if method is \code{ks_score}, specifies a vector of weights 
#' to perform a weighted-KS testing. Default is \code{NULL}.
#' @param top_N an integer specifies the number of features to start the
#' search over. By default, it starts from the top best feature (top_N = 1).
#' NOTE: If \code{top_N} is provided, then \code{search_start} parameter
#' will be ignored.
#' @param search_start a list of character strings (separated by commas)
#' which specifies feature names within the FS object to start
#' the search with. If \code{search_start} is provided, then \code{top_N}
#' parameter will be ignored. Default is \code{NULL}.
#' @param search_method a character string specifies an algorithm to filter out
#' the best candidates (\code{"forward"} or \code{"both"}). Default is
#' \code{both} (i.e., backward and forward).
#' @param max_size an integer specifies a maximum size that a meta-feature can
#' extend to do for a given search. Default is \code{7}.
#' @param n_perm an integer specifies the number of permutations to perform.
#' Default is \code{1000}.
#' @param smooth a logical value indicates whether or not to smooth the p-value
#' calculation to avoid p-value of 0. Default is \code{TRUE}.
#' @param obs_best_score a numeric value corresponding to the best observed
#' score. This value is used to compare against the permuted best scores.
#' Default is \code{NULL}. If set to NULL, we will compute the observed
#' best score based on the given parameters.
#' @param smooth a logical value indicates whether or not to smooth the p-value
#' calculation to avoid p-value of 0. Default is \code{TRUE}.
#' @param plot a logical value indicates whether or not to plot the empirical
#' null distribution of the permuted best scores. Default is \code{TRUE}.
#' @param ncores an integer specifies the number of cores to perform
#' parallelization for permutation-based testing. Default is \code{1}.
#' @param cache_path a full path uses to cache the permuted best scores.
#' We recycle these scores instead of re-computing them to save time.
#' Default is \code{NULL}. If NULL, the cache path is set to \code{~/.Rcache}
#' for future loading.
#' @param verbose a logical value indicates whether or not to print the
#' diagnostic messages. Default is \code{FALSE}.
#'
#' @return a list of key parameters that are used to cache the result of
#' permutation-based testing, a vector of permuted best scores for a given
#' \code{n_perm}, an observed best score, and a permutation p-value.
#'
#' @examples
#'
#' # Load pre-computed feature set
#' data(sim_FS)
#'
#' # Load pre-computed input-score
#' data(sim_Scores)
#'
#' # Define additional parameters and start the function
#' cadra_result <- CaDrA(
#'   FS = sim_FS, input_score = sim_Scores, method = "ks_pval", 
#'   weight = NULL, alternative = "less", top_N = 1,
#'   search_start = NULL, search_method = "both", max_size = 7, 
#'   n_perm = 10, plot = FALSE, smooth = TRUE, obs_best_score = NULL,
#'   ncores = 1, cache_path = NULL
#' )
#'
#' @export
#' @import R.cache doParallel ggplot2 plyr methods
#'
#' @author Reina Chau
#'
CaDrA <- function(
    FS,
    input_score,
    method = c("ks_pval", "ks_score", "wilcox_pval", "wilcox_score", "revealer", "custom"),
    custom_function = NULL,
    custom_parameters = NULL,
    alternative = c("less", "greater", "two.sided"),
    weight = NULL,
    top_N = 1,
    search_start = NULL,
    search_method = c("both", "forward"),
    max_size = 7,
    n_perm = 1000,
    obs_best_score = NULL,
    smooth = TRUE,
    plot = TRUE,
    ncores = 1,
    cache_path = NULL,
    verbose = FALSE
){

  # Set up verbose option
  options(verbose = verbose)

  # Match arguments
  method <- match.arg(method)
  alternative <- match.arg(alternative)
  search_method <- match.arg(search_method)

  # Check n_perm
  stopifnot("invalid number of permutations (nperm)"=
              (length(n_perm)==1 && !is.na(n_perm) &&
                 is.numeric(n_perm) && n_perm > 0) )

  # Check ncores
  stopifnot("invalid number of CPU cores (ncores)"=
              (length(ncores)==1 && !is.na(ncores) &&
                 is.numeric(ncores) && ncores > 0) )

  ####### CACHE CHECKING #######
  if(!is.null(cache_path)){
    R.cache::setCacheRootPath(cache_path)
    message("Using provided cache root path: ", cache_path, "")
  } else{
    R.cache::setCacheRootPath()
    message("Setting cache root path as: ", getCacheRootPath(), "\n")
  }

  # Define the key for each cached result
  key <- list(FS = FS,
              input_score = if(method %in% c("revealer", "custom"))
              { input_score } else { NULL },
              method = method,
              custom_function = custom_function,
              custom_parameters = custom_parameters,
              alternative = alternative,
              weight = weight,
              top_N = top_N,
              search_start = search_start,
              search_method = search_method,
              max_size = max_size)

  # Load perm_best_scores with the given key parameters
  perm_best_scores <- R.cache::loadCache(key)

  # Start the 'clock' to see how long the process takes
  ptm <- proc.time()

  # Check if, given the dataset and search-specific parameters,
  # there is already a cached null distribution available
  n_perm <-  as.integer(n_perm)

  if(!is.null(perm_best_scores) & (length(perm_best_scores) >= n_perm)){

    if(length(perm_best_scores) == n_perm){
      message("Found ", length(perm_best_scores),
              " permutated scores for the specified dataset",
              " and search parameters in cache path\n")
      message("LOADING PERMUTATED SCORES FROM CACHE\n")
    }else{
      message("n_perm is set to ", n_perm, " but found ",
              length(perm_best_scores),
              " permutated scores for the specified dataset",
              " and search parameters in cache path\n")
      message("LOADING LARGER PERMUTATED SCORES FROM CACHE\n")
    }

  }else{

    if(is.null(perm_best_scores)){
      verbose("No permutated scores for the specified dataset and ",
              "search parameters were found in cache path\n")
      verbose("BEGINNING PERMUTATION-BASED TESTINGS\n")
    }else if (length(perm_best_scores) < n_perm) {
      verbose("n_perm is set to ", n_perm, " but found only ",
              length(perm_best_scores),
              " permutated scores for the specified dataset",
              " and search parameters in cache path\n")
      verbose("RE-COMPUTE PERMUTATION-BASED TESTINGS WITH LARGER NUMBER OF PERMUTATIONS\n")
    }

    #######################################################################

    # Check ncores
    ncores <-  as.integer(ncores)

    # Sets up the parallel backend which will be utilized by Plyr.
    parallel <- FALSE
    progress <- "text"

    if(ncores > 1){
      doParallel::registerDoParallel(cores = ncores)
      parallel <- TRUE
      progress <- "none"
      verbose("Running tests in parallel...")
    }

    # Generate matrix of permuted input_score
    perm_labels_matrix <- generate_permutations(
      input_score = input_score, 
      n_perm = n_perm
    )
    
    # Run permutation-based testing
    perm_best_scores_l <- plyr::alply(
      perm_labels_matrix,
      1,
      function(x){
        perm_input_score <- x
        names(perm_input_score) <- colnames(perm_labels_matrix)
        best_score <- candidate_search(
          FS = FS,
          input_score = perm_input_score,
          method = method,
          custom_function = custom_function,
          custom_parameters = custom_parameters,
          alternative = alternative,
          weight = weight,
          top_N = top_N,
          search_start = search_start,
          search_method = search_method,
          max_size = max_size,
          best_score_only = TRUE,
          do_plot = FALSE,
          do_check = FALSE,
          verbose = FALSE
        ) 
        return(best_score)
      },
      .parallel = parallel,
      .progress = progress)
    
    # Extract the permuted best scores
    perm_best_scores <- lapply(
      seq_along(perm_best_scores_l), 
      function(l){ perm_best_scores_l[[l]] }) |> unlist()
    
    # Save computed scores to cache
    verbose("Saving to cache...\n")
    R.cache::saveCache(perm_best_scores, key=key, comment="null_scores()")
    
  } # end caching else statement block
  
  # Return to using just a single core
  doParallel::registerDoParallel(cores = 1)
  
  verbose("FINISHED\n")
  verbose("Time elapsed: ", round((proc.time()-ptm)[3]/60, 2), " mins \n\n")

  #########################################################################

  if(is.null(obs_best_score)){

    verbose("Computing observed best score...\n\n")

    obs_best_score <- candidate_search(
      FS = FS,
      input_score = input_score,
      method = method,
      custom_function = custom_function,
      custom_parameters = custom_parameters,
      alternative = alternative,
      weight = weight,
      top_N = top_N,
      search_start = search_start,
      search_method = search_method,
      max_size = max_size,
      best_score_only = TRUE,
      do_plot = FALSE,
      do_check = FALSE,
      verbose = FALSE
    ) |> unlist()

  }else{

    # Check obs_best_score
    stopifnot("invalid observed best score (obs_best_score)"=
                (length(obs_best_score)==1 && !is.na(obs_best_score) &&
                   is.numeric(obs_best_score)))
    
    verbose("Using provided value of observed best score...\n\n")
    obs_best_score <- as.numeric(obs_best_score)
    
  }

  verbose("Observed score: ", obs_best_score, "\n")

  ########### PERMUTATION P-VALUE COMPUTATION ############

  #Add a smoothing factor of 1 if smooth is specified
  #This is just to not return a p-value of 0
  c <- 0
  if(smooth) c <- 1

  perm_pval <- (sum(perm_best_scores > obs_best_score) + c)/
    (length(perm_best_scores) + c)

  verbose("Permutation p-value: ", perm_pval, "\n")
  verbose("Number of permutations: ", length(perm_best_scores), "\n")
  
  ########### END PERMUTATION P-VALUE COMPUTATION ############
  perm_res <- list(
    key = key,
    perm_best_scores = perm_best_scores,
    obs_best_score = obs_best_score,
    perm_pval = perm_pval
  )

  # If plot = TRUE, produce the permutation plot
  if(plot == TRUE){
    permutation_plot(perm_res = perm_res)
  }

  return(perm_res)

}
