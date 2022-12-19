#' CaDrA Search
#' 
#' Performs a permutation-based testing on sample permutation of observed input 
#' scores using candidate_search() as the main iterative function for each run.
#'
#' @param ES an expression set of binary features (required). It must be a 
#' \code{BioBase expressionSet} object. The rownames of the expression set must 
#' contain unique features which are used to search for best features.   
#' @param input_score a vector of continuous values of a response of 
#' interest (required). The \code{input_score} object must have names or 
#' labels that match the colnames of the expression matrix.
#' @param method a character string specifies a scoring method that is 
#' used in the search. There are 4 options: (\code{"ks"} or \code{"wilcox"} or 
#' \code{"revealer"} (conditional mutual information from REVEALER) or 
#' \code{"custom"} (a user customized scoring method)). Default is \code{ks}.
#' @param custom_function if method is \code{"custom"}, specifies the 
#' customized function here. Default is \code{NULL}.
#' @param custom_parameters if method is \code{"custom"}, specifies a list of 
#' arguments to be passed to the custom_function(). Default is \code{NULL}.
#' @param alternative a character string specifies an alternative hypothesis 
#' testing (\code{"two.sided"} or \code{"greater"} or \code{"less"}). 
#' Default is \code{less} for left-skewed significance testing.
#' @param metric a character string specifies a metric to search for 
#' best features. \code{"pval"} or \code{"stat"} may be used which 
#' corresponding to p-value or score statistic. Default is \code{pval}. 
#' NOTE: \code{Revealer} method only utilized score statistics values 
#' (no p-values).
#' @param weights if method is \code{ks}, specifies a vector of weights 
#' will perform a weighted-KS testing. Default is \code{NULL}.   
#' @param top_N an integer specifies the number of features to start the 
#' search over, starting from the top 'N' features in each case. If \code{top_N} 
#' is provided, then \code{search_start} parameter will be ignored. Default is 
#' \code{1}.
#' @param search_start a list of character strings (separated by commas) 
#' which specifies feature names within the expression set object to start 
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
#' score or p-value and later use to compare against permutation based scores or
#' p-values. Default is \code{NULL}. If set to NULL, we will compute the observed 
#' best score based on the given \code{input_score} and \code{ES} variables.
#' @param plot a logical value indicates whether or not to plot the empirical 
#' null distribution with the observed scores and permutation p-values. Default is
#' \code{TRUE}.
#' @param ncores an integer specifies the number of cores to perform 
#' parallelization for permutation testing. Default is \code{1}.
#' @param cache_path a full path uses to cache permutation-based score 
#' distributions. If the permutation for a given ES and its dependent search 
#' variables such as 'top_N' exist, we recycle these values instead of 
#' re-computing them to save time. Default is \code{NULL}. If NULL, the 
#' cache path is set to \code{~/.Rcache} for future loading.
#' @param verbose a logical value indicates whether or not to print the 
#' diagnostic messages. Default is \code{FALSE}. 
#' 
#' @return a list of key parameters that were used to cache the permutation 
#' test, permutation best scores for n_perm, a permutation p-value, and 
#' observed best score
#' @examples
#'\donttest{
#'
#' # Load R library
#' library(Biobase)
#'
#' # Load pre-computed expression set
#' data(sim.ES)
#' 
#' # Set seed
#' set.seed(123)
#' 
#' # Provide a vector of continuous scores for a target profile with names 
#' # to each score value 
#' input_score = rnorm(n = ncol(sim.ES))
#' names(input_score) <- colnames(sim.ES)
#' 
#' # Define additional parameters and start the function
#' # This function takes some time to run
#' cadra_result <- CaDrA(
#'   ES = sim.ES, input_score = input_score, method = "ks", weights = NULL,
#'   alternative = "less", metric = "pval", top_N = 1, 
#'   search_start = NULL, search_method = "both", max_size = 7, n_perm = 1000, 
#'   plot = TRUE, smooth = TRUE, obs_best_score = NULL, 
#'   ncores = 1, cache_path = NULL
#' )
#'  
#'}
#' 
#' @export
#' @import Biobase R.cache doParallel ggplot2 plyr methods
#'
#' @author Reina Chau
#' 
CaDrA <- function(
    ES,
    input_score,
    method = c("ks","wilcox","revealer", "custom"),
    custom_function = NULL,
    custom_parameters = NULL,
    alternative = c("less", "greater", "two.sided"),
    metric = c("pval", "stat"),
    weights = NULL,
    top_N = 1,
    search_start = NULL,
    search_method = c("both", "forward"),
    max_size = 7,  
    n_perm = 1000,
    smooth = TRUE,
    obs_best_score = NULL,
    plot = TRUE,
    ncores = 1,
    cache_path = NULL,
    verbose = FALSE
){
  
  # Set up verbose option
  options(verbose = verbose)
  
  method <- match.arg(method)  
  alternative <- match.arg(alternative)  
  metric <- match.arg(metric)  
  search_method <- match.arg(search_method)  
  
  # Stop if input values are not valid
  cadra_check_input( ES, input_score, n_perm, ncores)
  

  ## Samples to keep based on the overlap between the two inputs
  overlap <- intersect(names(input_score), Biobase::sampleNames(ES))
  ES <- ES[,overlap]
  input_score <- input_score[overlap]
  
  # Check if the dataset has any all 0 or 1 features 
  # (these are to be removed since they are not informative)
  if(any(rowSums(exprs(ES)) == 0) || 
     any(rowSums(exprs(ES)) == ncol(exprs(ES)))){
    warning("Provided dataset has features that are either all 0s or 1s. ",
            "These features will be removed from the computation.")
    ES <- ES[!(rowSums(exprs(ES))==0 | rowSums(exprs(ES)) == ncol(exprs(ES))),]
    
    # Make sure matrix is not empty after removing uninformative features
    stopifnot("There are no more features remained for downsteam computation."=
                (nrow(exprs(ES)) != 0))
  }
 
  # Select the appropriate method to compute scores based on 
  # skewness of a given binary matrix
  mat <- exprs(ES)
  
  # Compute the row-wise scoring
  s <- switch(
    method,
    ks = ks_genescore_mat(                # Kolmogorov-Smirnov Scoring Method    
        mat = mat,
        alternative = alternative, 
        weights = weights),
    wilcox = wilcox_genescore_mat(        # Wicox Scoring Method
      mat = mat,
      alternative = alternative,
      ranks = NULL),
    revealer = revealer_genescore_mat(    # Revealer Scoring Method
      mat = mat,                                   
      input_score = input_score,      
      seed_names = NULL,
      target_match = "positive",
      assoc_metric = "IC"),
    custom = custom_genescore_mat(        # Custom Scoring Method
      mat = mat,
      input_score = input_score,      
      custom_function = custom_function,
      custom_parameters = custom_parameters)
  )
  
  # Check if the returning result has one or two columns: 
  # score or p_value or both
  if(ncol(s) == 1 && colnames(s)[1] == "score" && metric == "pval"){
      warning("metric = 'pval' is specified but the ", method, 
              " method only returns score values. ",
              "Thus, using 'stat' as metric to search for best features.")
      metric <- "stat"
  }else if(ncol(s) == 1 && colnames(s)[1] == "p_value" & metric == "stat"){
      warning("metric requested is 'stat' but the ", method, 
              "method only returns p-values. ",
              "Thus, using 'pval' as metric to search for best features.")
      metric <- "pval"
  }
  
  ####### CACHE CHECKING #######
  options(verbose = verbose)

  if(!is.null(cache_path)){
    message("Using provided cache root path: ", cache_path, "")
    setCacheRootPath(cache_path)
  } else{
    setCacheRootPath()
    message("Setting cache root path as: ", getCacheRootPath(), "\n")
  }
  
  # We use the ES, top N (or search_start), score metric, 
  # scoring method as the key for each cached result  
  key <- list(ES=ES, 
              input_score=if(method %in% c("revealer", "custom"))
                { input_score }else { NULL },
              method=method, 
              custom_function=custom_function, 
              custom_parameters=custom_parameters, 
              alternative=alternative, 
              metric=metric, 
              weights=weights, 
              top_N=top_N, 
              search_start=search_start, 
              search_method=search_method, 
              max_size=max_size, 
              best_score_only=TRUE)
  
  perm_best_scores <- loadCache(key)
  
  #Start the 'clock' to see how long the process takes
  ptm <- proc.time()
  
  # Check if, given the dataset and search-specific parameters,
  # there is already a cached null distribution available
  n_perms <-  as.integer(n_perm)
  if (!is.null(perm_best_scores) & (length(perm_best_scores) >= n_perm)){
    
    if(length(perm_best_scores) == n_perm){
      message("Found ", length(perm_best_scores), 
              " cached permutation-based scores ",
              "for the specified dataset and search parameters...\n")      
      message("LOADING PERMUTATION SCORES FROM CACHE\n")
    }else{
      message("n_perm is set to ", n_perm, " but found ", 
              length(perm_best_scores), 
              " cached permutation-based scores ",
              "for the specified dataset and search parameters...\n")
      message("LOADING LARGER PERMUTATION SCORES FROM CACHE\n")
    }
    
  } else{
    
    if (is.null(perm_best_scores)){
      message("No permutation scores for the specified dataset and ",
              "search parameters were found in cache path...\n")
      message("BEGINNING PERMUTATION-BASED TESTINGS\n")
    } else if (length(perm_best_scores) < n_perm) {
      message("n_perm is set to ", n_perm, " but found only ", 
              length(perm_best_scores), 
              " cached permutation-based scores for the specified dataset ",
              "and search parameters...\n")
      message("RE-COMPUTE PERMUTATION-BASED TESTINGS WITH LARGER NUMBER OF PERMUTATIONS")
    }
    
    #######################################################################
    
    # Sets up the parallel backend which will be utilized by Plyr.
    parallel <- FALSE
    progress <- "text"
    
    if(ncores > 1){
      registerDoParallel(cores = ncores)
      parallel <- TRUE
      progress <- "none"
      verbose("Running tests in parallel...")
    } 
    
    message("Using ", ncores, " core(s)...")
    
    # Generate matrix of permutated input_score  
    perm_labels_matrix <- generate_permutations(ord=input_score, 
                                                n_perm = n_perm, 
                                                verbose = FALSE)
    
    verbose("Computing permutation-based scores for N = ", n_perm, "...\n")
    
    perm_best_scores <- plyr::alply(
      perm_labels_matrix, 
      1, 
      function(x){ 
        perm_input_score <- x;
        names(perm_input_score) <- names(input_score); 
        candidate_search(ES=ES, 
                         input_score=perm_input_score, 
                         method=method, 
                         custom_function=custom_function, 
                         custom_parameters=custom_parameters, 
                         alternative=alternative, 
                         metric=metric, 
                         weights=weights, 
                         top_N=top_N, 
                         search_start=search_start, 
                         search_method=search_method, 
                         max_size=max_size, 
                         best_score_only=TRUE, 
                         do_plot = FALSE, 
                         verbose = FALSE) }, 
      .parallel=parallel, 
      .progress=progress) %>% unlist()
    
    #Save computed scores to cache 
    message("Saving to cache...\n")
    saveCache(perm_best_scores, key=key, comment="null_ks()")
    
  } # end caching else statement block
  
  registerDoParallel(cores = 1) #Return to using just a single core
  
  message("FINISHED\n")
  message("Time elapsed: ", round((proc.time()-ptm)[3]/60,2), " mins \n\n")
  
  #########################################################################
  
  if(is.null(obs_best_score)){
    
    message("Computing observed best score...\n\n")
    
    obs_best_score <- candidate_search(
      ES = ES,
      input_score = input_score, 
      method = method,
      custom_function = custom_function,
      custom_parameters = custom_parameters,
      alternative = alternative,
      metric = metric,
      weights = weights,
      top_N = top_N,
      search_start = search_start,
      search_method = search_method,
      max_size = max_size,
      best_score_only = TRUE,
      do_plot = FALSE,
      verbose = FALSE
    ) %>% unlist()
    
  } else{
    
    message("Using provided value of observed best score...\n\n")
    
  }
  
  message("Observed score: ", obs_best_score, "\n")
  
  ########### PERMUTATION P-VALUE COMPUTATION ############
  message("Number of permutation-based scores being considered: ", 
          length(perm_best_scores), "\n")
  
  #Add a smoothening factor of 1 if smooth is specified
  #This is just to not return a p-value of 0
  c <- 0
  if(smooth) c <- 1
  
  if(metric == "pval"){
    
    # Use negative log transform of returned search score 
    # (either computed above, or passed to the null_ks function 
    # if previously computed)
    obs_best_score <- -(log(obs_best_score))
    
    # Use negative log transform on the permuted scores 
    # (either computed above or loaded from Cache)
    # NOTE: there is a very small chance some signed observed scores (p-values) 
    # are anti-correlated (meaning negative)
    # To avoid NaNs, we remove these. 
    # Keep in mind this is uncommon and will contribute very few 
    # permutations (n<10) if running N=1000
    perm_best_scores <- perm_best_scores[perm_best_scores > 0]
    perm_best_scores <- -(log(perm_best_scores))
    
  }
  
  perm_pval <- (sum(perm_best_scores > obs_best_score) + c)/
    (length(perm_best_scores) + c) 
  
  message("Permutation p-value: ", perm_pval, "\n\n")
  
  ########### END PERMUTATION P-VALUE COMPUTATION ############
  # Plot result
  
  
  if(plot == TRUE)cadra_plot( top_N,
                              search_start,
                              obs_best_score,
                              perm_pval,
                              perm_best_scores)
  
  perm_res <- list(
    key = key,
    perm_best_scores = perm_best_scores, 
    perm_pval = perm_pval, 
    obs_best_score = obs_best_score
  )
  
  return(perm_res)
  
}


#' Random permutation matrix generator
#' 
#' Produces a random permutation rank matrix given a vector of values
#' @param ord vector to be permuted. This determines the number of columns 
#' in the permutation matrix
#' @param n_perm number of permutations to generate. This determines 
#' the number of rows in the permutation matrix
#' @param verbose a logical value indicates whether or not to print the 
#' diagnostic messages. Default is \code{FALSE}. 
#' @return A row matrix of permuted values (i.e. ranks) where each row is a 
#' single permutation result
#' 
#' @examples
#' 
#' # Load pre-computed input score
#' data(TAZYAP_BRCA_ACTIVITY)
#' input_score = TAZYAP_BRCA_ACTIVITY
#' 
#' # Set seed
#' set.seed(123)
#' 
#' # Define number of permutations
#' n_perm = 1000
#'  
#' # Define additional parameters and start the function
#' perm_labels_matrix <- generate_permutations(ord=input_score, n_perm=n_perm)
#'   
#' @export
generate_permutations <- function(
    ord,                  # These are the sample orderings to be permuted
    n_perm,              # Number of permutations to produce
    verbose = FALSE
){
  
  options(verbose = verbose)
  
  # Get number of samples
  m  <- length(ord)
  
  # Create permutation matrix
  perm <- matrix(NA, nrow=n_perm, ncol=length(ord))
  
  verbose("Generating ", n_perm," permuted sample observed input scores...")
  
  # Sample the input scores
  for (i in seq_len(n_perm) ) perm[i,] <- sample(ord, m)
  
  verbose("Are all generated permutations unique?..")
  verbose(nrow(perm)==nrow(unique.matrix(perm)))
  
  stopifnot("invalid nperm value" = 
           (nrow(perm) == nrow(unique.matrix(perm))) )
  
  return(perm) 
  
}
