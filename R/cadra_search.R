
#' CaDrA Search
#' 
#' Performs a permutation-based testing using results from \code{candidate_search()} function.
#'
#' @param ES an expression set of binary features (required). It must be a BioBase expressionSet object. The rownames of the expression set must contain unique features which are used to search for best features.   
#' @param input_score a vector of continuous values for a target profile (required). The \code{input_score} object must have names or labels that match the colnames of the expression matrix.
#' @param method a character string specifies a method to compute the score for each feature (\code{"ks"} or \code{"wilcox"} or \code{"revealer"} (conditional mutual information from REVEALER) or \code{"custom"} (a customized method)). Default is \code{ks}.
#' @param custom_function if method is \code{"custom"}, specifies the customized function here. Default is \code{NULL}.
#' @param custom_parameters if method is \code{"custom"}, specifies a list of arguments to be passed to the custom_function(). Default is \code{NULL}.
#' @param alternative a character string specifies an alternative hypothesis testing (\code{"two.sided"} or \code{"greater"} or \code{"less"}). Default is \code{less} for left-skewed significance testing.
#' @param metric a character string specifies a metric to use for candidate search criteria. \code{"pval"} or \code{"stat"} may be used, corresponding to the score p-value or statistic. Default is \code{pval}.
#' @param weights a vector of weights use to perform a weighted-KS testing. Default is \code{NULL}.   
#' @param target_match a direction of target matching (\code{"negative"} or \code{"positive"}) from REVEALER. Use \code{"positive"} to match the higher values of the target, \code{"negative"} to match the lower values. Default is \code{positive}. 
#' @param top_N an integer specifies the number of features to start the search over, starting from the top 'N' features in each case. Default is \code{NULL}.
#' @param search_start an integer specifies an index within the expression set object of which feature to start the candidate search with. Default is \code{NULL}. If NULL, then the search starts with the top ranked feature. If an integer is specified (N, where N < nrow(ES)), the search starts with the Nth best feature. If a string is specified, the search starts with the feature with this name (must be a valid rowname in ES)
#' @param search_method a character string specifies a method to filter out the best candidates (\code{"forward"} or \code{"both"}). Default is \code{both} (backward and forward).
#' @param max_size an integer specifies a maximum size that a meta-feature can extend to do for a given search. Default is \code{7}.
#' @param n_perm an integer specifies the number of permutations to perform. Default is \code{1000}.
#' @param seed a seed set for permutation. Default is \code{123}.
#' @param smooth a logical value indicates whether or not to smooth the p-value calculation to avoid p-value of 0. Default is \code{TRUE}.
#' @param obs_best_score a numeric value corresponding to the observed (best) candidate search score to use for permutation based p-value computation. Default is \code{NULL}. If set to NULL, we compute the observed score given the \code{input_score} and \code{ES} variables.
#' @param plot a logical value indicates whether or not to plot the empirical null distribution with the observed score and permutation p-value. Default is \code{TRUE}.
#' @param ncores an integer specifies the number of cores to perform parallelization for permutation testing. Default is \code{1}.
#' @param cache_path a full path uses to cache permutation-based score distributions. If the permutation for a given ES and its dependent search variables such as 'N' exist, we recycle these values instead of re-computing them to save time. Default is \code{NULL}. If NULL, the cache path is set to \code{~/.Rcache} for future loading.
#' @param return_perm_pval a logical value indicates whether or not to return the permutation-based p-value computed by the function. Default is \code{TRUE}.
#' 
#' @return If \code{return_perm_pval} is set to \code{TRUE}, it will return a permutation p-value.
#' @examples
#' 
#' # Load R library
#' library(Biobase)
#'
#' # Load pre-computed expression set
#' data(sim.ES)
#' 
#' # set seed
#' set.seed(123)
#' 
#' # Provide a vector of continuous scores for a target profile with names to each score value 
#' input_score = rnorm(n = ncol(sim.ES))
#' names(input_score) <- colnames(sim.ES)
#' 
#' # Define additional parameters and start the function
#' # Not run as this would take some time to run
#' # candidate_search_result <- cadra_search(
#' #  ES = sim.ES, input_score = input_score, method = "ks", weights = NULL,
#' #  alternative = "less", metric = "pval", top_N = NULL, 
#' #  search_start = NULL, search_method = "both", max_size = 7, n_perm = 1000, 
#' #  seed = 123, plot = TRUE, smooth = TRUE, obs_best_score = NULL, 
#' #  ncores = 1, cache_path = NULL, return_perm_pval = TRUE
#' #)
#' 
#' @export
#' @import Biobase R.cache doParallel ggplot2 plyr
#'
#' @author Reina Chau
#' 
cadra_search <- function(
  ES,
  input_score,
  method = "ks",
  custom_function = NULL,
  custom_parameters = NULL,
  alternative = "less",
  metric = "pval",
  weights = NULL,
  target_match = "positive",
  top_N = NULL,
  search_start = NULL,
  search_method = "both",
  max_size = 7,  
  n_perm = 1000,
  seed = 123,
  smooth = TRUE,
  obs_best_score = NULL,
  plot = TRUE,
  ncores = 1,
  cache_path = NULL,
  return_perm_pval = TRUE
){
  
  # Check if the ES is provided and is a BioBase ExpressionSet object
  if(length(ES) == 0 || class(ES)[1] != "ExpressionSet") 
    stop("'ES' must be an ExpressionSet class argument (required).")
  
  # Check if the dataset has only binary 0 or 1 values 
  if(!all(exprs(ES) %in% c(0,1))){
    stop("The expression matrix (ES) must contain only binary values with no NAs.\n")
  }
  
  # Make sure the input ES has rownames for features tracking
  if(is.null(rownames(ES)))
    stop("The ES object does not have rownames or featureData to track the features by. Please provide unique features or rownames for the expression matrix.\n")
  
  # Check input_score is provided and is a continuous values with no NAs
  if(length(input_score) == 0 || any(!is.numeric(input_score)) || any(is.na(input_score)))
    stop("input_score must be a vector of continous values (with no NAs) where the vector names match the colnames of the expression matrix (required).\n")
  
  # Make sure the input_score has names or labels that are the same as the colnames of ES
  if(is.null(names(input_score)))
    stop("The input_score object must have names or labels to track the samples by. Please provide unique sample names or labels that matches the colnames of the expression matrix.\n")
  
  # Make sure the input_score has the same length as number of samples in ES
  if(length(input_score) != ncol(ES)){
    stop("The input_score must have the same length as the number of columns in ES.\n")
  }else{
    if(any(!names(input_score) %in% colnames(ES))){
      stop("The input_score object must have names or labels that match the colnames of the expression matrix.\n")
    }
    # match colnames of expression matrix with names of provided input_score
    ES <- ES[,names(input_score)]
  }
  
  # Check if the dataset has any all 0 or 1 features (these are to be removed since they are not informative)
  if(any(rowSums(exprs(ES)) == 0) || any(rowSums(exprs(ES)) == ncol(exprs(ES)))){
    warning("Provided dataset has features that are either all 0 or 1. These features will be removed from the computation.\n")
    ES <- ES[!(rowSums(exprs(ES)) == 0 | rowSums(exprs(ES)) == ncol(exprs(ES))),]
  }
  
  # Make sure matrix is not empty after removing uninformative features
  if(nrow(exprs(ES)) == 0){
    stop("After removing features that are either all 0 or 1. There are no more features remained for downsteam computation.\n")
  }
  
  # Check the method 
  if(length(method) == 1 & method %in% c("ks", "wilcox", "revealer", "custom")){
    
    # Compute row-wise directional KS scores for binary features in ES
    if(method == "ks"){
      verbose("Using Kolmogorov-Smirnov method for features scoring.\n")
      # Re-order the samples by input_score sorted from highest to lowest values
      ES <- ES[,names(sort(input_score, decreasing=T))]
    }
    
    # Compute row-wise Wilcox rank sum scores for binary features in ES 
    if(method == "wilcox"){
      verbose("Using Wilcoxon method for features scoring.\n")
      # Ranking the samples by input_score sorted from highest to lowest values
      ES <- ES[,names(sort(input_score, decreasing=T))]
    }
    
    # Compute mutually exclusive method for binary features in ES
    if(method == "revealer"){
      verbose("Using Revealer's Mutually Exclusive method for features scoring.\n")
    }
    
    # Compute row-wise directional scores using user's customized function for binary features in ES
    if(method == "custom"){
      verbose("Using a customized method for features scoring.\n")
    }
    
  } else {
    
    stop(paste0("Invalid method specified. The method can be ", paste0(c("ks", "wilcox", "revealer", "custom"), collapse="/"), "."))
    
  } 
  
  # Define scores based on specified metric of interest
  if(!metric %in% c('stat', 'pval'))
    stop("Please specify metric parameter as either 'stat' or 'pval' to use for candidate_search().\n")  
  
  # Select the appropriate method to compute scores based on skewdness of a given binary matrix
  mat <- exprs(ES)
  
  # Compute the row-wise scoring
  s <- switch(
    method,
    ks = ks_gene_score_mat(
      mat = mat,
      alternative = alternative, 
      weights = weights
    ),
    wilcox = wilcox_genescore_mat(
      mat = mat,
      alternative = alternative,
      ranks = NULL
    ),
    revealer = revealer_genescore_mat(
      mat = mat,                                   
      input_score = input_score,      
      seed_names = NULL,
      target_match = target_match,
      assoc_metric = "IC"
    ),
    custom = custom_genescore_mat(
      mat = mat,
      input_score = input_score,      
      custom_function = custom_function,
      custom_parameters = custom_parameters
    )
  )
  
  # Check if the returning result has one or two columns: score or p_value or both
  if(ncol(s) == 1){
    if(colnames(s) == "score" & metric == "pval"){
      warning(paste0("metric = 'pval' is provided but the ", method, " method only return score values. Thus, using 'stat' as metric to search for best features."))
      metric <- "stat"
    }
    if(colnames(s) == "p_value" & metric == "stat"){
      warning(paste0("metric provided is 'stat' but the ", method, "method only return p-values. Thus, using 'pval' as metric to search for best features."))
      metric <- "pval"
    }
    metric <- metric
  }
  
  ####### CACHE CHECKING #######
  if(!is.null(cache_path)){
    cat("Using provided cache root path: ", cache_path, "\n")
    setCacheRootPath(cache_path)
  } else{
    setCacheRootPath()
    cat("Setting cache root path as: ", getCacheRootPath(), "\n")
  }
  
  # We use the ES, top N (or search_start), score metric, scoring method and seed for random permutation as the key for each cached result  
  if(!is.null(top_N)){ # If N is defined here, we will use it as part of the key (topn_eval is called)
    key <- list(ES=ES, input_score=input_score, method=method, custom_function=custom_function, custom_parameters=custom_parameters, alternative=alternative, metric=metric, weights=weights, target_match=target_match, top_N=top_N, search_method=search_method, max_size=max_size, best_score_only=TRUE, seed=seed)
  }else{ # If N is not defined, we will use the search_start parameter as part of the key instead (candidate_search is called)
    key <- list(ES=ES, input_score=input_score, method=method, custom_function=custom_function, custom_parameters=custom_parameters, alternative=alternative, metric=metric, weights=weights, target_match=target_match, search_start=search_start, search_method=search_method, max_size=max_size, best_score_only=TRUE, seed=seed)
  }
  
  cat("Using the following as the key for saving/loading cached permutation values:\n")
  print(key)
  cat("\n\n")
  perm_best_scores <- loadCache(key)
  
  #Start the 'clock' to see how long the process takes
  ptm <- proc.time()
  
  ####### CACHE CHECKING #######
  
  # Check if, given the dataset and search-specific parameters, there is already a cached null distribution available 
  if (!is.null(perm_best_scores) & (length(perm_best_scores) >= n_perm)){
    
    cat("Found ", length(perm_best_scores), " cached permutation-based scores for the specified dataset and search parameters..\n")
    cat("Loading permutation scores from cache..\n")
    
  }  else{
    
    if (is.null(perm_best_scores)){
      cat("No permutation scores for the specified dataset and search parameters were found in cache path...\n")
    } else if (length(perm_best_scores) < n_perm) {
      cat("n_perm is set to ", n_perm, " but found ", length(perm_best_scores), " cached permutation-based scores for the specified dataset and search parameters...\n")
    }
    
    cat("\n\n\nBEGINNING PERMUTATION-BASED SIGNIFICANCE TESTING\n\n\n")
    
    ##############################################################################################
    # Sets up the parallel backend which will be utilized by Plyr.
    parallel = FALSE
    progress = "text"
    
    if(ncores > 1){
      registerDoParallel(cores = ncores)
      parallel = TRUE
      progress = "none"
      cat("Running tests in parallel..\n")
    } 
    
    cat("Using ", ncores, " core(s)..\n")
    
    # Generate matrix of permutated input_score  
    perm_labels_matrix <- generate_permutations(ord=input_score, n_perms=n_perm, seed=seed)
    
    #Set verbose to FALSE (override parameter specification) since we don't want to print any diagnostic statements
    options(verbose=FALSE)
    
    cat("Computing permutation-based scores for N = ", n_perm, "...\n\n")
    
    if(!is.null(top_N)){ # Run top N evaluation if N is specified
      perm_best_scores <- unlist(alply(perm_labels_matrix, 1, function(x){ perm_input_score=x; names(perm_input_score) <- names(input_score); topn_eval(ES=ES, input_score=perm_input_score, method=method, custom_function=custom_function, custom_parameters=custom_parameters, alternative=alternative, metric=metric, weights=weights, target_match=target_match, top_N=top_N, search_method=search_method, max_size=max_size, best_score_only=TRUE) },.parallel=parallel,.progress=progress))
    } else { # Run basic candidate search otherwise
      perm_best_scores <- unlist(alply(perm_labels_matrix, 1, function(x){ perm_input_score=x; names(perm_input_score) <- names(input_score); candidate_search(ES=ES, input_score=perm_input_score, method=method, custom_function=custom_function, custom_parameters=custom_parameters, alternative=alternative, metric=metric, weights=weights, target_match=target_match, search_start=search_start, search_method=search_method, max_size=max_size, best_score_only=TRUE) },.parallel=parallel,.progress=progress))  
    }
    
    #Save computed scores to cache 
    cat("Saving to cache ..\n")
    saveCache(perm_best_scores, key=key, comment="null_ks()")
    
  } # end caching else statement block
  
  registerDoParallel(cores = 1) #Return to using just a single core
  
  cat("FINISHED\n")
  cat("Time elapsed: ", round((proc.time()-ptm)[3]/60,2), " mins \n\n")
  ############################################################################################## 
  
  if(is.null(obs_best_score)){
    
    cat("Computing observed best score ..\n\n")
    
    if(!is.null(top_N)){
      
      obs_best_score <- topn_eval(ES = ES,
                                  input_score = input_score, 
                                  method = method,
                                  custom_function = custom_function,
                                  custom_parameters = custom_parameters,
                                  alternative = alternative,
                                  metric = metric,
                                  weights = weights,
                                  target_match = target_match,
                                  top_N = top_N,
                                  search_method = search_method,
                                  max_size = max_size,
                                  best_score_only = TRUE) %>% unlist()
      
    }else{
      
      obs_best_score <- candidate_search(ES = ES,
                                         input_score = input_score, 
                                         method = method,
                                         custom_function = custom_function,
                                         custom_parameters = custom_parameters,
                                         alternative = alternative,
                                         metric = metric,
                                         weights = weights,
                                         target_match = target_match,
                                         search_start = search_start,
                                         search_method = search_method,
                                         max_size = max_size,
                                         best_score_only = TRUE) %>% unlist()
      
    } 
    
  } else{
    
    cat("Using provided value of observed best score...\n\n")
    
  }
  
  cat("Observed score: ", obs_best_score, "\n\n")
  
  ########### PERMUTATION P-VALUE COMPUTATION ############
  cat("Number of permutation-based scores being considered: ", length(perm_best_scores), "\n")
  
  #Add a smoothening factor of 1 if smooth is specified
  #This is just to not return a p-value of 0
  c=0
  
  if(smooth)
    c=1
  
  if(metric == "pval"){
    
    #Use negative log transform of returned search score (either computed above, or passed to the null_ks function if previously computed)
    obs_best_score <- -(log(obs_best_score))
    
    # Use negative log transform on the permuted scores (either computed above or loaded from Cache)
    # NOTE: there is a very small chance some signed observed scores (p-values) are anti-correlated (meaning negative)
    # To avoid NaNs, we remove these. Keep in mind this is uncommon and will contribute very few permutations (n<10) if running N=1000
    perm_best_scores <- perm_best_scores[perm_best_scores > 0]
    perm_best_scores <- -(log(perm_best_scores))
    
  }
  
  perm_pval <- (sum(perm_best_scores > obs_best_score) + c)/(n_perm + c) 
  
  cat("Permutation p-value: ", perm_pval, "\n\n")
  
  ########### END PERMUTATION P-VALUE COMPUTATION ############
  
  if(plot == TRUE){
    
    plot_title <- paste("Emperical Null distribution (N = ", length(perm_best_scores), ")\n Permutation p-val <= ", round(perm_pval, 5), "\nBest observed score: ", round(obs_best_score, 5), sep="")
    
    if(!is.null(top_N)){
      plot_title <- paste(plot_title,"\n Top N: ", top_N, sep="")
    }else{
      plot_title <- paste(plot_title,"\n Seed: ", search_start, sep="")
    }
    
    #Here, let us plot the absolute values of the permutation p-values, for simplicity
    #You only consider absolute values when calculating the permutation p-values
    #Katia: adding ".data" to avoid a warning during check:
    # no visible binding for global variable 
    g <- ggplot(data = data.frame("x" = perm_best_scores), aes(x = .data$x)) +
      geom_histogram(fill = "black", color = "gray") +
      theme_classic() +
      theme(
        axis.line.x=element_line(color = "black"),
        axis.line.y=element_line(color = "black")
      )
    
    g <- g + geom_vline(xintercept = obs_best_score, linetype = "longdash", size = 1.5, colour = "red") +
      labs(
        title = plot_title,
        x = "Score",
        y = "Count"
      ) +
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0))
    
    print(g)
    
  }
  
  if(return_perm_pval){ return(list(top_N=top_N, search_start=search_start, perm_best_scores=perm_best_scores, perm_pval=perm_pval, obs_best_score=obs_best_score)) }
  
}


#' Random permutation matrix generator
#' 
#' Produces a random permutation rank matrix given a vector of values
#' @param ord vector to be permuted. This determines the number of columns in the permutation matrix
#' @param n_perms number of permutations to generate. This determines the number of rows in the permutation matrix
#' @param seed seed which can be set for reproducibility of 'random' results. Default is 123
#' @return A row matrix of permuted values (i.e. ranks) where each row is a single permutation result
#' @export
generate_permutations<-function(
  ord,        #These are the sample orderings to be permuted
  n_perms,    #Number of permutations to produce
  seed=123    #Seed which can be set for reproducibility of results
){
  
  m  <- length(ord)
  perm <- matrix(NA, nrow=n_perms, ncol=length(ord) )
  if ( !is.null(seed) ){
    set.seed(seed)
    verbose("Seed set: ", seed,"\n")
  }
  verbose("Generating ",n_perms," permuted sample ranks..\n")
  for ( i in (1:n_perms) ) {
    perm[i,] <- sample(ord,m)
  }
  verbose("Are all generated permutations unique?..\n")
  verbose(nrow(perm)==nrow(unique.matrix(perm)))
  
  if(nrow(perm)!=nrow(unique.matrix(perm)))
    stop("Not enough unique sample permutations for the permutation number specified.. Please provide a reasonabl nperm value ..\n")
  
  return(perm) 
  
}
