#' Permutation-based step-wise searching
#' 
#' Performs permutation-based significance testing of step-wise search results.
#'
#' @param ES an expression set object of binary features (required). It must be a BioBase expressionSet object. The rownames or featureData of the expression set must contain the names of the corresponding features which are used in the search.   
#' @param input_score a vector of continuous values (required). 
#' @param method a character string specifying the method used to compute scores for features, must be one of "ks" or "wilcox" or "mi" (mutually exclusive method from REVEALER) or "custom" (a personal customization method). If input_score contains ranked scores, then 'ks' method is used by default. Otherwise, 'mi" is the default method
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided", "greater" or "less". Default is "less" for left-skewed significance testing.
#' @param metric a character string specifying which metric to use for candidate search. One of either 'pval' or 'stat' may be used, corresponding to the score p-value or statistic. Default is 'pval'
#' @param search_method a character string specifying which method to perform or filter out the best candidates. Default is 'forward'.
#' @param search_start an integer specifying a specific index within the expression set object of the feature to start the step-wise search with. Default is NULL. If NULL, then the search starts with the top ranked feature. If an integer is specified (N, where N < nrow(dataset)), then the search starts with the Nth best feature. If a string is specified, then the search starts with the feature with this name (must be a valid rowname in the dataset)
#' @param max_size an integer specifying the maximum size a meta-feature can extend do for a given search. Default is 7.
#' @param top_N an integer specifying the number of features to start the search over, starting from the top 'N' features in each case. Default is 1.
#' @param n_perm an integer specifying the number of permutations to perform. Default = 1000.
#' @param plot logical indicating whether or not to plot the emperical null distribution with the observed score and permutation p-value
#' @param obs_best_score a numeric value corresponding to the observed (best) stepwise search score to use for permutation based p-value computation. Default is NULL. If set to NULL, we compute the observed score given the ranking variable and ESet
#' @param cache_path full path to permutation-based (null) score distributions cache files. If a permutation for a given dataset (and dependent search variables such as 'N') exist, we recycle values instead of re-computing them to save time. Default is NULL. If NULL, cache path is set to the default ~/.Rcache for future cache loading.
#' @param smooth logical indicating whether or not to smoothen the p-value calculation to avoid p-value of 0. Default is TRUE
#' @param return_perm_pval logical indicating whether or not to return the permutation-based p-value computed by the function. Default is FALSE 
#' @param seed seed set for permutation. Default = 123.
#' @param ncores number of cores to use, if using parallelization for permutation testing. Default = 1..
#' 
#' @return If \code{return_perm_pval} is set to \code{TRUE} will return a permutation p-value.
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
  alternative = "less",
  metric = "pval",
  search_method = "both",
  search_start = NULL,
  max_size = 7,  
  top_N = NULL,
  n_perm = 1000,
  plot = TRUE,
  obs_best_score = NULL,
  cache_path = NULL,
  smooth = TRUE,
  return_perm_pval = FALSE,
  seed = 123,
  ncores = 1
){
  
  # Check if the ES is provided
  if(length(ES) == 0 || class(ES)[1] != "ExpressionSet") 
    stop("'ES' must be an  ExpressionSet class argument (required).")
  
  # Check input_score is provided
  if(length(input_score) == 0 || !is.numeric(input_score))
    stop("input_score must be a vector of continous values where the vector names matched colnames of ExpressionSet (required).\n")
  
  # Make sure the input ES has rownames for features tracking
  if(is.null(rownames(ES)))
    stop("The ES object does not have rownames or featureData to track the features by. Please provide unique features or rownames for the expression matrix.\n")
  
  # Make sure the input_score has names as the colnames of ES
  if(is.null(names(input_score)))
    stop("The input_score object must have names or labels to track the samples by. Please provide unique sample names or labels that matches the colnames of the expression matrix.\n")
  
  # Make sure the input_score has the same length as number of samples in ES
  if(length(input_score) != ncol(ES)){
    stop("The input_score must have the same length as the number of columns in ES.\n")
  }else{
    if(any(names(input_score) != colnames(ES))){
      stop("The input_score object must have names or labels that matches colnames of the expression matrix.\n")
    }
  }
  
  # Check if the dataset has only binary 0 or 1 values 
  if(!all(exprs(ES) %in% c(0,1))){
    stop("The expression matrix (ES) must contain only binary values (0 or 1).\n")
  }
  
  # Check if the dataset has any all 0 or 1 features (these are to be removed since they are not informative)
  if(any(rowSums(exprs(ES)) == 0) || any(rowSums(exprs(ES)) == ncol(exprs(ES)))){
    warning("Provided dataset has features that are either all 0 or 1. These features will be removed from the computation.\n")
    ES <- ES[!(rowSums(exprs(ES)) == 0 | rowSums(exprs(ES)) == ncol(exprs(ES))),]
  }
  
  # Check the method 
  if(length(method)==1 & method %in% c("ks", "wilcox", "revealer", "custom")){
    
    # Compute row-wise directional KS scores for binary features in ES
    if(method == "ks"){
      verbose("Using Kolmogorov-Smirnov method for features scoring.\n")
      ES <- ES[,names(sort(input_score, decreasing=T))]
    }
    
    # Compute row-wise Wilcox rank sum scores for binary features in ES 
    if(method == "wilcox"){
      verbose("Using Wilcoxon method for features scoring.\n")
    }
    
    # Compute mutually exclusive method for binary features in ES 
    if(method == "revealer"){
      verbose("Using Revealer's Mutually Exclusive method for features scoring.\n")
    }
    
    # Other future methods can be implemented here and add its verbose message here
    if(method == "custom"){
      verbose("Using a customized method for features scoring.\n")
    }
    
  } else {
    
    stop(paste0("Invalid method specified. The method can be ", paste0(c("ks", "wilcox", "revealer", "custom"), collapse="/"), "."))
    
  } 
  
  ####### CACHE CHECKING #######
  if(!is.null(cache_path)){
    cat("Using provided cache root path: ",cache_path,"\n")
    setCacheRootPath(cache_path)
  } else{
    setCacheRootPath()
    cat("Setting cache root path as: ",getCacheRootPath(),"\n")
  }
  
  # We use the ESet, top N (or search_start), score metric, scoring method and seed for random permutation as the key for each cached result  
  if(!is.null(top_N)) # If N is defined here, we will use it as part of the key (topn_eval is called)
    key <- list(ESet=ES,input_score=input_score,method=method,alternative=alternative,metric=metric,search_method=search_method,max_size=max_size,N=top_N,seed=seed)
  else # If N is not defined, we will use the search_start parameter as part of the key instead (candidate_search is called)
    key <- list(ESet=ES,input_score=input_score,method=method,alternative=alternative,metric=metric,search_method=search_method,max_size=max_size,search_start=search_start,seed=seed)
  
  cat("Using the following as the key for saving/loading cached permutation values:\n")
  print(key)
  cat("\n\n")
  perm_best_scores <- loadCache(key)
  
  #Start the 'clock' to see how long the process takes
  ptm <- proc.time()
  
  ####### CACHE CHECKING #######
  
  # Check if, given the dataset and search-specific parameters, there is already a cached null distribution available 
  if (!is.null(perm_best_scores) & (length(perm_best_scores) >= n_perm)){
    
    cat("Found ",length(perm_best_scores), " cached permutation-based scores for the specified dataset and search parameters..\n")
    cat("Loading permutation scores from cache..\n")
    
  }  else{
    
    if (is.null(perm_best_scores)){
      
      cat("No permutation scores for the specified dataset and search parameters were found in cache path ..\n")
      
    } else if (length(perm_best_scores) < n_perm) {
      
      cat("n_perm is set to ",n_perm," but found ",length(perm_best_scores), " cached permutation-based scores for the specified dataset and search parameters..\n")
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
    
    cat("Using ",ncores," core(s)..\n")
    
    # Generate matrix of permutated labels  
    perm_labels_matrix <- generate_permutations(ord=seq(1,ncol(ES)), n_perms=n_perm, seed=seed)
    
    #Set verbose to FALSE (override parameter specification) since we don't want to print any diagnostic statements
    options(verbose=FALSE)
    
    cat("Computing permutation-based scores for N = ",n_perm," ..\n\n")
    
    if(!is.null(top_N)){ # Run top N evaluation if N is specified
      perm_best_scores <- unlist(alply(perm_labels_matrix,1,topn_eval,ESet=ES,input_score=input_score,method=method,alternative=alternative,metric=metric,search_method=search_method,max_size=max_size,best_score_only=TRUE,N=top_N,verbose=FALSE,.parallel=parallel,.progress=progress))
    } else { # Run basic stepwise search otherwise
      perm_best_scores <- unlist(alply(perm_labels_matrix,1,candidate_search,ES=ES,input_score=input_score,method=method,alternative=alternative,metric=metric,search_start=search_start,search_method=search_method,max_size=max_size,best_score_only=TRUE,verbose=FALSE,.parallel=parallel,.progress=progress))  
    }
    
    #Save computed scores to cache 
    cat("Saving to cache ..\n")
    saveCache(perm_best_scores, key=key, comment="null_ks()")
    
  } # end caching else statement block
  
  registerDoParallel(cores = 1) #Return to using just a single core
  
  cat("FINISHED\n")
  cat("Time elapsed: ", round((proc.time()-ptm)[3]/60,2)," mins \n\n")
  ############################################################################################## 
  
  if(is.null(obs_best_score)){
    cat("Computing observed best score ..\n\n")
    
    if(!is.null(top_N)){
      obs_best_score <- topn_eval(ESet = ES,
                                  input_score = input_score, 
                                  method = method,
                                  alternative = alternative,
                                  metric = metric,
                                  search_method = search_method,
                                  max_size = max_size,
                                  best_score_only = TRUE,
                                  N = top_N,
                                  verbose = FALSE) 
    } else {
      obs_best_score <- candidate_search(ES = ES,
                                         input_score = input_score, 
                                         method = method,
                                         alternative = alternative,
                                         metric = metric,
                                         search_start = search_start,
                                         search_method = search_method,
                                         max_size = max_size,
                                         best_score_only = TRUE,
                                         verbose = FALSE) 
      
    } 
  } else{
    cat("Using provided value of observed best score ..\n\n")
  }
  
  cat("Observed score: ",unlist(obs_best_score),"\n\n")
  
  ########### PERMUTATION P-VALUE COMPUTATION ############
  cat("Number of permutation-based scores being considered: ",length(perm_best_scores), "\n")
  #Add a smoothening factor of 1 if specified
  #This is just to not return a p-value of 0
  c=0
  
  if(smooth)
    c=1
  
  if(metric == "pval"){
    
    #Use negative log transform of returned search score (either computed above, or passed to the null_ks function if previously computed)
    obs_best_score <- -(log(unlist(obs_best_score)))
    
    # Use negative log transform on the permuted scores (either computed above or loaded from Cache)
    # NOTE: there is a very small chance some signed observed scores (p-values) are anti-correlated (meaning negative)
    # To avoid NaNs, we remove these. Keep in mind this is uncommon and will contribute very few permutations (n<10) if running N=1000
    perm_best_scores <- perm_best_scores[perm_best_scores > 0]
    perm_best_scores <- -(log(perm_best_scores))
    
  }
  
  perm_pval <- (sum(perm_best_scores > obs_best_score) + c)/(n_perm + c) 
  
  
  cat("Permutation p-value: ",perm_pval,"\n\n")
  
  
  ########### END PERMUTATION P-VALUE COMPUTATION ############
  
  if(plot==TRUE){
    
    plot_title <- paste("Emperical Null distribution (N = ",length(perm_best_scores),")\n Permutation p-val <= ",round(perm_pval,5),"\nBest observed score: ",round(obs_best_score,5),sep="")
    
    if(!is.null(top_N)){
      plot_title <- paste(plot_title,"\n Top N: ",top_N,sep="")
    }else{
      plot_title <- paste(plot_title,"\n Seed: ",search_start,sep="")
    }
    
    #Here, let us plot the absolute values of the permutation p-values, for simplicity
    #You only consider absolute values when calculating the permutation p-values
    #Katia: adding ".data" to avoid a warning during check:
    # no visible binding for global variable 
    g <- ggplot(data=data.frame("x"=perm_best_scores),aes(x=.data$x))+
      geom_histogram(fill="black",color="gray")+
      theme_classic()+
      theme(
        axis.line.x=element_line(color="black"),
        axis.line.y=element_line(color="black")
      )
    
    g <- g + geom_vline(xintercept=obs_best_score,linetype="longdash",size=1.5,colour="red") +
      labs(
        title=plot_title,
        x="Score",
        y="Count"
      )+
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0))
    
    print(g)
    
  }
  
  if(return_perm_pval)
    return(perm_pval)
  
}
