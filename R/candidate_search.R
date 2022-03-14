
#' Candidate Search
#' 
#' Performs heuristic search on a set of binary features to determine whether there are features whose union is more skewed (enriched at the extremes) than either features alone. This is the main functionality of the CaDrA package.
#' @param ES an expression set of binary features (required). It must be a BioBase expressionSet object. The rownames of the expression set must contain unique features which are used in the search.   
#' @param input_score a vector of continuous values for a target profile (required). The input_score must have names or labels that matches the colnames of the expression matrix.
#' @param method a character string specifies a method to compute the score for each feature (\code{"ks"} or \code{"wilcox"} or \code{"revealer"} (conditional mutual information from REVEALER) or \code{"custom"} (a customized method)). Default is \code{ks}.
#' @param custom_function if method is \code{"custom"}, specifies the customized function here. Default is \code{NULL}.
#' @param custom_parameters if method is \code{"custom"}, specifies a list of arguments to be passed to the custom_function(). Default is \code{NULL}.
#' @param alternative a character string specifies an alternative hypothesis testing (\code{"two.sided"} or \code{"greater"} or \code{"less"}). Default is \code{less} for left-skewed significance testing.
#' @param metric a character string specifies a metric to search for best features. \code{"pval"} or \code{"stat"} may be used, corresponding to p-value or score statistic. Default is \code{pval}. Note: \code{Revealer} method only return score statistics values.
#' @param weights a vector of weights use to perform a weighted-KS testing. Default is \code{NULL}.   
#' @param ranks a vector of sample rankings use to perform Wilcoxon rank sum testing. Default is \code{NULL}. If NULL, samples are assumed to be ordered by increasing rankings.
#' @param target_match a direction of target matching (\code{"negative"} or \code{"positive"}) from REVEALER. Use \code{"positive"} to match the higher values of the target, \code{"negative"} to match the lower values. Default is \code{positive}. 
#' @param search_start an integer specifies an index within the expression set object of which feature to start the candidate search with. Default is \code{NULL}. If NULL, then the search starts with the top ranked feature. If an integer is specified (N, where N < nrow(ES)), the search starts with the Nth best feature. If a string is specified, the search starts with the feature with this name (must be a valid rowname in ES)
#' @param search_method a character string specifies a method to filter out the best candidates (\code{"forward"} or \code{"both"}). Default is \code{both} (backward and forward).
#' @param max_size an integer specifies a maximum size that a meta-feature can extend to do for a given search. Default is \code{7}.
#' @param best_score_only a logical value indicates whether or not the function should return only the score corresponding to the search results. Default is \code{FALSE}.
#' @param verbose a logical value indicates whether or not to print the diagnostic messages. Default is \code{FALSE}. 
#'
#' @return If \code{best_score_only} is set to \code{TRUE}, this function returns a list object with the score corresponding to the union of the search meta-feature. If \code{best_score_only} is set to \code{FALSE}, a list containing both the ES object pertaining to the returned meta-feature as well as the corresponding score is returned. 
#' @examples
#' # Load R library
#' library(Biobase)
#'
#' # Load pre-computed expression set
#' data(sim.ES)
#' 
#' # Provide a vector of continuous scores for a target profile with names to each score value 
#' input_score = rnorm(n = ncol(sim.ES))
#' names(input_score) <- colnames(sim.ES)
#' 
#' # Define additional parameters and run the function
#' candidate_search_result <- candidate_search(
#'   ES = sim.ES, input_score = input_score, method = "ks", 
#'   alternative = "less", weights = NULL, metric = "pval", 
#'   search_start = NULL, search_method = "both", 
#'   max_size = 7, best_score_only = FALSE
#' )
#' 
#' @export
candidate_search <- function(
  ES, 
  input_score, 
  method = "ks",
  custom_function = NULL,
  custom_parameters = NULL,
  alternative = "less", 
  metric = "stat",
  weights = NULL,
  ranks = NULL,  
  target_match = "positive",          
  search_start = NULL,
  search_method = "both",
  max_size = 7,
  best_score_only = FALSE,
  verbose = FALSE
){
  
  # Set up verbose option
  options(verbose=FALSE)
  
  # Check if the ES is provided
  if(length(ES) == 0 || class(ES)[1] != "ExpressionSet") 
    stop("'ES' must be an ExpressionSet class argument (required).")
  
  # Check if the dataset has only binary 0 or 1 values 
  if(!all(exprs(ES) %in% c(0,1))){
    stop("The expression matrix (ES) must contain only binary values with no NAs.\n")
  }
  
  # Make sure the input ES has rownames for features tracking
  if(is.null(rownames(ES)))
    stop("The ES object does not have rownames or featureData to track the features by. Please provide unique features or rownames for the expression matrix.\n")
    
  # Check input_score is provided
  if(length(input_score) == 0 || any(!is.numeric(input_score)) || any(is.na(input_score)))
    stop("input_score must be a vector of continous values (with no NAs) where the vector names match the colnames of the expression matrix (required).\n")
  
  # Make sure the input_score has same names as the colnames of ES
  if(is.null(names(input_score)))
    stop("The input_score object must have names or labels to track the samples by. Please provide unique sample names or labels that matches the colnames of the expression matrix.\n")
  
  # Make sure the input_score has the same length as number of samples in ES
  if(length(input_score) != ncol(ES)){
    stop("The input_score must have the same length as the number of columns in ES.\n")
  }else{
    if(any(!names(input_score) %in% colnames(ES))){
      stop("The input_score object must have names or labels that match the colnames of the expression matrix.\n")
    }
    # match colnames of expression matrix with names of provided input_score values
    if(nrow(ES) == 1){
      ES <- t(ES[,names(input_score)]) 
    }else{
      ES <- ES[,names(input_score)]
    }
  }
  
  # Check if the dataset has any all 0 or 1 features (these are to be removed since they are not informative)
  if(any(rowSums(exprs(ES)) == 0) || any(rowSums(exprs(ES)) == ncol(exprs(ES)))){
    warning("Provided dataset has features that are either all 0 or 1. These features will be removed from the computation.\n")
    ES <- ES[!(rowSums(exprs(ES)) == 0 | rowSums(exprs(ES)) == ncol(exprs(ES))),]
  }
  
  if(nrow(exprs(ES)) == 0){
    stop("After removing features that are either all 0 or 1. There are no more features remained for downsteam computation.\n")
  }
  
  # Check the method 
  if(length(method) == 1 & method %in% c("ks", "wilcox", "revealer", "custom")){
    
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
  
  # Define scores based on specified metric of interest
  if(!metric %in% c('stat', 'pval'))
    stop("Please specify metric parameter as either 'stat' or 'pval' to use for candidate_search().\n")  
  
  # Select the appropriate method to compute scores based on skewdness of a given binary matrix
  mat <- exprs(ES)
    
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
      ranks = ranks
    ),
    revealer = revealer_genescore_mat(
      mat = mat,                                   
      target = input_score,      
      seed_names = NULL,
      target_match = target_match,
      assoc_metric = "IC"
    ),
    custom = custom_genescore_mat(
      mat = mat,
      custom_function = custom_function,
      custom_parameters = custom_parameters
    )
  )
  
  # Check if the returning result has one or two columns: score or p_value or both
  if(ncol(s) == 1){
    if(colnames(s) == "score" & metric == "pval"){
      warning("metric = 'pval' is provided but the method function only return score values. Thus, using 'stat' as metric to search for best features.")
      metric <- "stat"
    }
    if(colnames(s) == "p_value" & metric == "stat"){
      warning("metric provided is 'stat' but the method function only return p-values. Thus, using 'pval' as metric to search for best features.")
      metric <- "pval"
    }
    metric <- metric
  }
  
  # Score returned by the given method
  s.stat <- if("score" %in% colnames(s)){ s[,"score"] }
  s.pval <- if("p_value" %in% colnames(s)){ s[,"p_value"] }
  
  # compute the scores according to the provided metric 
  score <- ifelse(rep(metric,nrow(ES)) %in% "pval", sign(s.stat)*s.pval, s.stat)
  
  verbose("Using ", metric, " as measure of improvement measure...\n\n")
  
  ###### FEATURE PRE-RANKING #####
  ##################################
  
  # Re-order ES in decreasing order of user-defined score (stat or pval)
  # This comes in handy when doing the top-N evaluation of the top N 'best' features
  score.rank <- if(metric != "pval") order(score) else order(-sign(score), score)
  
  verbose("Ranking ES features by metric..\n")
  
  ES <- ES[score.rank,]
  
  score <- score[score.rank]
  
  if(is.null(search_start)){ 
    
    verbose("Starting with feature having best ranking...\n")
    best.s.index <- 1  
    
  } else {
    
    if(is.numeric(search_start)){ 
      # User-specified feature index (has to be an integer from 1:nrow(ES))
      verbose("Starting with specified sorted feature index ..\n")
      
      if(search_start > nrow(ES)) # Index out of range
        stop("Invalid starting index specified... Please specify a valid starting index within the range of the existing ES...\n")
      
      best.s.index <- search_start 
    }
    
    if(is.character(search_start)){
      # User-specified feature name (has to be a character from rownames(1:nrow(ES)))
      verbose("Starting with specified feature name ..\n")
      
      if(!(search_start %in% rownames(ES))) #provided feature name not in rownames
        stop("Provided starting feature does not exist among ES's rownames.\n\n")
      
      best.s.index <- which(rownames(ES) == search_start)  
    } # end if is.character 
    
  } # end else (!is.null)
  
  ###### INITIALIZE VARIABLES ###########
  #######################################
  
  # Print the featureData for this starting point feature so that we are aware
  # Here we assume the ES's fData is included as rownames
  #start.feature <- as.character(fData(ES)[best.s.index,1])
  start.feature <- rownames(ES)[best.s.index]
  best.feature <- start.feature
  best.s <- score[best.s.index]
  
  verbose("Feature: ", start.feature, "\n")
  verbose("Score: ", best.s, "\n")
  
  #Fetch the vector corresponding to best score
  #Set this as the initial 'meta-feature'
  best.meta <- as.numeric(exprs(ES)[best.s.index,])
  
  #counter variable for number of iterations
  i=0
  
  #Variable to store best score attained over all iterations
  #initialize this to the starting best score
  global.best.s <- best.s
  
  # Vector of features in the (growing) obtained meta-feature. Begin with just the starting feature
  global.best.s.features <- c()
  
  ###### BEGIN ITERATIONS ###############
  #######################################
  
  verbose("\n\nBeginning candidate search...\n\n")
  
  if(search_method == "both"){ 
    back_search = TRUE 
  }else if(search_method == "forward"){ 
    back_search = FALSE
  }else {
    stop("Incorrect search_method parameter provided. The search_method can be either 'both' or 'forward'.")
  }
  
  while((ifelse(metric=="pval", (sign(best.s) > 0 & (abs(best.s) < abs(global.best.s))), best.s > global.best.s) | i == 0) & (length(global.best.s.features) < max_size)){
    
    verbose("\n\n")
    verbose("Iteration number ", (i+1), " ..\n")
    verbose("Global best score: ", global.best.s, "\n")
    verbose("Previous score: ", best.s, "\n")
    
    # Update scores and feature set since since entry into the loop means there is an improvement (iteration > 0)
    global.best.s <- best.s
    global.best.s.features <- c(global.best.s.features, best.feature)
    
    verbose("Current feature set: ", global.best.s.features, "\n")
    
    if(i != 0){
      
      verbose("Found feature that improves score!\n")
      # Update the new best meta feature (from meta mat)
      best.meta <- meta.mat[hit.best.s.index,]  
      #Add that index to the group of indices to be excluded for subsequent checks
      #Here we go off the rownames in the original matrix to find which index to exclude from the ES in subsequent iterations
      best.s.index <- c(best.s.index,which(rownames(ES)==best.feature))
      
    } 
    
    # Perform a backward check on the list of existing features and update global scores/feature lists accordingly  
    if(length(global.best.s.features) > 3 & back_search == TRUE){
      
      backward_search.results <- forward_backward_check(ES = ES,
                                                        input_score = input_score,
                                                        glob.f = global.best.s.features, # Global feature set so far
                                                        glob.f.s = global.best.s,        # score corresponding to this global feature set
                                                        method = method, 
                                                        custom_function = custom_function,
                                                        custom_parameters = custom_parameters,
                                                        alternative = alternative,
                                                        metric = metric,
                                                        weights = weights,       
                                                        ranks = ranks,
                                                        target_match = target_match)    
      # Update globlal features, scores 
      global.best.s.features <- backward_search.results[[1]]
      global.best.s <- backward_search.results[[2]]
      # Update best.meta based on feature set
      best.meta <- as.numeric(ifelse(colSums(exprs(ES)[global.best.s.features,])==0,0,1))
      
    }
    
    #Take the OR function between that feature and all other features, to see which gives the best  score
    #Keep in mind, the number of rows in meta.mat keeps reducing by one each time we find a hit that improves the  score
    verbose("Forming meta-feature matrix with all other features in dataset..\n")
    # Here "*1" is used to convert the boolean back to integer 1's and 0's
    # Notice we remove anything in best.s.index from the original matrix first, to form the meta matrix.
    meta.mat <- sweep(exprs(ES)[-best.s.index,],2,best.meta,`|`)*1
    
    # Check if there are any features that are all 1's generated on taking the union
    # We cannot compute statistics for such features and they thus need to be filtered out
    if(any(rowSums(meta.mat) == ncol(meta.mat))){
      warning("Features with all 1's generated upon taking matrix union .. removing such features before progressing..\n")
      meta.mat <- meta.mat[rowSums(meta.mat) != ncol(meta.mat),]
    }
    
    #With the newly formed 'meta-feature' matrix, compute directional  scores and choose the feature that gives the best score
    #Compute row-wise directional  scores for existing (raw/starting) binary features in ES
    s <- switch(
      method,
      ks = ks_gene_score_mat(
        mat = meta.mat,
        alternative = alternative, 
        weights = weights
      ),
      wilcox = wilcox_genescore_mat(
        mat = meta.mat,
        alternative = alternative,
        ranks = ranks
      ),
      revealer = revealer_genescore_mat(
        mat = meta.mat,                                   
        target = input_score,      
        seed_names = NULL,
        target_match = target_match,
        assoc_metric = "IC"
      ),
      custom = custom_genescore_mat(
        mat = meta.mat,
        custom_function = custom_function,
        custom_parameters = custom_parameters
      )
    ) 

    # Check if the returning result has one or two columns: score or p_value or both
    if(ncol(s) == 1){
      if(colnames(s) == "score" & metric == "pval"){
        warning("metric = 'pval' is provided but the method function only return score values. Thus, using 'stat' as metric to search for best features.")
        metric <- "stat"
      }
      if(colnames(s) == "p_value" & metric == "stat"){
        warning("metric provided is 'stat' but the method function only return p-values. Thus, using 'pval' as metric to search for best features.")
        metric <- "pval"
      }
      metric <- metric
    }
    
    # Score returned by the given method
    s.stat <- if("score" %in% colnames(s)){ s[,"score"] }
    s.pval <- if("p_value" %in% colnames(s)){ s[,"p_value"] }
    
    # Take signed pval or stat depending on user-defined metric
    # This will be the same length as nrow(meta.mat)
    scores <- ifelse(rep(metric,nrow(meta.mat)) %in% "pval", sign(s.stat)*s.pval, s.stat)
    
    #Find index of feature that gives lowest s score when combined with chosen starting feature
    if(metric!="pval"){
      hit.best.s.index <- which.max(scores) #This is the index within the meta matrix
    } else { #If signed pvalues
      hit.best.s.index <- order(-sign(scores),scores)[1] #Top p-value ordered by sign and numerical value; #This is the index within the meta matrix
    }
    
    best.s <- scores[hit.best.s.index] #This is the best score from the meta matrix
    
    # Find which feature produced that score, in combination with meta feature used
    # We go from index to rowname space here in the meta matrix
    # We can do this because rownames are preserved between the original and meta features on using sweep()
    best.feature <- rownames(meta.mat)[hit.best.s.index]
    verbose("Feature that produced best score in combination with previous meta-feature: ", best.feature, "\n")
    verbose("Score: ",best.s,"\n")
    
    # If no improvement (exiting loop)
    if(ifelse(metric == "pval", sign(best.s) < 0 | (abs(best.s) >= abs(global.best.s)), best.s <= global.best.s)){
      verbose("No further improvement in score has been found...\n")
    }
    
    #Increment counter
    i = i+1
    
  } #########End of while loop
  
  verbose("\n\n")
  verbose("\n\nFinished!\n\n")
  verbose("Number of iterations covered: ",i,"\n")
  verbose("Best  score attained over iterations: ",global.best.s,"\n")
  if(length(global.best.s.features)==1){
    warning("No meta-feature that improves the enrichment was found ..\n") 
  }
  
  verbose("Features returned in ES: ",global.best.s.features,"\n")
  verbose("\n\n")
  
  if(best_score_only==FALSE){
    
    #We don't just want the combination (meta-feature) at the end. We want all the features that make up the meta-feature
    #This can be obtained using the list of indices that were progressively excluded (if at all) in the step-wise procedure
    #If returning only those features that led to the best global score
    ES.best <- ES[global.best.s.features,]
    #Here, give the returned ES an annotation based on the starting feature that gave these best results
    annotation(ES.best) <- start.feature
    colnames(ES.best) <- colnames(ES)
    
    #Make a list contaning two elements. 
    #The first will be the ES with the features that gave the best meta-feature
    #The second will be the score corresponding to the meta-feature (named by the starting feature that led to that score)
    #Assign the name of the best meta-feature score to be the starting feature that gave that score
    names(global.best.s) <- start.feature 
    
    return(list("ES"= ES.best,"Score"= global.best.s))
    
  } else{
    #Just return the score. Here we put this in the list just to support permutation-based apply functionality
    return(list(global.best.s)) 
  }
  
}  

# Performance backward selection 
forward_backward_check <- function
(
  ES,                      # an Expression Set object with the same sample ordering and features as processed by the stepwise.search() function 
  input_score,
  glob.f,                  # a vector containing the feature (row) names whose union gives the best score (so far) in the search. Feature names should match those of the provided expression set object
  glob.f.s,                # score corresponding to the union of the specified vector of features 
  method,
  custom_function,
  custom_parameters,
  alternative,
  metric,                  # a character string specifying which metric to use for stepwise search criteria. One of either 'pval' or 'stat' may be used, corresponding to the  p-value or statistic. Uses value passed in the candidate_search() function
  weights,
  ranks,
  target_match
){ 
  
  verbose("Performing backward search...\n")
  verbose("Iterating over ", length(glob.f), " chosen features...\n")
  
  # Matrix of only global best features so far
  gmat <- exprs(ES[glob.f,])
  rownames(gmat) <- glob.f 
  
  # Here, we make a list that should store the features and their corresponding meta-feature  score for each leave-one-out run
  f.names <- list()
  f.scores <- c()
  
  # We want to see if leaving anyone feature out improves the overall meta-feature  score
  # Here we only consider previous features in the meta-feature to remove (i.e. not the last one which was just added)
  for(n in 1:(length(glob.f)-1)){
    #n=1;
    f.names[[n]] <- glob.f[-n]
    
    #Take leave-one-out union of features from matrix
    # This will result in a single vector to compute the scores on
    u <- ifelse(colSums(gmat[-n,]) == 0, 0, 1)
    
    # Compute scores for this meta feature
    # Here we suprress warnings just to avoid messages warning-related single vector score computation (nrow(mat) < 2)
    u.mat <- matrix(t(matrix(u)), nrow=1, byrow=T, dimnames=list(c("sum"), colnames(ES)))

    s <- suppressWarnings(
      switch(
        method,
        ks = ks_gene_score_mat(
          mat = u.mat,
          alternative = alternative, 
          weights = weights
        ),
        wilcox = wilcox_genescore_mat(
          mat = u.mat,
          alternative = alternative,
          ranks = ranks
        ),
        revealer = revealer_genescore_mat(
          mat = u.mat,                                   
          target = input_score,      
          seed_names = NULL,
          target_match = target_match,
          assoc_metric = "IC"
        ),
        custom = custom_genescore_mat(
          mat = u.mat,
          custom_function = custom_function,
          custom_parameters = custom_parameters
        )
      )
    )
    
    # Check if the returning result has one or two columns: score or p_value or both
    if(ncol(s) == 1){
      if(colnames(s) == "score" & metric == "pval"){
        warning("metric = 'pval' is provided but the method function only return score values. Thus, using 'stat' as metric to search for best features.")
        metric <- "stat"
      }
      if(colnames(s) == "p_value" & metric == "stat"){
        warning("metric provided is 'stat' but the method function only return p-values. Thus, using 'pval' as metric to search for best features.")
        metric <- "pval"
      }
      metric <- metric
    }
    
    # Score returned by the given method
    s.stat <- if("score" %in% colnames(s)){ s[,"score"] }
    s.pval <- if("p_value" %in% colnames(s)){ s[,"p_value"] }
    
    score <- ifelse(metric %in% "pval", sign(s.stat)*s.pval, s.stat)
    
    f.scores <- c(f.scores, score)
    
  } # end for loop
  
  if(metric != "pval"){
    f.best.index <- which.max(f.scores) #This is the index within the meta matrix
  } else { #If signed pvalues
    f.best.index <- order(-sign(f.scores),f.scores)[1] #Top p-value ordered by sign and numerical value; #This is the index within the meta matrix
  }  
  
  f.best.score <- f.scores[f.best.index]
  
  # Check if any one of the computed scores has a better score than the entire meta-feature's score
  if(ifelse(metric == "pval", (sign(f.best.score) > 0 & abs(f.best.score) < abs(glob.f.s)), f.best.score > glob.f.s)){
    
    verbose("Found improvement on removing existing feature..\n")
    verbose("New feature set: ",f.names[[f.best.index]],"\n")
    verbose("New global best score: ",f.best.score,"\n")
    
    #Return the set of features that gave a better score than the existing best score, and the score as well
    return(list(f.names[[f.best.index]],f.best.score))  
    
  } else{
    
    #Don't change anything. Return the existing best set of features and the corresponding  score
    return(list(glob.f, glob.f.s))
    
  } 
  
}


