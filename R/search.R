verbose <- function(...){
  
  #Fetch verbose option set in the stepwise.search() function
  opt <- getOption("verbose",FALSE)
  if(!opt) return(invisible(NULL))
  msgs <- list(...)
  #msgs <- do.call(paste, c(msgs))
  message(msgs)
  
}

#' Pre-filter features
#' 
#' Pre-filter a dataset prior to running step-wise heuristic search in order to avoid testing features that are too prevalent or too sparse across samples in the dataset
#' @param ES an expression set object containing binary features used for step-wise search
#' @param max.cutoff a numeric value between 0 and 1 describing the absolute prevalence of a feature across all samples in the dataset above which the feature will be filtered out. Default is 0.6 (feature that occur in 60 percent or more of the samples will be removed)
#' @param min.cutoff a numeric value between 0 and 1 describing the absolute prevalence of a feature across all samples in the dataset below which the feature will be filtered out. Default is 0.03 (feature that occur in 3 percent or less of the samples will be removed)
#' @return An expression set object with only the filtered-in features given the filter thresholds specified
#' @examples
#' data(sim.ES)
#' 
#' # Filter out features having < 3 and > 60% prevalence across all samples (default)
#' sim.ES.filt1 <- prefilter_data(sim.ES)
#' 
#' # Change the min cut-off to 1% prevalence, instead of the default 3%
#' sim.ES.filt2 <- prefilter_data(sim.ES,min.cutoff=0.01)
#' 
#' # Change the max cut-off to 65% prevalence, instead of the default 60%
#' sim.ES.filt3 <- prefilter_data(sim.ES,max.cutoff=0.65) 
#' 
#' @export
#' @import Biobase
prefilter_data <- function(ES, 
                         max.cutoff=0.6,
                         min.cutoff=0.03){
  # Compute the frequency of feature occurence across all samples  (i.e. fraction of samples having the feature)
  frac <- round(rowSums(exprs(ES))/ncol(ES),2)
  
  cat("Pre-filtering features ..\n\n")
  cat("Removing features having < ",min.cutoff*100, "and > ",max.cutoff*100, " % occurence in sample set..\n")
  
  ES <- ES[ (frac >= min.cutoff) & (frac <= max.cutoff) , ]
  
  cat(nrow(ES)," features retained out of ",length(frac)," supplied features in dataset\n\n")
  return(ES)
  
}


#' Top 'N' evaluate
#' 
#' Generates and evaluates stepwise search results for the top 'N' starting indices, checking for overlapping resulting features from each case. This function is mainly used to evaluate search results over the top 'N' best starting features for a given dataset.
#' @param ranking integer vector specifying if and how samples should be re-ordered. Default is NULL. If NULL, we assume the ESet is already ordered  
#' @param ESet an ordered expression set object with the same sample ordering and features as processed by the stepwise.search() function when performing a step-wise heuristic search
#' @param N an integer specifying the number of features to start the search over, starting from the top 'N' features in each case. Default is 1
#' @param do.plot a logical indicating whether you want to plot the resulting evaluation matrix. Default is TRUE
#' @param best.score.only a logical indicating whether or not to only return the best meta-feature score over the top 'N' evaluation. Default is TRUE
#' @param ... additional parameters passed to the stepwise.search() function, which will be applied to each top 'N' run 
#' @return Default is a list of lists, where each list entry is one that is returned by the stepwise search run for a given starting index (See stepwise.search()). If best.score.only is set to TRUE, only the best score over the top N space is returned (useful for permutation-based testing)
#' 1's and 0's represent whether a feature in any given row is present in a meta-feature along with a starting feature in the corresponding column.
#' @export
#' @import Biobase gplots
topn.eval <- function(ranking=NULL,
                      ESet,
                      N=1,
                      do.plot=TRUE,
                      best.score.only=TRUE,
                      ...){
  
  if (N > nrow(ESet))
    stop("Please specify an N value that is less than the number of features in the ESet..\n")
  
  if (N > 10)
    warning("N value specified is greater than 10. This may result in longer search time..\n")

  verbose("Evaluating search over top features: ",1:N, "\n\n")
  
  #Performs stepwise search over top N indices
  topn.l <- sapply(1:N, function(x) stepwise.search(ranking = ranking,ES = ESet,cust_start = x,...),simplify = FALSE) 
  
  
  if(best.score.only==TRUE){
    scores.l <- lapply(topn.l, "[[", 2)
    
    # Working with scores for each top N run
    s <- unlist(scores.l)
    #Fetch the best score from the iterations
    # This ASSUMES you're using metric = "pval"
    # NEEDS UPDATING TO ACCOMODATE STATISTIC 
    best_score <- s[order(s)][1] #Based on the p-values, the lowest value will be the most sig 
    
    return(best_score)
  } #best.score.only
  
  
  if(do.plot){
    topn.plot(topN.list=topn.l)  
  } #do.plot
  
  return(topn.l) #Default is to return the top N stepwise search results as a list of lists
}

#' Top 'N' best
#' 
#' Takes the resulting list of meta-features returned from topn.eval() function and fetches the meta-feature with best (local) score
#' @param topn.list The nested list object that is returned by topn.eval() function when best.score.only is set to FALSE. This contains both the ESets as well as the scores for each resulting meta-feature in the top 'N'  search mode.
#' @return A list containing the (local) best meta-feature ESet, as well as its corresponding search score
#' @export
#' @import Biobase
topn.best <- function(topn.list){
  
  # Fetch the index housing the best ESet (this wil be the one with the best score)
  n <- which.min(sapply(topn.list,"[[",2))
  
  # Also store the score
  top.score <- min(sapply(topn.list,"[[",2))
  
  # Corresponding ESet object
  best.meta <- topn.list[[n]]$ESet
  
  return(list("ESet"=best.meta,
              "Score"=top.score))

}

backward_check <- function
(
  ESet,                    # an Expression Set object with the same sample ordering and features as processed by the stepwise.search() function 
  glob.f,                  # a vector containing the feature (row) names whose union gives the best score (so far) in the search. Feature names should match those of the provided expression set object
  glob.f.s,                # score corresponding to the union of the specified vector of features 
  m,                       # a character string specifying which metric to use for stepwise search criteria. One of either 'pval' or 'stat' may be used, corresponding to the  p-value or statistic. Uses value passed in the stepwise.search() function
  ...                      # additional parameters passed to the compute_score() function (method, alternative etc.)
  ){ 
  verbose("Performing backward search step..\n")
  verbose("Iterating over ",length(glob.f)," chosen features..\n")
  
  # Matrix of only global best features so far
  gmat <- exprs(ESet[glob.f,])
  rownames(gmat) <- glob.f 
  # Here, we make a list that should store the features and their corresponding meta-feature  score for each leave-one-out run
  f.names <- list()
  f.scores <- c()
  
  # We want to see if leaving anyone feature out improves the overall meta-feature  score
  # Here we only consider previous features in the meta-feature to remove (i.e. not the last one which was just added)
  for(n in 1:(length(glob.f)-1)){
    
    f.names[[n]] <- glob.f[-n]
    #Take leave-one-out union of features from matrix
    # This will result in a single vector to compute the scores on
    u <- ifelse(colSums(gmat[-n,])==0,0,1)
    
    # Compute  scores for this meta feature
    # Here we suprress warnings just to avoid messages warning-related single vector score computation (nrow(mat) < 2)
    u.s <- suppressWarnings(compute_score(mat=t(matrix(u)),...) )
    score <- ifelse(m %in% "pval",sign(u.s[1,])*u.s[2,], u.s[1,]) # Assuming bare=TRUE in compute_score call
    f.scores <- c(f.scores,score)
  } # end for loop
  
  if(m!="pval"){
    f.best.index <- which.max(f.scores) #This is the index within the meta matrix
  } else { #If signed pvalues
    f.best.index <- order(-sign(f.scores),f.scores)[1] #Top p-value ordered by sign and numerical value; #This is the index within the meta matrix
 }  
  f.best.score <- f.scores[f.best.index]
 
  # Check if any one of the computed scores has a better score than the entire meta-feature's score
  if(ifelse(m=="pval", (sign(f.best.score) > 0 & abs(f.best.score) < abs(glob.f.s)),f.best.score > glob.f.s)){
    verbose("Found improvement on removing existing feature..\n")
    verbose("New feature set: ",f.names[[f.best.index]],"\n")
    verbose("New global best score: ",f.best.score,"\n")
    
    #Return the set of features that gave a better score than the existing best score, and the score as well
    return(list(f.names[[f.best.index]],f.best.score))  
  } else{
    #Don't change anything. Return the existing best set of features and the corresponding  score
    return(list(glob.f,glob.f.s))
  } 
}

#' Step-wise search
#' 
#' Performs step-wise heuristic search using an ordered set of binary features to see whether there are features whose union is more skewed (enriched at the extremes) than either features alone. This is the main functionality of the CaDrA package.
#' @param ranking vector containing rankings for sample ordering. Default is NULL. If NULL, we assume the ESet is already ranked  
#' @param ES an expression set object of binary features. The first column of the featureData for the expression set must contain the names of the corresponding features, which are used in the search   
#' @param max.size an integer specifying the maximum size a meta-feature can extend do for a given search. Default is 7
#' @param method a character string specifying the method used to compute scores for features, must be one of "ks" or "wilcox"
#' @param metric a character string specifying which metric to use for stepwise search criteria. One of either 'pval' or 'stat' may be used, corresponding to the score p-value or statistic. Default is 'pval'
#' @param back_search a logical indicating whether or not to perform a forward-backward search (i.e. remove features along the search if it improves the best score). Default is TRUE. Uses function backward_check() 
#' @param cust_start an integer specifying a specific index within the expression set object of the feature to start the step-wise search with. Default is NULL. If NULL, then the search starts with the top ranked feature. If an integer is specified (N, where N < nrow(dataset)), then the search starts with the Nth best feature. If a string is specified, then the search starts with the feature with this name (must be a valid rowname in the dataset)
#' @param best.score.only a logical indicating whether or not the function should return only the score corresponding to the search results. Default is FALSE
#' @param alt a character string specifying the alternative hypothesis, must be one of "two.sided","greater" or "less". Default is "less" for left-skewed significance testing.
#' @param wts an integer vector of weights to use if performing weighted-KS testing. Default is NULL. Value passed to compute_score() function  
#' @param rnks an integer vector of sample rankings to use if performing Wilcoxon rank sum testing. Default is NULL. If NULL, then samples are assumed to be ordered by increasing ranking. Value passed to compute_score() function 
#' @param verb a logical indicating whether or not to print diagnostic messages. Default is FALSE 
#' @return If best.score.only is set to TRUE, this function returns a list object with the score corresponding to the union of the search meta-feature. If this is set to FALSE, a list containing both the ESet object pertaining to the returned meta-feature as well as the corresponding score  is returned. 
#' @export
#' @import Biobase 
stepwise.search <- function(ranking=NULL,  
                      ES, 
                      max.size=7,
                      metric="pval",
                      method=c("ks","wilcox"),
                      back_search=TRUE, 
                      cust_start=NULL,
                      best.score.only=FALSE, 
                      alt="less", 
                      wts=NULL,
                      rnks=NULL,
                      verb=FALSE
){
  
  # Setup verbose option definition
  options(verbose=verb)
  
  # We will rely on matrix row names, so make sure the input ES has rownames for features
  if(is.null(rownames(ES)))
    stop("ESet object provided does not have rownames/featureData to identify/track features by .. please provide unique feature names as rownames of the input ESet..\n")
  
  ###### SAMPLE PRE-RANKING #####
  #################################
  
  #Use ranking to re-order samples 
  if(!(is.null(ranking))){
    if(length(ranking)!=ncol(ES))
      stop("Ranking variable has to be of the same length as the number of samples\n\n")
    verbose("Using provided ordering to re-rank samples..\n")
    ES <- ES[,ranking]
  }
  
  if(length(method)==1 & method %in% c("ks","wilcox"))
  {
    #Compute row-wise directional KS scores for existing (raw/starting) binary features in ESet
    if(method=="ks"){
      verbose("Using KS method for feature scoring..\n")
      if(!(is.null(wts)))
        verbose("Using weighted method for KS testing using provided weights..\n")
    }
    #Compute row-wise Wilcox rank sum scores for existing (raw/starting) binary features in ESet 
    if(method=="wilcox"){
      verbose("Using Wilcoxon method for feature scoring..\n")
      if(!(is.null(rnks)))
        verbose("Using provided ranks for Wilcoxon rank sum testing..\n")
    } 
  } else {
    stop("Invalid method specification to compute scores.. please specify either 'ks' or 'wilcox'..\n")
  }
  
  ###### INVALID FEATURE REMOVAL #####
  #####################################
  # Check if the dataset has any all 0 or 1 features (these are to be removed since they are not informative)
  
  if(any(rowSums(exprs(ES))==0) | any(rowSums(exprs(ES))==ncol(ES))){
    warning("Provided dataset has features that are either all 0 or 1.. removing these features prior to beginning search\n\n")
    ES <- ES[!(rowSums(exprs(ES))==0 | rowSums(exprs(ES))==ncol(ES)),]
  }
  
  # Compute initial scores per feature given the dataset
  
  s <- compute_score(mat=exprs(ES),method=method,alt=alt,weight=wts,ranks=rnks)
  
  s.stat <- s[1,]
  s.pval <- s[2,]
  
  # Define scores based on specified metric of interest
  
  if(!metric %in% c('stat','pval'))
    stop("Please specify metric parameter as either 'stat' or 'pval' to use for search..\n")
  
  score <- ifelse(rep(metric,nrow(ES)) %in% "pval",sign(s.stat)*s.pval, s.stat)
  verbose("Using ",metric," as measure of improvement measure ..\n\n")
  
  
  ###### FEATURE PRE-RANKING #####
  ##################################
  
  
  # Re-order ESet in decreasing order of user-defined score (s stat or p-val)
  # This comes in handy when doing the top-N evaluation of the top N 'best' features
  
  score.rank <- if (metric!="pval") order(score) else order(-sign(score),score)
  verbose("Ranking ESet features by metric..\n")
  
  ES <- ES[score.rank,]
  score <- score[score.rank]
  
  if(is.null(cust_start)){ 
    verbose("Starting with feature having best ranking ..\n")
    best.s.index <- 1  
  } else {
    if(is.numeric(cust_start)){ 
      # User-specified feature index (has to be an integer from 1:nrow(ES))
      verbose("Starting with specified sorted feature index ..\n")
      
      if(cust_start > nrow(ES)) # Index out of range
        stop("Invalid starting index specified.. Please specify a valid starting index within the range of the existing ESet..\n")
      
      best.s.index <- cust_start 
    }
    
    if(is.character(cust_start)){
      # User-specified feature name (has to be a character from rownames(1:nrow(ES)))
      verbose("Starting with specified feature name ..\n")
      
      if(!(cust_start %in% rownames(ES))) #provided feature name not in rownames
        stop("Provided starting feature does not exist among ESet's rownames .. please check stepwise.search cust_start parameter options for more details ..\n\n")
      
      best.s.index <- which(rownames(ES)==cust_start)  
    } # end if is.character 
  } # end else (!is.null)
  ###### INITIALIZE VARIABLES ###########
  #######################################
  
  # Print the featureData for this starting point feature so that we are aware
  # Here we assume the ESet's fData is included as rownames
  #start.feature <- as.character(fData(ES)[best.s.index,1])
  start.feature <- rownames(ES)[best.s.index]
  best.feature <- start.feature
  best.s <- score[best.s.index]
  
  
  verbose("Feature: ",start.feature,"\n")
  verbose("Score: ",best.s,"\n")
  
  #Fetch the vector corresponding to best score
  #Set this as the initial 'meta-feature'
  best.meta <- as.numeric(exprs(ES)[best.s.index,])
  
  
  #counter variable for number of iterations
  i=0
  
  #Variable to store best score attained over all iterations
  #initialize this to the starting best score
  global.best.s <- best.s
  
  # Vector of features in the (growing) obtained metafeature. Begin with just the starting feature
  global.best.s.features <- c()
  
  ###### BEGIN ITERATIONS ###############
  #######################################
  
  verbose("\n\nBeginning stepwise search..\n\n")
  
  while ((ifelse(metric=="pval", (sign(best.s) > 0 & (abs(best.s) < abs(global.best.s))),best.s > global.best.s) | i == 0) & (length(global.best.s.features) < max.size)){
    verbose("\n\n")
    verbose("Iteration number ",(i+1)," ..\n")
    verbose("Global best score: ",global.best.s,"\n")
    verbose("Previous score: ",best.s,"\n")
    
    
    
    # Update scores and feature set since since entry into the loop means there is an improvement (iteration > 0)
    global.best.s <- best.s
    global.best.s.features <- c(global.best.s.features,best.feature)
    
    verbose("Current feature set: ",global.best.s.features,"\n")
    
    if(i!=0){
      verbose("Found feature that improves  score!\n")
      # Update the new best meta feature (from meta mat)
      best.meta <- meta.mat[hit.best.s.index,]  
      #Add that index to the group of indices to be excluded for subsequent checks
      #Here we go off the rownames in the original matrix to find which index to exclude from the ESet in subsequent iterations
      best.s.index <- c(best.s.index,which(rownames(ES)==best.feature))
    } 
    
    # Perform a backward check on the list of existing features and update global scores/feature lists accordingly  
    if(length(global.best.s.features) > 3 & back_search==TRUE){
      
      backward_search.results <- backward_check(ESet=ES,
                                                glob.f=global.best.s.features, #Global feature set so far
                                                glob.f.s=global.best.s, # score corresponding to this global feature set
                                                m=metric,
                                                method=method,  # passed to compute_score() function
                                                alt=alt,        # passed to compute_score() function
                                                wts=wts,        # passed to compute_score() function
                                                rnks=rnks)      # passed to compute_score() function
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
    if(any(rowSums(meta.mat)==ncol(meta.mat))){
      warning("Features with all 1's generated upon taking matrix union .. removing such features before progressing..\n")
      meta.mat <- meta.mat[rowSums(meta.mat) != ncol(meta.mat),]
    }
    
    
    
    #With the newly formed 'meta-feature' matrix, compute directional  scores and choose the feature that gives the best  score
    #Compute row-wise directional  scores for existing (raw/starting) binary features in ESet
    s <- compute_score(mat=meta.mat,method=method,alt=alt,weight=wts,ranks=rnks)
    s.stat <- s[1,]
    s.pval <- s[2,]
    
    
    # Take signed pval or stat depending on user-defined metric
    # This will be the same length as nrow(meta.mat)
    scores <- ifelse(rep(metric,nrow(meta.mat)) %in% "pval",sign(s.stat)*s.pval,s.stat)
    
    
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
    verbose("Feature that produced best score in combination with previous meta-feature: ",best.feature,"\n")
    verbose("Score: ",best.s,"\n")
    
    # If no improvement (exiting loop)
    if(ifelse(metric=="pval", sign(best.s) < 0 | (abs(best.s) >= abs(global.best.s)),best.s <= global.best.s)){
      verbose("No further improvement in score has been found..\n")
    }
    
    #Increment counter
    i=i+1
    
  } #########End of while loop
  
  verbose("\n\n")
  verbose("\n\nFinished!\n\n")
  verbose("Number of iterations covered: ",i,"\n")
  verbose("Best  score attained over iterations: ",global.best.s,"\n")
  if(length(global.best.s.features)==1){
    warning("No meta-feature that improves the enrichment was found ..\n") 
  }
  
  verbose("Features returned in ESet: ",global.best.s.features,"\n")
  verbose("\n\n")
  
  if(best.score.only==FALSE){
    
    #We don't just want the combination (meta-feature) at the end. We want all the features that make up the meta-feature
    #This can be obtained using the list of indices that were progressively excluded (if at all) in the step-wise procedure
    #If returning only those features that led to the best global score
    ES.best <- ES[global.best.s.features,]
    #Here, give the returned ESet an annotation based on the starting feature that gave these best results
    annotation(ES.best) <- start.feature
    colnames(ES.best) <- colnames(ES)
    
    #Make a list contaning two elements. 
    #The first will be the ESet with the features that gave the best meta-feature
    #The second will be the score corresponding to the meta-feature (named by the starting feature that led to that score)
    #Assign the name of the best meta-feature score to be the starting feature that gave that score
    names(global.best.s) <- start.feature 
    
    return(list("ESet"= ES.best,"Score"= global.best.s))
    
  } else{
    #Just return the score. Here we put this in the list just to support permutation-based apply functionality
    return(list(global.best.s)) 
  }
} # end stepwise.search function

#' Random permutation matrix generator
#' 
#' Produces a random permutation rank matrix given a vector of values
#' @param ord vector to be permuted. This determines the number of columns in the permutation matrix
#' @param n_perms number of permutations to generate. This determines the number of rows in the permutation matrix
#' @param seed seed which can be set for reproducibility of 'random' results. Default is 123
#' @return A row matrix of permuted values (i.e. ranks) where each row is a single permutation result
#' @export
generate_permutations <- function(ord, #These are the sample orderings to be permuted
                                n_perms, #Number of permutations to produce
                                seed=123 #Seed which can be set for reproducibility of results
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

#' Permutation-based step-wise searching
#' 
#' Performs permutation-based significance testing of step-wise search results.
#' @param ranking a vector containing rankings for sample ordering. Default is NULL. If NULL, we assume the ESet is already ranked. This is used only to compute the observed stepwise search score, if not provided  
#' @param ES an ordered expression set object with the same sample ordering and features as processed by the stepwise.search() function when performing a step-wise heuristic search
#' @param nperm an integer specifying the number of permutations to perform. Default = 1000
#' @param plot logical indicating whether or not to plot the emperical null distribution with the observed score and permutation p-value
#' @param obs.best.score a numeric value corresponding to the observed (best) stepwise search score to use for permutation based p-value computation. Default is NULL. If set to NULL, we compute the observed score given the ranking variable and ESet
#' @param cache.path full path to permutation-based (null) score distributions cache files. If a permutation for a given dataset (and dependent search variables such as 'N') exist, we recycle values instead of re-computing them to save time. Default is NULL. If NULL, cache path is set to the default ~/.Rcache for future cache loading.
#' @param smooth logical indicating whether or not to smoothen the p-value calculation to avoid p-value of 0. Default is TRUE
#' @param return.perm.pval logical indicating whether or not to return the permutation-based p-value computed by the function. Default is FALSE 
#' @param seed seed set for permutation. Default = 123
#' @param ncores number of cores to use, if using parallelization for permutation testing. Default = 1
#' @param ... additional parameters passed to the stepwise.seach() function, which will be applied to each permutation run (called within topn.eval())
#' @return If return.perm.pval is set to TRUE, will return the permutation p-value
#' @export
#' @import Biobase R.cache doParallel ggplot2 plyr
null.search <- function(ranking=NULL,
                  ES,
                  nperm=1000,
                  plot=TRUE,
                  obs.best.score=NULL,
                  cache.path=NULL,
                  smooth=TRUE,
                  return.perm.pval=FALSE,
                  seed=123,
                  ncores=1,
                  ...){
  
  #Here we fetch the additional arguments passed to stepwise.search() so that it can be accessed
  null.args <- list(...)
  topN = null.args$N # This will either be null (if not specified) or user-defined for the top-N evaluation
  m <- null.args$metric
  met <- null.args$method
  m.size <- null.args$max.size # This is the max size specification parameter for the search
  
  # If there isn't any overriding of the default max meta-feature size
  if(is.null(m.size))
    m.size <- 7 # Assign the default (this is only for the key assignment)
  
  # If there isn't any overriding of the default metric to use for the search
  if(is.null(m))
    m <- "pval" # Assign the default (this is only for the key assignment)
  
  # If there isn't any specification for the method to use (CaDrA needs one of "ks" or "wilcox" specified)
  if(is.null(met)) # Throw an error
    stop('CaDrA requires a statistical method specification. Please see ?stepwise.search for options when specifying this "method" parameter, which will be applied to all null searches\n\n')
  
  ####### CACHE CHECKING #######
  
  if(!is.null(cache.path)){
    cat("Using provided cache root path: ",cache.path,"\n")
    setCacheRootPath(cache.path)
    
  } else{
    setCacheRootPath()
    cat("Setting cache root path as: ",getCacheRootPath(),"\n")
  }
  
  # We use the ESet, top N (or cust_start), score metric, scoring method and seed for random permutation as the key for each cached result  
  if(!is.null(topN)) # If N is defined here, we will use it as part of the key (topn.eval is called)
   key <- list(ESet=ES,topN=topN,method=met,metric=m,max_size=m.size,seed=seed)
  else # If N is not defined, we will use the cust_start parameter as part of the key instead (stepwise.search is called)
   key <- list(ESet=ES,cust_start=null.args$cust_start,method=met,metric=m,max_size=m.size,seed=seed)
  
  cat("Using the following as the key for saving/loading cached permutation values:\n")
  print(key)
  cat("\n\n")
  perm.best.scores <- loadCache(key)
  
  #Start the 'clock' to see how long the process takes
  ptm<-proc.time()
  
  ####### CACHE CHECKING #######
  
  # Check if, given the dataset and search-specific parameters, there is already a cached null distribution available 
  if (!is.null(perm.best.scores) & (length(perm.best.scores) >= nperm)){
    
    cat("Found ",length(perm.best.scores), " cached permutation-based scores for the specified dataset and search parameters..\n")
    cat("Loading permutation scores from cache..\n")
    
  }  else{
    
    if (is.null(perm.best.scores)){
      
      cat("No permutation scores for the specified dataset and search parameters were found in cache path ..\n")
      
    } else if (length(perm.best.scores) < nperm) {
      
      cat("nperm is set to ",nperm," but found ",length(perm.best.scores), " cached permutation-based scores for the specified dataset and search parameters..\n")
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
    perm_labels_matrix<-generate_permutations(ord=seq(1,ncol(ES)),
                                              n_perms=nperm,
                                              seed=seed)
    
    
    #Set verbose to FALSE (override parameter specification) since we don't want to print any diagnostic statements
    options(.null_ks.verb=FALSE)
    
    
    cat("Computing permutation-based scores for N = ",nperm," ..\n\n")
    if(!is.null(topN)){ # Run top N evaluation if N is specified
      perm.best.scores<-unlist(alply(perm_labels_matrix,1,topn.eval,ESet = ES,best.score.only=TRUE,verb=FALSE,...,.parallel=parallel,.progress = progress))
    } else { # Run basic stepwise search otherwise
      perm.best.scores<-unlist(alply(perm_labels_matrix,1,stepwise.search,ES = ES,best.score.only=TRUE,verb=FALSE,...,.parallel=parallel,.progress = progress))  
    }
    #Save computed scores to cache 
    cat("Saving to cache ..\n")
    saveCache(perm.best.scores,key=key,comment="null_ks()")
    
  } # end caching else statement block
  
  registerDoParallel(cores = 1) #Return to using just a single core
  
  cat("FINISHED\n")
  cat("Time elapsed: ",round((proc.time()-ptm)[3]/60,2)," mins \n\n")
  ############################################################################################## 
  
  
  if(is.null(obs.best.score)){
    cat("Computing observed best score ..\n\n")
    if(!is.null(topN)){
    obs.best.score <- topn.eval(ranking=ranking, #This is passed to the null_ks function
                              ESet = ES,
                              best.score.only = TRUE,
                              verb=FALSE,
                              ...) 
    } else {
    obs.best.score <- stepwise.search(ranking=ranking, #This is passed to the null_ks function
                                  ES = ES,
                                  best.score.only = TRUE,
                                  verb=FALSE,
                                  ...) 
      
    } 
  } else{
    cat("Using provided value of observed best score ..\n\n")
  }
  
  cat("Observed score: ",unlist(obs.best.score),"\n\n")

  ########### PERMUTATION P-VALUE COMPUTATION ############
  cat("Number of permutation-based scores being considered: ",length(perm.best.scores), "\n")
  #Add a smoothening factor of 1 if specified
  #This is just to not return a p-value of 0
  c=0
  if(smooth)
    c=1
  
  if(m=="pval"){
    
    #Use negative log transform of returned search score (either computed above, or passed to the null_ks function if previously computed)
    obs.best.score <- -(log(unlist(obs.best.score)))
    
    # Use negative log transform on the permuted scores (either computed above or loaded from Cache)
    # NOTE: there is a very small chance some signed observed scores (p-values) are anti-correlated (meaning negative)
    # To avoid NaNs, we remove these. Keep in mind this is uncommon and will contribute very few permutations (n<10) if running N=1000
    perm.best.scores <- perm.best.scores[perm.best.scores > 0]
    perm.best.scores <- -(log(perm.best.scores))
    
  }
  
  perm.pval <- (sum(perm.best.scores > obs.best.score) + c)/(nperm + c) 
  
  
  cat("Permutation p-value: ",perm.pval,"\n\n")
  
  
  ########### END PERMUTATION P-VALUE COMPUTATION ############
  
  if(plot==TRUE){
    plot.title <- paste("Emperical Null distribution (N = ",length(perm.best.scores),")\n Permutation p-val <= ",round(perm.pval,5),"\nBest observed score: ",round(obs.best.score,5),sep="")
    if(!is.null(topN))
      plot.title <- paste(plot.title,"\n Top N: ",topN,sep="")
    else
      plot.title <- paste(plot.title,"\n Seed: ",null.args$cust_start,sep="")
    
    #Here, let us plot the absolute values of the permutation p-values, for simplicity
    #You only consider absolute values when calculating the permutation p-values
    #Katia: adding ".data" to avoid a warning during check:
    # no visible binding for global variable 
    g <- ggplot(data=data.frame("x"=perm.best.scores),aes(x=.data$x))+
      geom_histogram(fill="black",color="gray")+
      theme_classic()+theme(axis.line.x=element_line(color="black"),
                            axis.line.y=element_line(color="black"))
    
    g <- g+geom_vline(xintercept=obs.best.score,linetype="longdash",size=1.5,colour="red")+
      labs(title=plot.title,
           x="Score",
           y="Count")+
      theme(plot.title = element_text(hjust = 0.5))+
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0))
    
    print(g)
  }
  
  if(return.perm.pval)
    return(perm.pval)
}
