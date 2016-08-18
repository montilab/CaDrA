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
#' @export
prefilter_data<-function(ES, 
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
    best_score <- s[order(s)][1] #Based on the p-values, the lowest value will be the most sig
    
    return(best_score)
  } #best.score.only
  
  
  if(do.plot){
    topn.plot(topN.list=topn.l)  
  } #do.plot
  
  return(topn.l) #Default is to return the top N stepwise search results as a list of lists
}

backward_check <- function
(
  ESet,                    # an Expression Set object with the same sample ordering and features as processed by the stepwise.search() function 
  glob.f,                  # a vector containing the feature names whose union gives the best score (so far) in the search. Feature names should match those of the provided expression set object
  glob.f.s,                # score corresponding to the union of the specified vector of features 
  m,                       # a character string specifying which metric to use for stepwise search criteria. One of either 'pval' or 'stat' may be used, corresponding to the  p-value or statistic. Uses value passed in the stepwise.search() function
  ...                      # additional parameters passed to the compute_score() function (method, alternative etc.)
  ){ 
  verbose("Performing backward search step..\n")
  verbose("Iterating over ",length(glob.f)," chosen features..\n")
  #print(glob.f)
  
  # Matrix of only global best features so far
  gmat <- exprs(ESet[glob.f,])
  rownames(gmat) <- glob.f 
  # Here, we make a list that should store the features and their corresponding meta-feature  score for each leave-one-out run
  f.names <- list()
  f.scores <- c()
  
  # We want to see if leaving anyone feature out improves the overall meta-feature  score
  for(n in 1:(length(glob.f)-1)){
    
    f.names[[n]] <- glob.f[-n]
    #Take leave-one-out union of features from matrix
    u <- ifelse(colSums(gmat[-n,])==0,0,1)
    
    #Compute  scores for this meta feature
    u.s <- compute_score(mat=t(matrix(u)),...) 
    score <- ifelse(m %in% "pval",sign(u.s[1,])*u.s[2,], u.s[1,]) # Assuming bare=TRUE in compute_score call
    f.scores <- c(f.scores,score)
  }
  if(m!="pval"){
    f.best.index <- which.max(f.scores) #This is the index within the meta matrix
  } else { #If signed pvalues
    f.best.index <- order(-sign(f.scores),f.scores)[1] #Top p-value ordered by sign and numerical value; #This is the index within the meta matrix
 }  
  f.best.score <- f.scores[f.best.index]
 
  # Check if any one of the computed scores has a better score than the entire meta-feature's score
  if(ifelse(m=="pval", (sign(f.best.score)>0 & abs(f.best.score) < abs(glob.f.s)),f.best.score > glob.f.s)){
    verbose("Found improvement on removing existing feature..\n")
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
#' @param metric a character string specifying which metric to use for stepwise search criteria. One of either 'pval' or 'stat' may be used, corresponding to the score p-value or statistic
#' @param back_search a logical indicating whether or not to perform a forward-backward search (i.e. remove features along the search if it improves the best score). Default is TRUE. Uses function backward_check() 
#' @param cust_start an integer specifying a specific index within the expression set object of the feature to start the step-wise search with. Default is NULL 
#' @param best_score_only a logical indicating whether or not the function should return only the score corresponding to the search results. Default is FALSE
#' @param alt a character string specifying the alternative hypothesis, must be one of "two.sided","greater" or "less". Default is "two.sided"
#' @param wts an integer vector of weights to use if performing weighted-KS testing. Default is NULL. Value passed to compute_score() function  
#' @param rnks an integer vector of sample rankings to use if performing Wilcoxon rank sum testing. Default is NULL. If NULL, then samples are assumed to be ordered by increasing ranking. Value passed to compute_score() function 
#' @param verb a logical indicating whether or not to print diagnostic messages. Default is FALSE 
#' @return If best_score_only is set to TRUE, this function returns a list object with the score corresponding to the union of the search meta-feature. If this is set to FALSE, an expression set object containing the features whose union gave the best score is returned. 
#' @export
stepwise.search<-function(ranking=NULL,  
                      ES, 
                      max.size=7,
                      metric,
                      method=c("ks","wilcox"),
                      back_search=TRUE, 
                      cust_start=NULL,
                      best_score_only=FALSE, 
                      alt="two.sided", 
                      wts=NULL,
                      rnks=NULL,
                      verb=FALSE
){
  
  #Setup verbose option definition
  options(verbose=verb)
  
  
  ###### SAMPLE PRE-RANKING #####
  
  #Use ranking to re-order samples 
  if(!(is.null(ranking))){
    if(length(ranking)!=ncol(ES))
      stop("Ranking variable has to be of the same length as the number of samples\n\n")
    verbose("Using provided ordering to re-rank samples..\n")
    ES<-ES[,ranking]
  }
  
  if(length(method)==1 & method %in% c("ks","wilcox"))
  {
  #Compute row-wise directional KS scores for existing (raw/starting) binary features in ESet
  if(method=="ks"){
    verbose("Using KS method for feature scoring..\n")
    if(!(is.null(wts)))
      verbose("Using weighted method for KS testing using provided weights..\n")
  } 
  if(method=="wilcox"){
      verbose("Using Wilcoxon method for feature scoring..\n")
      if(!(is.null(rnks)))
        verbose("Using provided ranks for Wilcoxon rank sum testing..\n")
    } 
  } else {
      stop("Invalid method specification to compute scores.. please specify either 'ks' or 'wilcox'..\n")
    }
  
  s <- compute_score(mat=exprs(ES),method=method,alt=alt,weight=wts,ranks=rnks)
  s.stat <- s[1,]
  s.pval <- s[2,]
  
  #Define scores based on specified metric of interest
  
  if(!metric %in% c('stat','pval'))
    stop("Please specify metric parameter as either 'stat' or 'pval' to use for search..\n")
  score <- ifelse(rep(metric,nrow(ES)) %in% "pval",sign(s.stat)*s.pval, s.stat)
  
  
  verbose("Using ",metric," as measure of improvement measure ..\n\n")
  
  
  ###### FEATURE PRE-RANKING #####
  
  #ADDED
  #Order ESet in increasing/decreasing order of user-defined score (s stat or p-val)
  #This comes in handy when doing the top-N evaluation of the top N 'best' features
  
  score.rank <- if (metric!="pval") order(score) else order(-sign(score),score)
  verbose("Ranking ESet features by metric..\n")
  ES <- ES[score.rank,]
  score <- score[score.rank]
  
  
  ###### FEATURE PRE-RANKING #####
  
  
  #Here, we will assume ASSIGN scores have samples ranked in decreasing order
  #This means an higher/positive stat is associated with a higher ASSIGN score 
  
  #Let us start with the first (top ranked) feature
  #Fetch index of feature having best score. We start here
  if(is.null(cust_start)){
    
    verbose("Starting with feature having best ranking ..\n")
    #MODIFIED
    #best.s.index<-ifelse(metric!="pval",which.max(score),order(-sign(score),score)[1]) #This assumes that samples are ordered in decreasing order of ASSIGN score
    best.s.index <- 1 
    
    
  } else {
    
    verbose("Starting with specified sample feature ..\n")
    if(!(cust_start <= nrow(ES)))
      stop("Invalid starting index specified.. Please specify a valid starting index within the range of the existing ESet..\n")
    best.s.index <- cust_start
    
  } 
  
  best.s <- score[best.s.index]
  
  #Print the featureData for this starting point feature so that we are aware
  start_feature <- as.character(fData(ES)[best.s.index,1])
  verbose("Feature: ",start_feature,"\n")
  verbose("Score: ",best.s,"\n")
  
  #Fetch the vector corresponding to best score
  #Set this as the initial 'meta-feature'
  best.meta <- as.numeric(exprs(ES)[best.s.index,])
  
  #This is just so that it enters the while loop for the first iteration
  new.best.s <- best.s
  
  ###### INITIALIZE VARIABLES ###########
  #######################################
  
  #counter variable for number of iterations
  i=0
  #Parameter for number of continuous 'mistakes' (no improvement in score)
  b=0
  
  #Variable to store best score attained over all iterations
  #initialize this to the starting best score
  global.best.s <- best.s
  #Variable to store indices to make meta-feature with best score attained over all iterations
  #initialize this to the starting feature's index
  global.best.s.index <- best.s.index
  global.best.s.features <- c(as.character(fData(ES)[best.s.index,1]))
  
  
  ###### BEGIN ITERATIONS ###############
  #######################################
  
  verbose("\n\nBeginning stepwise search..\n\n")
  
  #This condition defines the overall search criteria
  while ((ifelse(metric=="pval", (sign(new.best.s)>0 & abs(new.best.s) < abs(best.s)),new.best.s > best.s) | b < 2) & (length(global.best.s.features) < max.size)){
    verbose("\n\n")
    verbose("Iteration number ",i+1," ..\n")
    
    
    if(i!=0){
      #Not the first iteration	
      
      verbose("New best s: ",new.best.s,"\n")
      verbose("Best s: ",best.s,"\n")
      
      if(ifelse(metric=="pval", (sign(new.best.s)>0 & abs(new.best.s) < abs(best.s)),new.best.s > best.s))
      {
        verbose("Found feature that improves  score!\n")
        #Now that we have an improvement, let us reset the 'continuous mistake counter' to 0
        b=0
      }
      new.best.meta <- meta.mat[hit.best.s.index,]
      
      #Add that index to the group of indices to be excluded for subsequent checks
      #Here we go off the rownames to find which index to exclude from the ESet
      best.s.index <- c(best.s.index,which(rownames(ES)==best.feature))
      
      #If performing a forward-backward search, we need to check if adding this last feature
      #works better when leaving any of the existing best features out
      #This is only useful if you have at least 4 features
      if(length(global.best.s.features) > 3 & back_search==T){
        
        backward_search.results <- backward_check(ESet=ES,
                                                glob.f=global.best.s.features, #Global feature set so far
                                                glob.f.s=global.best.s, # score corresponding to this global feature set
                                                m=metric,
                                                method=method,  # passed to compute_score() function
                                                alt=alt,        # passed to compute_score() function
                                                wts=wts,        # passed to compute_score() function
                                                rnks=rnks)      # passed to compute_score() function
        global.best.s.features <- backward_search.results[[1]]
        global.best.s <- backward_search.results[[2]]
      }
      
      #Reset current minimum values and feature to new minimum values and feature
      best.meta <- new.best.meta
      best.s <- new.best.s
      verbose(" score: ", best.s,"\n")
    }
    
    
    #Take the OR function between that feature and all other features, to see which gives the best  score
    #Keep in mind, the number of rows in meta.mat keeps reducing by one each time we find a hit that improves the  score
    verbose("Forming meta-feature matrix with all other features in dataset..\n")
    meta.mat <- sweep(exprs(ES)[-best.s.index,],2,best.meta,`|`)*1
    
    #verbose("Number of rows in meta-feature matrix: ",nrow(meta.mat),"\n")
    #print(rownames(meta.mat))
    
    
    #####<Still need to implement this functionality>######
    
    
    #Let us see if any of the newly formed meta features are the same as the feature itself
    #Leave out these metafeatures as they add no new information
    
    #####<Still need to implement this functionality>######
    
    
    #With the newly formed 'meta-feature' matrix, compute directional  scores and choose the feature that gives the best  score
    #Compute row-wise directional  scores for existing (raw/starting) binary features in ESet
    s <- compute_score(mat=meta.mat,method=method,alt=alt,weight=wts,ranks=rnks)
    s.stat <- s[1,]
    s.pval <- s[2,]
    
    
    scores <- ifelse(rep(metric,nrow(meta.mat)) %in% "pval",sign(s.stat)*s.pval,s.stat)
    
    
    
    #Find index of feature that gives lowest s score when combined with chosen starting feature
    if(metric!="pval")
      hit.best.s.index <- which.max(scores) #This is the index within the meta matrix
    else #If signed pvalues
      hit.best.s.index <- order(-sign(scores),scores)[1] #Top p-value ordered by sign and numerical value; #This is the index within the meta matrix
    new.best.s <- scores[hit.best.s.index] #This is within the meta matrix
    
    #Diagnostic
    #verbose("Best score from meta matrix: ",new.best.s,"\n")
    
    
    #Find which feature produced that score, in combination with meta feature used
    best.feature <- rownames(meta.mat)[hit.best.s.index]
    verbose("Feature that produced best score in combination with previous meta-feature: ",best.feature,"\n")
    #Diagnostics
    #verbose("Meta mat index: ",hit.best.s.index,"\n")
    #verbose("Estimated ES index for feature: ", which(rownames(ES)==best.feature),"\n")
    #verbose("Actual index of feature in ES: ",which(rownames(ES)==best.feature),"\n")
    
    
    #Checking whether it's a new best score
    if(ifelse(metric=="pval", (sign(new.best.s)>0 & abs(new.best.s) < abs(global.best.s)),new.best.s > global.best.s))
      #if(new.best.s<global.best.s)
      #Update best score
    { global.best.s <- new.best.s
      global.best.s.index <- c(global.best.s.index,hit.best.s.index+(nrow(ES)-nrow(meta.mat))) #Here we add the difference in rows because it is relative 
      global.best.s.features <- c(global.best.s.features,best.feature)
    }
    
    verbose("Global best score so far: ",global.best.s,"\n")
    #to the entire ES which will have x more rows than the meta matrix 
    #Let us set a parameter that can allow for a few mistakes
    #If the score isn't improving for the metafeature, add a counter (constrained in the while loop entry)
    if(ifelse(metric=="pval", sign(new.best.s) < 0 | (abs(new.best.s) >= abs(best.s)),new.best.s <= best.s)){
      verbose("Allowing for one non-minimizing step..\n")
      b=b+1
    }
    #Increment counter
    i=i+1
  } #########End of while loop
  
  
  verbose("\n\n")
  verbose("Number of iterations covered: ",i,"\n")
  verbose("Number of continuous mistakes made: ",b,"\n")
  
  
  verbose("\n\nFinished!\n\n")
  verbose("Best  score attained over iterations: ",global.best.s,"\n")
  verbose("Features returned in ESet: ",global.best.s.features,"\n")
  verbose("\n\n")
  
  if(best_score_only==F){
    
    #We don't just want the combination (meta-feature) at the end. We want all the features that make up the meta-feature
    #This can be obtained using the list of indices that were progressively excluded (if at all) in the step-wise procedure
    #If returning only those features that led to the best global score
    ES.best <- ES[global.best.s.features,]
    #Here, give the returned ESet an annotation based on the starting feature that gave these best results
    annotation(ES.best) <- start_feature
    
    if(length(global.best.s.features)==1){
      verbose("No meta-feature that improves the enrichment was found ..\n") 
    }
    
    #Make a list contaning two elements. 
    #The first will be the ESet with the features that gave the best meta-feature
    #The second will be the score corresponding to the meta-feature (named by the starting feature that led to that score)
    #Assign the name of the best meta-feature score to be the starting feature that gave that score
    names(global.best.s) <- start_feature 
    
    return(list("ESet"= ES.best,"Score"= global.best.s))
    
  } else{
    #Just return the score. Here we put this in the list just to support permutation-based apply functionality
    return(list(global.best.s)) }
}

#' Random permutation matrix generator
#' 
#' Produces a random permutation rank matrix given a vector of values
#' @param ord vector to be permuted. This determines the number of columns in the permutation matrix
#' @param n_perms number of permutations to generate. This determines the number of rows in the permutation matrix
#' @param seed seed which can be set for reproducibility of 'random' results. Default is 123
#' @return A row matrix of permuted values (i.e. ranks) where each row is a single permutation result
#' @export
generate_permutations<-function(ord, #These are the sample orderings to be permuted
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
null_ks<-function(ranking=NULL,
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
  
  ####### CACHE CHECKING #######
  
  if(!is.null(cache.path)){
    cat("Using provided cache root path: ",cache.path,"\n")
    setCacheRootPath(cache.path)
    
  } else{
    setCacheRootPath()
    cat("Setting cache root path as: ",getCacheRootPath(),"\n")
  }
  
  # We use the ESet, top N, score metric and seed for random permutation as the key for each cached result  
  key <- list(ES,N,null.args$metric,seed)
  
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
    if(ncores > 1 && require(doMC)){
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
    
    perm.best.scores<-unlist(alply(perm_labels_matrix,1,topn.eval,ES = ES,best.score.only=TRUE,verb=FALSE,...,.parallel=parallel,.progress = progress))
    
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
    obs.best.score <- topn.eval(ranking=ranking, #This is passed to the null_ks function
                              ES = ES,
                              best.score.only = TRUE,
                              verb=FALSE,
                              ...)
  } else{
    cat("Using provided value of observed best score ..\n\n")
  }
  
  ########### PERMUTATION P-VALUE COMPUTATION ############
  cat("Number of permutation-based scores being considered: ",length(perm.best.scores), "\n")
  #Add a smoothening factor of 1 if specified
  #This is just to not return a p-value of 0
  c=0
  if(smooth)
    c=1
  
  if(null.args$metric=="pval"){
    
    #Use negative log transform of returned search score (either computed above, or passed to the null_ks function if previously computed)
    obs.best.score <- -(log(obs.best.score))
    
    #Use negative log transform on the permuted scores (either computed above or loaded from Cache)
    perm.best.scores <- -(log(perm.best.scores))
    
  }
  
  perm.pval <- (sum(perm.best.scores > obs.best.score) + c)/(nperm + c) 
  
  
  cat("Permutation p-value: ",perm.pval,"\n\n")
  
  
  ########### END PERMUTATION P-VALUE COMPUTATION ############
  
  if(plot==T){
    #Here, let us plot the absolute values of the permutation p-values, for simplicity
    #You only consider absolute values when calculating the permutation p-values
    g <- ggplot(data=data.frame("x"=perm.best.scores),aes(x=x))+
      geom_histogram(fill="black",color="gray")+
      theme_classic()
    
    g <- g+geom_vline(xintercept=obs.best.score,linetype="longdash",size=1.5,colour="red")+
      labs(title=paste("Emperical Null distribution (N = ",length(perm.best.scores),")\n Permutation p-val <= ",round(perm.pval,5),"\nBest observed score: ",round(obs.best.score,5),sep=""),
           x="Score",
           y="Count")+
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0))
    
    print(g)
  }
  
  if(return.perm.pval)
    return(perm.pval)
}
