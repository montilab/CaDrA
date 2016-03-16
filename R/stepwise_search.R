#' Verbose control for diagnostic messages
#' 
#' Control verbosity of messages used for stepwise heuristic search tracking
#' @param ...
#' @return a message corresponding to a print statement if verbose is set to TRUE
#' @export 
verbose <- function(...){
  #Fetch verbose option set in the ks.stepwise function
  opt <- getOption("verbose",FALSE)
  if(!opt) return(invisible(NULL))
  msgs <- list(...)
  #msgs <- do.call(paste, c(msgs))
  message(msgs)
}

#' Pre-filter features
#' 
#' Pre-filter a dataset prior to running step-wise heuristic search in order to avoid testing features that are too prevalent across samples in the dataset
#' @param ES an expression set object containing binary features used for step-wise search
#' @param cutoff a fraction between 0 and 1 describing the absolute prevalence of a feature across all samples in the dataset above which the feature will be filtered out. Default is 0.6 (feature should occur in 60 percent or more of the samples to be removed)
#' @return An expression set object with the filtered-in features only
#' @export
prefilter_data<-function(ES, 
                         cutoff=0.6){
  frac<-round(rowSums(exprs(ES))/ncol(ES),2)
  cat("Pre-filtering features ..\n\n")
  cat("Retaining features having < ",cutoff*100, " % occurence in sample set..\n")
  ES<-ES[frac<cutoff,]
  cat(nrow(ES)," features retained out of ",length(frac)," supplied features\n\n")
  return(ES)
}

#' Top 'N' evaluate
#' 
#' Evaluates a list of stepwise search results from the ks.stepwise function for overlapping features. This function is mainly used to evaluate search results over the top 'N' best starting features for a given dataset.
#' @param l a list of expression sets returned by each stepwise search run of the top 'n' specified features
#' @param alt a character string specifying the alternative hypothesis, must be one of "two.sided","greater" or "less". Default is "two.sided". 
#' @param do.plot a logical indicating whether you want to only plot the resulting evaluation matrix instead of returning it. Default is TRUE
#' @param best.score.only a logical indicating whether or no to only return the best meta-feature score over the top 'n' evaluation. Default is FALSE 
#' @return Default is a binary overlap matrix (or a heatmap thereof) with rows being the union of all reported features and the column being the starting feature for each run. If best.score.only is set to TRUE, returns only the numeric score for the best meta-feature score over the top 'n' runs
#' 1's and 0's represent whether a feature in any given row is present in a meta-feature along with a starting feature in the corresponding column.
#' @export  
topn.eval<-function(l,
                    alt="two.sided", 
                    do.plot=TRUE,
                    best.score.only=FALSE){
  
  f_list<-lapply(l,featureNames)  #Get the list of feature names from each ESet
  topn_names<-names(f_list)  #Get the feature names for each top n KSS start 
  f_union<-Reduce(f = union,f_list) #Get the union of all features that were returned across all top n runs
  
  f_checklist<-lapply(f_list,function(x,ref=f_union){
    return(f_union %in% x)
  })
  
  #Make a matrix indicating which features are found across each top n run
  m<-do.call(cbind,f_checklist)*1   #Multiplying by 1 is just to convert boolean values into 1's and 0's
  rownames(m)<-f_union
  
  #Let us compute the KS scores for each top_n run so that we can arrange the results by score
  #This is just to prioritize results, in a way, instead of just clustering them heirarchically
  m.list<-lapply(l,FUN=function(x){t(matrix(ifelse(colSums(exprs(x))==0,0,1)))})
  
  #Fetch the p-values corresponding to the KS test for each top-n union
  s<-sapply(lapply(m.list,ks.genescore.mat,"less",NULL),"[[",2) # Here alternative is "less" and weight = NULL
  names(s)<-topn_names
  
  #Fetch the best score from the iterations
  best_score <- s[order(s)][1] #Based on the p-values, the lowest value will be the most sig
  
  #Order matrix in increasing order of KS score p-values
  #Add labels of which rank it was originally, and what the meta-feature p-value is
  # Here I take the logit() transform of the p-value just to avoid 0s (if p-values are too small)
  colnames(m)<-paste(colnames(m)," [",seq(1,ncol(m)),"] ",round(logit(s),3),sep="")
  m<-m[,order(s)]
  
  if(do.plot){
    heatmap.2(x = m,
              col=c("white","firebrick2"),
              Colv=FALSE,
              dendrogram="none",
              margins=c(10,10),
              cexRow=0.7,
              cexCol=0.7,
              cex.main=0.8,
              key=F,
              trace="none",
              sepwidth=c(0.1,0.1),
              sepcolor="grey90",
              colsep=1:ncol(m),
              rowsep=1:nrow(m))
    
    legend("topleft",
           legend=c("Present","Absent"),
           fill=c("firebrick2","white"),
           #legend=c("Present","Absent","Starting"),
           #fill=c("grey","white","firebrick3"),
           bty="n")
  }
  
  if(best.score.only==TRUE)
    return(best_score)
  else
    return(m)
}

#' Backward search
#' 
#' Performs a 'leave-one-out' evaluation on the existing best feature-set to see if the removal of any pre-selected feature gives a better score than that of the existing metafeature obtained in the step-wise search
#' @param ESet an ordered expression set object with the same sample ordering and features as processed by the ks.stepwise() function when performing a step-wise heuristic search
#' @param glob.f a vector containing the feature names whose union gives the best score (so far) in the search. Feature names should match those of the provided expression set object 
#' @param glob.f.ks KS score (p-value or statistic depending on chosen metric of interest) corresponding to the union of the specified vector of features 
#' @param m a character string specifying which metric to use for stepwise search criteria. One of either 'pval' or 'stat' may be used, corresponding to the KS p-value or statistic. Uses value passed in the ks.stepwise() function
#' @param alt a character string specifying the alternative hypothesis, must be one of "two.sided","greater" or "less". Uses value passed in the ks.stepwise() function
#' @param wts a vector of weights to use if performing a weighted-KS test. Default is NULL, passed in the ks.stepwise() function 
#' @return A list object with the first element being the set of features that gave a better score than the existing best score, and second being the corresponding score (to update the KS search results)
#' @export
backward_check<-function(ESet, 
                         glob.f, 
                         glob.f.ks, 
                         m,
                         alt,
                         wts){ 
  verbose("Performing backward search step..\n")
  verbose("Iterating over ",length(glob.f)," chosen features..\n")
  #print(glob.f)
  #Matrix of only global best features so far
  gmat<-exprs(ESet[glob.f,])
  rownames(gmat)<-glob.f #Not sure if this is required (might already be the case)
  #Here, we make a list that should store the features and their corresponding meta-feature KS score for each leave-one-out run
  f.names<-list()
  f.scores<-c()
  #We want to see if leaving anyone feature out improves the overall meta-feature KS score
  for(n in 1:(length(glob.f)-1)){
    
    f.names[[n]]<-glob.f[-n]
    #Take leave-one-out union of features from matrix
    u<-ifelse(colSums(gmat[-n,])==0,0,1)
    
    #Compute KS scores for this meta feature
    u.ks<-ks.genescore.mat(mat=t(matrix(u)),alt=alt,weight=wts) #Need to check if this function works for single row matrix
    score<-ifelse(m %in% "pval",sign(u.ks[1,])*u.ks[2,], u.ks[1,])
    f.scores<-c(f.scores,score)
  }
  if(m!="pval")
    f.best.index<-which.max(f.scores) #This is the index within the meta matrix
  else #If signed pvalues
    f.best.index<-order(-sign(f.scores),f.scores)[1] #Top p-value ordered by sign and numerical value; #This is the index within the meta matrix
  f.best.score<-f.scores[f.best.index]
  #Check if any one of the computed scores has a better score than the entire meta-feature's score
  if(ifelse(m=="pval", (sign(f.best.score)>0 & abs(f.best.score) < abs(glob.f.ks)),f.best.score > glob.f.ks)){
    verbose("Found improvement on removing existing feature..\n")
    #Return the set of features that gave a better score than the existing best score, and the score as well
    return(list(f.names[[f.best.index]],f.best.score))  
  } else{
    #Don't change anything. Return the existing best set of features and the corresponding KS score
    return(list(glob.f,glob.f.ks))
  }
  
}

#' Kolmogorov-Smirnov step-wise search
#' 
#' Performs step-wise heuristic search using an ordered set of binary features to see whether there are features whose union is more skewed (enriched at the extremes) than either features alone. This is the main functionality of the CaDrA package.
#' @param ranking vector containing rankings for sample ordering. Default is NULL. If NULL, we assume the ESet is already ranked  
#' @param ES an expression set object of binary features. The first column of the featureData for the expression set must contain the names of the corresponding features, which are used in the search   
#' @param metric a character string specifying which metric to use for stepwise search criteria. One of either 'pval' or 'stat' may be used, corresponding to the KS p-value or statistic
#' @param top_score_start a logical indicating whether or not to start the step-wise search with the feature having the best score based on the specified metric of interest. Default is TRUE
#' @param back_search a logical indicating whether or not to perform a forward-backward search (i.e. remove features along the search if it improves the best score). Default is TRUE. Uses function backward_check() 
#' @param cust_start an integer specifying a specific index within the expression set object of the feature to start the step-wise search with. Default is NULL 
#' @param decr a logical indicating whether samples are sorted in decreasing order of phenotype of interest. This affects how the ks statistic is interpreted
#' @param best_score_only a logical indicating whether or not the function should return only the score corresponding to the search results
#' @param alt a character string specifying the alternative hypothesis, must be one of "two.sided","greater" or "less". Default is "two.sided"
#' @param wts a vector of weights to use if performing weighted-KS testing. Default is NULL
#' @param verb a logical indicating whether or not to print diagnostic messages. Default is FALSE 
#' @return If best_score_only is set to TRUE, this function returns a list object with the score corresponding to the union of the search meta-feature. If this is set to FALSE, an expression set object containing the features whose union gave the best score is returned. 
#' @export
ks.stepwise<-function(ranking=NULL,  
                      ES, 
                      metric, 
                      top_score_start=TRUE, 
                      back_search=TRUE, 
                      cust_start=NULL,
                      decr=TRUE, #Set this to true if ASSIGN score sorting is done in decreasing order of ASSIGN scores. False otherwise
                      #Reverses direction of KS statistic, and hence, direction of feature sorting 
                      best_score_only=FALSE, 
                      alt="two.sided", 
                      wts=NULL,
                      verb=FALSE
){
  
  #Setup verbose option definition
  options(verbose=verb)
  
  #Use ranking to re-order samples 
  if(!(is.null(ranking))){
    if(length(ranking)!=ncol(ES))
      stop("Ranking variable has to be of the same length as the number of samples\n\n")
    verbose("Using provided ordering to re-rank samples..\n")
    ES<-ES[,ranking]
  }
  #Compute row-wise directional KS scores for existing (raw/starting) binary features in ESet
  if(!(is.null(wts)))
    verbose("Using weighted method for KS testing using provided weights..\n")
  
  ks<-ks.genescore.mat(mat=exprs(ES),alt=alt,weight=wts)
  ks.stat<-ks[1,]
  ks.pval<-ks[2,]
  
  #Define scores based on specified metric of interest
  
  if(!metric %in% c('stat','pval'))
    stop("Please specify metric parameter as either 'stat' or 'pval' to use for search..\n")
  score<-ifelse(rep(metric,nrow(ES)) %in% "pval",sign(ks.stat)*ks.pval, ks.stat)
  
  
  verbose("Using ",metric," as measure of improvement measure ..\n\n")
  
  
  
  #ADDED
  #Order ESet in increasing/decreasing order of user-defined score (ks stat or p-val)
  #This comes in handy when doing the top-N evaluation of the top N 'best' features
  
  score.rank<-if (metric!="pval") order(score) else order(-sign(score),score)
  verbose("Ranking ESet features by metric..\n")
  ES<-ES[score.rank,]
  score<-score[score.rank]
  ##
  
  #Here, we will assume ASSIGN scores have samples ranked in decreasing order
  #This means an higher/positive KS stat is associated with a higher ASSIGN score 
  
  #Let us start with the first (top ranked) feature
  #Fetch index of feature having best KS score. We start here
  if(top_score_start==T){
    verbose("Starting with feature having best ranking ..\n")
    #MODIFIED
    #best.ks.index<-ifelse(metric!="pval",which.max(score),order(-sign(score),score)[1]) #This assumes that samples are ordered in decreasing order of ASSIGN score
    best.ks.index<-1
    ##
    best.ks<-score[best.ks.index]
    
  } else {
    if(!is.null(cust_start)) {
      verbose("Starting with specified sample feature ..\n")
      best.ks.index<-cust_start
      best.ks<-score[best.ks.index]		
    } else {
      #Throw an error message saying you have to specify some starting criteria for the search
      stop("top_score_start is set to FALSE but no starting index specified.\n Please specify starting index of the feature to begin search..\n")
    }
  }
  #Print the featureData for this starting point feature so that we are aware
  start_feature <- as.character(fData(ES)[best.ks.index,1])
  verbose("Feature: ",start_feature,"\n")
  verbose("Score: ",best.ks,"\n")
  
  #Fetch the vector corresponding to best score
  #Set this as the initial 'meta-feature'
  best.meta<-as.numeric(exprs(ES)[best.ks.index,])
  
  #This is just so that it enters the while loop for the first iteration
  new.best.ks<-best.ks
  
  ###### INITIALIZE VARIABLES ###########
  #######################################
  
  #counter variable for number of iterations
  i=0
  #Parameter for number of continuous mistakes
  b=0
  #Parameter for total number of mistakes
  #c=0 #We don't really care about this (set it to abnormally large number)
  
  #Variable to store best KS score attained over all iterations
  #initialize this to the starting best KS score
  global.best.ks<-best.ks
  #Variable to store indices to make meta-feature with best KS score attained over all iterations
  #initialize this to the starting feature's index
  global.best.ks.index<-best.ks.index
  global.best.ks.features<-c(as.character(fData(ES)[best.ks.index,1]))
  
  
  ###### BEGIN ITERATIONS ###############
  #######################################
  
  verbose("\n\nBeginning stepwise search..\n\n")
  
  #This condition defines the search criteria
  while (ifelse(metric=="pval", (sign(new.best.ks)>0 & abs(new.best.ks) < abs(best.ks)),new.best.ks > best.ks) | b < 2){
    verbose("\n\n")
    verbose("Iteration number ",i+1," ..\n")
    
    
    if(i!=0){
      #Not the first iteration	
      
      verbose("New best ks: ",new.best.ks,"\n")
      verbose("Best ks: ",best.ks,"\n")
      
      if(ifelse(metric=="pval", (sign(new.best.ks)>0 & abs(new.best.ks) < abs(best.ks)),new.best.ks > best.ks))
      {
        verbose("Found feature that improves KS score!\n")
        #Now that we have an improvement, let us reset the 'continuous mistake counter' to 0
        b=0
      }
      new.best.meta<-meta.mat[hit.best.ks.index,]
      
      #Add that index to the group of indices to be excluded for subsequent checks
      #Here we go off the rownames to find which index to exclude from the ESet
      best.ks.index<-c(best.ks.index,which(rownames(ES)==best.feature))
      
      #If performing a forward-backward search, we need to check if adding this last feature
      #works better when leaving any of the existing best features out
      #This is only useful if you have at least 4 features
      if(length(global.best.ks.features) > 3 & back_search==T){
        
        backward_search.results<-backward_check(ESet=ES,
                                                glob.f=global.best.ks.features, #Global feature set so far
                                                glob.f.ks=global.best.ks, #KS score corresponding to this global feature set
                                                m=metric,
                                                alt=alt,
                                                wts=wts)
        global.best.ks.features<-backward_search.results[[1]]
        global.best.ks<-backward_search.results[[2]]
      }
      
      #Reset current minimum values and feature to new minimum values and feature
      best.meta<-new.best.meta
      best.ks<-new.best.ks
      verbose("KS score: ", best.ks,"\n")
    }
    
    
    #Take the OR function between that feature and all other features, to see which gives the best KS score
    #Keep in mind, the number of rows in meta.mat keeps reducing by one each time we find a hit that improves the KS score
    verbose("Forming meta-feature matrix with all other features in dataset..\n")
    meta.mat<-sweep(exprs(ES)[-best.ks.index,],2,best.meta,`|`)*1
    
    #verbose("Number of rows in meta-feature matrix: ",nrow(meta.mat),"\n")
    #print(rownames(meta.mat))
    
    
    #Let us see if any of the newly formed meta features are the same as the feature itself
    #Leave out these metafeatures as they add no new information
    
    #####<Still need to implement this functionality>######
    
    
    #With the newly formed 'meta-feature' matrix, compute directional KS scores and choose the feature that gives the best KS score
    #Compute row-wise directional KS scores for existing (raw/starting) binary features in ESet
    ks<-ks.genescore.mat(mat=meta.mat,alt=alt,weight=wts)
    ks.stat<-ks[1,]
    ks.pval<-ks[2,]
    
    
    scores<-ifelse(rep(metric,nrow(meta.mat)) %in% "pval",sign(ks.stat)*ks.pval,ks.stat)
    
    
    
    #Find index of feature that gives lowest ks score when combined with chosen starting feature
    if(metric!="pval")
      hit.best.ks.index<-which.max(scores) #This is the index within the meta matrix
    else #If signed pvalues
      hit.best.ks.index<-order(-sign(scores),scores)[1] #Top p-value ordered by sign and numerical value; #This is the index within the meta matrix
    new.best.ks<-scores[hit.best.ks.index] #This is within the meta matrix
    
    #Diagnostic
    #verbose("Best score from meta matrix: ",new.best.ks,"\n")
    
    
    #Find which feature produced that score, in combination with meta feature used
    best.feature<-rownames(meta.mat)[hit.best.ks.index]
    verbose("Feature that produced best score in combination with previous meta-feature: ",best.feature,"\n")
    #Diagnostics
    #verbose("Meta mat index: ",hit.best.ks.index,"\n")
    #verbose("Estimated ES index for feature: ", which(rownames(ES)==best.feature),"\n")
    #verbose("Actual index of feature in ES: ",which(rownames(ES)==best.feature),"\n")
    
    
    #Checking whether it's a new best score
    if(ifelse(metric=="pval", (sign(new.best.ks)>0 & abs(new.best.ks) < abs(global.best.ks)),new.best.ks > global.best.ks))
      #if(new.best.ks<global.best.ks)
      #Update best score
    { global.best.ks<-new.best.ks
      global.best.ks.index<-c(global.best.ks.index,hit.best.ks.index+(nrow(ES)-nrow(meta.mat))) #Here we add the difference in rows because it is relative 
      global.best.ks.features<-c(global.best.ks.features,best.feature)
    }
    
    verbose("Global best score so far: ",global.best.ks,"\n")
    #to the entire ES which will have x more rows than the meta matrix 
    #Let us set a parameter that can allow for a few mistakes
    #If the ks score isn't improving for the metafeature, add a counter (constrained in the while loop entry)
    if(ifelse(metric=="pval", sign(new.best.ks)<0 | (abs(new.best.ks) >= abs(best.ks)),new.best.ks <= best.ks)){
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
  verbose("Best KS score attained over iterations: ",global.best.ks,"\n")
  verbose("Features returned in ESet: ",global.best.ks.features,"\n")
  verbose("\n\n")
  
  if(best_score_only==F){
    
    #We don't just want the combination (meta-feature) at the end. We want all the features that make up the meta-feature
    #This can be obtained using the list of indices that were progressively excluded (if at all) in the step-wise procedure
    #If returning only those features that led to the best global score
    ES.best<-ES[global.best.ks.features,]
    #Here, give the returned ESet an annotation based on the starting feature that gave these best results
    annotation(ES.best) <- start_feature
    
    if(length(global.best.ks.features)==1){
      verbose("No meta-feature that improves the enrichment was found ..\n") 
    }
    
    return(ES.best)
  } else{
    return(list(global.best.ks)) }
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
  return(perm) 
}

#' Permutation-based step-wise searching
#' 
#' Performs permutation-based significance testing of step-wise KS search.
#' @param ranking a vector containing rankings for sample ordering. Default is NULL. If NULL, we assume the ESet is already ranked. This is used only to compute the observed KS stepwise search score  
#' @param ES an ordered expression set object with the same sample ordering and features as processed by the ks.stepwise() function when performing a step-wise heuristic search
#' @param metric a character string specifying which metric to use for stepwise search criteria. One of either 'pval' or 'stat' may be used, corresponding to the KS p-value or statistic
#' @param nperm an integer specifying the number of permutations to perform. Default = 1000
#' @param back_search a logical indicating whether or not to perform a forward-backward search (i.e. remove features along the search if it improves the best score). Default is TRUE. Uses function backward_check() 
#' @param cust_start an integer specifying a specific index within the expression set object of the feature to start the step-wise search with. Default is NULL
#' @param plot logical indicating whether or not to plot the emperical null distribution with the observed p-value
#' @param alt a character string specifying the alternative hypothesis, must be one of "two.sided","greater" or "less". Default is "two.sided"
#' @param wts a vector of weights to use if performing weighted-KS testing. Default is NULL 
#' @param smooth logical indicating whether or not to smoothen the p-value calculation to avoid p-value of 0. Default is TRUE 
#' @param seed seed set for permutation. Default = 123
#' @param ncores number of cores to use, if using parallelization. Default = 1
#' @return a list object containing the permutation-based p-value based on the generated empirical null distribution of stepwise KS searching
#' @export
null_ks<-function(ranking=NULL,
                  ES,
                  metric,
                  nperm=1000,
                  back_search=TRUE,
                  cust_start=NULL,
                  plot=TRUE,
                  alt="two.sided",
                  wts=NULL,
                  smooth=TRUE,
                  seed=123,
                  ncores=1){
  
  cat("\n\n\nBEGINNING PERMUTATION-BASED SIGNIFICANCE TESTING\n\n\n")
  
  
  if(!is.null(cust_start)) 
    tss<-FALSE
  else
    tss<-TRUE
  
  
  ################################# UNDER CONSTRUCTION #########################################
  # Sets up the parallel backend which will be used by Plyr.
  parallel = F
  progress = "none"
  if(ncores > 1 && require(doMC)){
    registerDoMC(ncores)
    parallel = T
    verbose("Running tests in parallel..\n")
  } else {
    registerDoSEQ()
    progress = "text"
  }
  
  perm_labels_matrix<-generate_permutations(ord=seq(1,ncol(ES)),
                                            n_perms=nperm,
                                            seed=seed)
  
  #Set verbose to FALSE since we don't want to print anything
  options(.null_ks.verb=FALSE)
  
  cat("Computing observed best score ..\n\n")
  obs.best_score<-unlist(ks.stepwise(ranking=ranking, #Here, we would pass the SAME ranks specified in the null_ks function, intended for the observed search results
                                     ES = ES,
                                     metric=metric,
                                     back_search=back_search,
                                     top_score_start = tss,
                                     cust_start=cust_start,
                                     best_score_only = T,
                                     alt=alt,
                                     wts=wts,
                                     verb=FALSE))
  #Use logit transform
  obs.best_score<-logit(obs.best_score)
  
  cat("Computing permutation-based scores ..\n\n")
  
  ptm<-proc.time()
  perm.best_scores<-unlist(alply(perm_labels_matrix,1,ks.stepwise,ES = ES,metric=metric,top_score_start = tss,cust_start=cust_start,best_score_only=T,alt=alt,wts=wts,verb=F,.parallel=parallel))

  #Use logit transform
  perm.best_scores <- logit(perm.best_scores)

  cat("FINISHED\n")
  cat("Time elapsed: ",round((proc.time()-ptm)[3]/60,2)," mins \n\n")
  ############################################################################################## 
  
  
  #Add a smoothening factor of 1 if specified
  #This is just to not return a p-value of 0
  c=0
  if(smooth)
    c=1
  if(metric=="pval")
    #perm.pval<-(sum(abs(perm.best_scores) < abs(obs.best_score)) + c)/(nperm + c) 
    perm.pval<-(sum(perm.best_scores < obs.best_score) + c)/(nperm + c) 
  else
    perm.pval<-(sum(perm.best_scores > obs.best_score) + c)/(nperm + c) 
  
  
  cat("Permutation p-value: ",perm.pval,"\n\n")

  if(plot==T){
    #Here, let us plot the absolute values of the permutation p-values, for simplicity
    #You only consider absolute values when calculating the permutation p-values
    g<-ggplot(data=data.frame("x"=perm.best_scores),aes(x=x))+
      geom_histogram(fill="black",color="gray")+
      theme_classic()
    
    g<-g+geom_vline(xintercept=obs.best_score,linetype="longdash",size=1.5,colour="red")+
      labs(title=paste("Emperical Null distribution (N = ",nperm,")\n Permutation p-val <= ",round(perm.pval,5),"\nBest observed score: ",round(obs.best_score,5),sep=""),
           x="Score",
           y="Count")+
      #scale_x_continuous(expand = c(0, 0),limits=c(0,1)) + scale_y_continuous(expand = c(0, 0))
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0))
    
    print(g)
  }
  
}
