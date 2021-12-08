
#' CaDrA Permutation-Based Candidate Searching
#' 
#' Performs permutation-based significance testing of candidate search results.
#' @param ES an expression set object of binary features (required). It can be a BioBase expressionSet object or an expression matrix. The rownames or featureData of the expression set must contain the names of the corresponding features which are used in the search.   
#' @param input_score a vector containing score for sample ordering (required). 
#' @param max_size an integer specifying the maximum size a meta-feature can extend to do for a given search. Default is 7
#' @param method a character string specifying the method used to compute scores for features, must be one of "ks" or "wilcox" or "mi" (mutually exclusive method from REVEALER) or "custom" (a personal customization method). If input_score contains ranked scores, then 'ks' method is used by default. Otherwise, 'mi" is the default method
#' @param metric a character string specifying which metric to use for candidate search. One of either 'pval' or 'stat' may be used, corresponding to the score p-value or statistic. Default is 'pval'
#' @param search_start an integer specifying a specific index within the expression set object of the features to start the candidate search. Default is 1 where the search starts with the top ranked feature. If an integer N is specified (where N < nrow(dataset), then the search starts with the Nth best feature. If a string is specified, then the search starts with the feature with this name (must be a valid rowname in the Eset)
#' @param search_method a character string specifying which method to perform or filter out the best candidates. Default is 'forward'.
#' @param best_score_only a logical indicating whether or not the function should return only the score corresponding to the search results. Default is FALSE
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided", "greater" or "less". Default is "less" for left-skewed significance testing.
#' @param weights an integer vector of weights to use if performing weighted-KS testing. Default is NULL. Value passed to compute_score() function  
#' @param ranks an integer vector of sample rankings to use if performing Wilcoxon rank sum testing. Default is NULL. If NULL, then samples are assumed to be ordered by increasing ranking. Value passed to compute_score() function 
#' @param verbose a logical indicating whether or not to print diagnostic messages. Default is TRUE 
#' @param outdir output the results to a desired directory if a path is specified. Default is NULL
#' @return If best.score.only is set to TRUE, this function returns a list object with the score corresponding to the union of the search meta-feature. If this is set to FALSE, a list containing both the ESet object pertaining to the returned meta-feature as well as the corresponding score  is returned. 
#' @export
#' @import Biobase 
CaDrA <- function(
  ES, 
  input_score, 
  method = c("ks", "wilcox", "mi"),
  metric = "pval",
  max_size = 7,
  search_start = 1,
  search_method = c("forward", 'backward', "both"), 
  best_score_only = FALSE, 
  alternative = "less", 
  weights = NULL,
  ranks = NULL,
  verbose = TRUE,
  outdir = NULL
){
  # # cust_function a customize function that computes the score. It will be used when the 'method' is set to 'custom'
  # 
  # # Set up verbose option
  # options(verbose=verbose)
  # 
  # #Check if the ES is provided
  # if(length(ES) == 0){
  # 
  #   stop("ES was not provided (required).\n")
  # 
  # }else{
  # 
  #   # Check if ES is an expression matrix or an expressionSet from BioBase
  #   es_check <- tryCatch({
  # 
  #     # Get expression matrix
  #     eset <- exprs(ES)
  # 
  #     return("ExpressionSet")
  # 
  #   }, error=function(err){
  # 
  #     # if not an expressionSet from BioBase, check if it is an expression matrix
  #     if(is.matrix(ES)){ return("ExpressionMatrix") } else { return(FALSE) }
  # 
  #   })
  # 
  # }
  # 
  # # If expressionSet from BioBase was provided, extract the expression matrix only
  # if(es_check == "ExpressionSet"){
  #   eset_matrix = exprs(ES)
  # }else if(es_check == "ExpressionMatrix"){
  #   eset_matrix = as.matrix(ES)
  # }else{
  #   stop("Invalid ES input. ES can be an expression set from BioBase or an expression matrix.\n")
  # }
  # 
  # # Make sure the input ES has rownames for features
  # if(is.null(rownames(eset_matrix)))
  #   stop("ES object provided does not have rownames or featureData to track the features by. Please provide unique feature names as rownames of the input ES.\n")
  # 
  # # Remove duplicated rownames or features
  # eset_matrix = eset_matrix[!duplicated(rownames(eset_matrix)),]
  # 
  # ###### SAMPLE RANKING OR NON-RANKING #####
  # 
  # # Check input_score is provided
  # if(length(input_score) == 0){
  # 
  #   stop("input_score was not provided (required).\n")
  # 
  # }else{
  # 
  #   # Make sure the input_score has the same length as number of samples in ES
  #   if(length(input_score) != ncol(eset_matrix))
  #     stop("input_score variable must have the same length as the number of samples in ES.\n")
  # 
  #   # check if input_score are ranked values or computed values
  #   ranking = ifelse(all(sort(input_score) == 1:ncol(eset_matrix)), TRUE, FALSE)
  # 
  #   # if(ranking){
  #   # 
  #   #   verbose("A vector of ranked input_score was provided. Using the ordering to re-rank the samples..\n")
  #   #   eset_matrix <- eset_matrix[,ranking]
  #   # 
  #   # }else{
  #   # 
  #   #   verbose("A vector of input_score values was provided. Using the computed values to re-ordering the samples..\n")
  #   # 
  #   #   ordering <- data.frame(position=1:length(input_score), input_score=input_score) %>% dplyr::arrange(desc(input_score))
  #   #   eset_matrix <- eset_matrix[,ordering]
  #   # 
  #   # }
  # 
  # }
  # 
  # # Check the method
  # if(length(method) > 0 && method %in% method_options){
  # 
  #   # Compute row-wise directional KS scores for existing (raw/starting) binary features in ES
  #   if(method == "ks"){
  #     verbose("Using KS method for feature scoring..\n")
  #   }
  # 
  #   # Compute row-wise Wilcox rank sum scores for existing (raw/starting) binary features in ES
  #   if(method == "wilcox"){
  #     verbose("Using Wilcoxon method for feature scoring..\n")
  #   }
  # 
  #   # Compute mutually exclusive method for existing (raw/starting) binary features in ES
  #   if(method == "mi"){
  #     verbose("Using Mutually Exclusive method for feature scoring.\n")
  #   }
  # 
  #   # Other future methods can be implemented here and its verbose message go here
  # 
  # 
  # } else {
  # 
  #   stop(paste0("Invalid method specified. The method can be ", paste0(paste0("'", method_options, "'"), collapse=","), "\n"))
  # 
  # }
  # 
  # # Select the appropriate method to compute scores based on skewdness of a given binary matrix
  # s <- switch(
  #   method,
  #   #ks = ks_genescore_mat(
  #   #  mat = eset_matrix,
  #   #  alt = alternative,
  #   #  weight = weights
  #   #),
  #   wilcox = wilcox_genescore_mat(
  #     mat = eset_matrix,
  #     alternative = alternative,
  #     rank = ranks
  #   )#,
  #   ##mi = revealer_genescore_mat(
  #   #  mat = eset_matrix,
  #   #  alt = alternative
  #   #)
  # )
  # 
  # # Score returned by either ks or wilcox-based functions
  # s.stat <- s[1,]
  # s.pval <- s[2,]
  # 
  # # Define scores based on specified metric of interest
  # 
  # if(!metric %in% c('stat','pval'))
  #   stop("Please specify metric parameter as either 'stat' or 'pval' to use for search..\n")
  # 
  # score <- ifelse(rep(metric,nrow(ES)) %in% "pval", sign(s.stat)*s.pval, s.stat)
  # verbose("Using ",metric," as measure of improvement measure ..\n\n")
  # 
  # ###### FEATURE PRE-RANKING #####
  # ##################################
  # 
  # # Re-order ESet in decreasing order of user-defined score (s stat or p-val)
  # # This comes in handy when doing the top-N evaluation of the top N 'best' features
  # 
  # score.rank <- if (metric!="pval") order(score) else order(-sign(score),score)
  # verbose("Ranking ESet features by metric..\n")
  # 
  # ES <- ES[score.rank,]
  # score <- score[score.rank]
  # 
  # if(is.null(cust_start)){
  #   verbose("Starting with feature having best ranking ..\n")
  #   best.s.index <- 1
  # } else {
  #   if(is.numeric(cust_start)){
  #     # User-specified feature index (has to be an integer from 1:nrow(ES))
  #     verbose("Starting with specified sorted feature index ..\n")
  # 
  #     if(cust_start > nrow(ES)) # Index out of range
  #       stop("Invalid starting index specified.. Please specify a valid starting index within the range of the existing ESet..\n")
  # 
  #     best.s.index <- cust_start
  #   }
  # 
  #   if(is.character(cust_start)){
  #     # User-specified feature name (has to be a character from rownames(1:nrow(ES)))
  #     verbose("Starting with specified feature name ..\n")
  # 
  #     if(!(cust_start %in% rownames(ES))) #provided feature name not in rownames
  #       stop("Provided starting feature does not exist among ESet's rownames .. please check stepwise.search cust_start parameter options for more details ..\n\n")
  # 
  #     best.s.index <- which(rownames(ES)==cust_start)
  #   } # end if is.character
  # } # end else (!is.null)
  # 
  # ###### INITIALIZE VARIABLES ###########
  # #######################################
  # 
  # # Print the featureData for this starting point feature so that we are aware
  # # Here we assume the ESet's fData is included as rownames
  # #start.feature <- as.character(fData(ES)[best.s.index,1])
  # start.feature <- rownames(ES)[best.s.index]
  # best.feature <- start.feature
  # best.s <- score[best.s.index]
  # 
  # verbose("Feature: ",start.feature,"\n")
  # verbose("Score: ",best.s,"\n")
  # 
  # #Fetch the vector corresponding to best score
  # #Set this as the initial 'meta-feature'
  # best.meta <- as.numeric(exprs(ES)[best.s.index,])
  # 
  # #counter variable for number of iterations
  # i=0
  # 
  # #Variable to store best score attained over all iterations
  # #initialize this to the starting best score
  # global.best.s <- best.s
  # 
  # # Vector of features in the (growing) obtained metafeature. Begin with just the starting feature
  # global.best.s.features <- c()
  # 
  # ###### BEGIN ITERATIONS ###############
  # #######################################
  # 
  # verbose("\n\nBeginning stepwise search..\n\n")
  # 
  # while ((ifelse(metric=="pval", (sign(best.s) > 0 & (abs(best.s) < abs(global.best.s))),best.s > global.best.s) | i == 0) & (length(global.best.s.features) < max.size)){
  #   verbose("\n\n")
  #   verbose("Iteration number ",(i+1)," ..\n")
  #   verbose("Global best score: ",global.best.s,"\n")
  #   verbose("Previous score: ",best.s,"\n")
  # 
  #   # Update scores and feature set since since entry into the loop means there is an improvement (iteration > 0)
  #   global.best.s <- best.s
  #   global.best.s.features <- c(global.best.s.features,best.feature)
  # 
  #   verbose("Current feature set: ",global.best.s.features,"\n")
  # 
  #   if(i!=0){
  #     verbose("Found feature that improves  score!\n")
  #     # Update the new best meta feature (from meta mat)
  #     best.meta <- meta.mat[hit.best.s.index,]
  #     #Add that index to the group of indices to be excluded for subsequent checks
  #     #Here we go off the rownames in the original matrix to find which index to exclude from the ESet in subsequent iterations
  #     best.s.index <- c(best.s.index,which(rownames(ES)==best.feature))
  #   }
  # 
  #   # Perform a backward check on the list of existing features and update global scores/feature lists accordingly
  #   if(length(global.best.s.features) > 3 & back_search==TRUE){
  # 
  #     backward_search.results <- backward_check(ESet=ES,
  #                                               glob.f=global.best.s.features, #Global feature set so far
  #                                               glob.f.s=global.best.s, # score corresponding to this global feature set
  #                                               m=metric,
  #                                               method=method,  # passed to compute_score() function
  #                                               alt=alt,        # passed to compute_score() function
  #                                               wts=wts,        # passed to compute_score() function
  #                                               rnks=rnks)      # passed to compute_score() function
  #     # Update globlal features, scores
  #     global.best.s.features <- backward_search.results[[1]]
  #     global.best.s <- backward_search.results[[2]]
  #     # Update best.meta based on feature set
  #     best.meta <- as.numeric(ifelse(colSums(exprs(ES)[global.best.s.features,])==0,0,1))
  #   }
  # 
  #   #Take the OR function between that feature and all other features, to see which gives the best  score
  #   #Keep in mind, the number of rows in meta.mat keeps reducing by one each time we find a hit that improves the  score
  #   verbose("Forming meta-feature matrix with all other features in dataset..\n")
  #   # Here "*1" is used to convert the boolean back to integer 1's and 0's
  #   # Notice we remove anything in best.s.index from the original matrix first, to form the meta matrix.
  #   meta.mat <- sweep(exprs(ES)[-best.s.index,],2,best.meta,`|`)*1
  # 
  # 
  #   # Check if there are any features that are all 1's generated on taking the union
  #   # We cannot compute statistics for such features and they thus need to be filtered out
  #   if(any(rowSums(meta.mat)==ncol(meta.mat))){
  #     warning("Features with all 1's generated upon taking matrix union .. removing such features before progressing..\n")
  #     meta.mat <- meta.mat[rowSums(meta.mat) != ncol(meta.mat),]
  #   }
  # 
  # 
  # 
  #   #With the newly formed 'meta-feature' matrix, compute directional  scores and choose the feature that gives the best  score
  #   #Compute row-wise directional  scores for existing (raw/starting) binary features in ESet
  #   s <- compute_score(mat=meta.mat,method=method,alt=alt,weight=wts,ranks=rnks)
  #   s.stat <- s[1,]
  #   s.pval <- s[2,]
  # 
  # 
  #   # Take signed pval or stat depending on user-defined metric
  #   # This will be the same length as nrow(meta.mat)
  #   scores <- ifelse(rep(metric,nrow(meta.mat)) %in% "pval",sign(s.stat)*s.pval,s.stat)
  # 
  # 
  #   #Find index of feature that gives lowest s score when combined with chosen starting feature
  #   if(metric!="pval"){
  #     hit.best.s.index <- which.max(scores) #This is the index within the meta matrix
  #   } else { #If signed pvalues
  #     hit.best.s.index <- order(-sign(scores),scores)[1] #Top p-value ordered by sign and numerical value; #This is the index within the meta matrix
  #   }
  # 
  #   best.s <- scores[hit.best.s.index] #This is the best score from the meta matrix
  # 
  #   # Find which feature produced that score, in combination with meta feature used
  #   # We go from index to rowname space here in the meta matrix
  #   # We can do this because rownames are preserved between the original and meta features on using sweep()
  #   best.feature <- rownames(meta.mat)[hit.best.s.index]
  #   verbose("Feature that produced best score in combination with previous meta-feature: ",best.feature,"\n")
  #   verbose("Score: ",best.s,"\n")
  # 
  #   # If no improvement (exiting loop)
  #   if(ifelse(metric=="pval", sign(best.s) < 0 | (abs(best.s) >= abs(global.best.s)),best.s <= global.best.s)){
  #     verbose("No further improvement in score has been found..\n")
  #   }
  # 
  #   #Increment counter
  #   i=i+1
  # 
  # } #########End of while loop
  # 
  # verbose("\n\n")
  # verbose("\n\nFinished!\n\n")
  # verbose("Number of iterations covered: ",i,"\n")
  # verbose("Best  score attained over iterations: ",global.best.s,"\n")
  # if(length(global.best.s.features)==1){
  #   warning("No meta-feature that improves the enrichment was found ..\n")
  # }
  # 
  # verbose("Features returned in ESet: ",global.best.s.features,"\n")
  # verbose("\n\n")
  # 
  # if(best.score.only==FALSE){
  # 
  #   #We don't just want the combination (meta-feature) at the end. We want all the features that make up the meta-feature
  #   #This can be obtained using the list of indices that were progressively excluded (if at all) in the step-wise procedure
  #   #If returning only those features that led to the best global score
  #   ES.best <- ES[global.best.s.features,]
  #   #Here, give the returned ESet an annotation based on the starting feature that gave these best results
  #   annotation(ES.best) <- start.feature
  #   colnames(ES.best) <- colnames(ES)
  # 
  #   #Make a list contaning two elements.
  #   #The first will be the ESet with the features that gave the best meta-feature
  #   #The second will be the score corresponding to the meta-feature (named by the starting feature that led to that score)
  #   #Assign the name of the best meta-feature score to be the starting feature that gave that score
  #   names(global.best.s) <- start.feature
  # 
  #   return(list("ESet"= ES.best,"Score"= global.best.s))
  # 
  # } else{
  #   #Just return the score. Here we put this in the list just to support permutation-based apply functionality
  #   return(list(global.best.s))
  # }

}  







