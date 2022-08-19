
verbose <- function(...){
  
  #Fetch verbose option set in the stepwise.search() function
  opt <- getOption("verbose",FALSE)
  if(!opt) return(invisible(NULL))
  msgs <- list(...)
  #msgs <- do.call(paste, c(msgs))
  message(msgs,"\n")
  
}

#' Pre-filter features
#' 
#' Pre-filter a dataset prior to running step-wise heuristic search in order to 
#' avoid testing features that are too prevalent or too sparse across samples in
#' the dataset
#' @param ES an expression set object containing binary features used for 
#' step-wise search
#' @param max.cutoff a numeric value between 0 and 1 describing the absolute 
#' prevalence of a feature across all samples in the dataset above which the 
#' feature will be filtered out. Default is 0.6 (feature that occur in 
#' 60 percent or more of the samples will be removed)
#' @param min.cutoff a numeric value between 0 and 1 describing the absolute 
#' prevalence of a feature across all samples in the dataset below which the 
#' feature will be filtered out. Default is 0.03 (feature that occur in 
#' 3 percent or less of the samples will be removed)
#' @param verbose a logical value indicates whether or not to print the 
#' diagnostic messages. Default is \code{FALSE}. 
#' @return An expression set object with only the filtered-in features 
#' given the filter thresholds specified
#' @examples
#' data(sim.ES)
#' 
#' # Filter out features having < 3 and > 60% prevalence across all samples 
#' # (default)
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
prefilter_data <- function(
  ES, 
  max.cutoff=0.6,
  min.cutoff=0.03,
  verbose=FALSE
){
  
  options(verbose = verbose)
  
  # Compute the frequency of feature occurence across all samples  
  # (i.e. fraction of samples having the feature)
  frac <- round(rowSums(exprs(ES))/ncol(ES),2)
  
  verbose("Pre-filtering features ..\n\n")
  verbose("Removing features having < ",min.cutoff*100, "and > ",
          max.cutoff*100, " % occurence in sample set..\n")
  
  ES <- ES[ (frac >= min.cutoff) & (frac <= max.cutoff) , ]
  
  verbose(nrow(ES)," features retained out of ",length(frac),
          " supplied features in dataset\n\n")
  
  return(ES)
  
}


#' ks_test_d_wrap_ wrapper
#'
#' Compute directional Kolmogorov-Smirnov scores
#' @param n_x length of ranked list
#' @param y positions of geneset items in ranked list (ranks)
#' @param alt alternative hypothesis for p-value calculation
#' @keywords internal
#' @useDynLib CaDrA ks_test_d_wrap_
#' 
#' @return need a return value here
ks_test_double_wrap <- function(n_x, y, alt="less") {
  
  if(length(alt) > 0){
    alt_int<- switch(alt, two.sided=0L, less=1L, greater=-1L, 1L)
  } else {
    alt_int <- 1L
  }
  
  # If input is an empty vector
  if( n_x != length(y) | length(y) < 1) return (NULL)
  
  res <- .Call(ks_test_d_wrap_,  as.integer(n_x), as.numeric(y), alt_int)
  res
  
}

#' ks_plot wrapper
#'
#' Return a dataframe from ks_genescore function
#' @param n_x length of ranked list
#' @param y positions of geneset items in ranked list (ranks)
#' @param weight a vector of weights 
#' @param alt alternative hypothesis for p-value calculation
#' @keywords internal
#' @useDynLib CaDrA ks_plot_wrap_
#' 
#' @return need a return value here
ks_plot_wrap <- function(n_x, y, weight, alt="less") {
  
  if(length(alt) > 0){
    alt_int<- switch(alt, two.sided=0L, less=1L, greater=-1L, 1L)
  } else {
    alt_int <- 1L
  }
  y <- as.integer(y)
  n_x <- as.integer(n_x)
  res <- .Call(ks_plot_wrap_, n_x, y, weight, alt_int)
  res <- res[!is.na(res$X), ]
  res
  
}



#' ks.genescore wrapper
#'
#' Compute directional Kolmogorov-Smirnov scores for each row of a given vector
#' @param n_x length of ranked list
#' @param y positions of geneset items in ranked list (ranks)
#' @param weight a vector of weights 
#' @param alt alternative hypothesis for p-value calculation
#' @keywords internal
#' @useDynLib CaDrA ks_genescore_wrap_
#' 
#' @return need a return value here
ks_genescore_wrap <- function(n_x, y, weight, alt="less") {
  
  if(length(alt) > 0){
    alt_int<- switch(alt, two.sided=0L, less=1L, greater=-1L, 1L)
  } else {
    alt_int <- 1L
  }
  
  # Ensure the right type of input
  y <- as.integer(y)
  n_x <- as.integer(n_x)
  if(length(weight) > 1) weight <- as.numeric(weight)
  res <- .Call(ks_genescore_wrap_, n_x, y, weight, alt_int)
  res
  
}



#' Compute KS scores for each row of a given matrix
#'
#' Compute directional Kolmogorov-Smirnov scores for each row of a 
#' given binary matrix
#' @param mat matrix of binary features to compute row-wise ks scores for
#' @param alt an integer value specifying the alternative hypothesis
#' @param weight a vector of weights to use if performing a weighted-KS test
#' @keywords internal
#' @useDynLib CaDrA ks_genescore_mat_
#' 
#' @return Two lists: score and p-value
ks_genescore_mat <- function(mat, alt="less", weight) {
  
  if(!is.matrix(mat)) 
    stop("Input argument to ks_genescore_mat function is not a matrix")
  
  if(length(alt) > 0){
    alt_int<- switch(alt, two.sided=0L, less=1L, greater=-1L, 1L)
  } else {
    alt_int <- 1L
  }
  
  # Ensure the right type of input
  mat.num <- matrix(as.numeric(mat), ncol=ncol(mat), nrow=nrow(mat))
  weight <- if( length(weight) > 1 ) as.numeric(weight)
  res <- .Call(ks_genescore_mat_, mat.num, weight, alt_int)
  res
  
}


