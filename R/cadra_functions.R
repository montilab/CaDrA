
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

#' cadra_check_input
#' 
#' Checks if the input values to CaDrA function are valid.
#'
#' @param ES an expression set of binary features (required). It must be a 
#' \code{BioBase expressionSet} object. The rownames of the expression set must 
#' contain unique features which are used to search for best features.   
#' @param input_score a vector of continuous values of a response of 
#' interest (required). The \code{input_score} object must have names or 
#' labels that match the colnames of the expression matrix.
#' @param n_perm an integer specifies the number of permutations to perform. 
#' @param ncores an integer specifies the number of CPU cores.
#' @noRd
cadra_check_input <- function(
    ES,
    input_score,
    method,
    n_perm,
    ncores
){
  
  # Check if the ES is provided and is a BioBase ExpressionSet object
  stopifnot("'ES' must be an ExpressionSet class argument"=
              (length(ES) > 0 && is(ES, "ExpressionSet") ) )
  
  # Check if the dataset has only binary 0 or 1 values 
  stopifnot("The expression matrix (ES) must contain only binary values"=
              (all(exprs(ES) %in% c(0,1))) )
  
  
  # Make sure the input ES has rownames for features tracking
  stopifnot("The expression matrix (ES) must have rownames or featureData"=
              (!is.null(rownames(ES))) )
  
  
  # Check input_score is provided and is a continuous values with no NAs
  stopifnot("invalid input_score"= (length(input_score) != 0 &&
                                      all(is.numeric(input_score)) &&
                                      sum(is.na(input_score))==0 &&
                                      !is.null(names(input_score))) )
  
  
  stopifnot("invalid number of permutations (nperm)"=
              (length(n_perm)==1 && !is.na(n_perm) &&  
                 is.numeric(n_perm) && n_perm > 0) )
  
  stopifnot("invalid number of CPU cores (ncores)"=
              (length(ncores)==1 && !is.na(ncores) &&  
                 is.numeric(ncores) && ncores > 0) )
  
}


#' cadra_plot
#' 
#' Plots result of CaDrA function
#'
#' @param top_N an integer specifies the number of features to start the 
#' search over, starting from the top 'N' features in each case. If \code{top_N} 
#' is provided, then \code{search_start} parameter will be ignored. Default is 
#' \code{1}.
#' @param search_start a list of character strings (separated by commas) 
#' which specifies feature names within the expression set object to start 
#' the search with. If \code{search_start} is provided, then \code{top_N}
#' parameter will be ignored. Default is \code{NULL}.
#' @param obs_best_score a numeric value corresponding to the best observed 
#' score or p-value and later use to compare against permutation based scores or
#' p-values. Default is \code{NULL}. If set to NULL, we will compute the observed 
#' best score based on the given \code{input_score} and \code{ES} variables.
#' @param perm_pval pval returned by CaDrA calculations
#' @param perm_best_scores permutations best scores
#' @noRd
cadra_plot <- function(top_N, search_start, obs_best_score,
                       perm_pval, perm_best_scores)
{
  
  plot_title <- paste0("Emperical Null distribution (N = ", 
                       length(perm_best_scores), ")\n Permutation p-val <= ", 
                       round(perm_pval, 5), "\nBest observed score: ", 
                       round(obs_best_score, 5))
  
  if(!is.null(top_N)){
    plot_title <- paste0(plot_title,"\n Top N: ", top_N)
  }else{
    plot_title <- paste0(plot_title,"\n Seed: ", search_start)
  }
  
  #Here, let us plot the absolute values of the permutation p-values, 
  # for simplicity
  # you only consider absolute values when calculating the permutation p-vals.
  g <- ggplot(data = data.frame("x" = perm_best_scores), aes(x = .data$x)) +
    geom_histogram(fill = "black", color = "gray") +
    theme_classic() +
    theme(
      axis.line.x=element_line(color = "black"),
      axis.line.y=element_line(color = "black")
    )
  
  g <- g + geom_vline(xintercept = obs_best_score, 
                      linetype = "longdash", size = 1.5, colour = "red") +
    labs(
      title = plot_title,
      x = "Score",
      y = "Count"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0))
  
  g
}


#' ks_test_d_wrap_ wrapper
#'
#' Compute directional Kolmogorov-Smirnov scores
#' @param n_x length of ranked list
#' @param y positions of geneset items in ranked list (ranks)
#' @param alt alternative hypothesis for p-value calculation 
#' (\code{"two.sided"} or \code{"greater"} or \code{"less"}). 
#' Default is \code{less} for left-skewed significance testing.
#' @noRd
#' @useDynLib CaDrA ks_test_d_wrap_
#' 
#' @return need a return value here
ks_test_double_wrap <- function(n_x, y, alt=c("less", "greater", "two.sided")) {
  
  
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
#' (\code{"two.sided"} or \code{"greater"} or \code{"less"}). 
#' Default is \code{less} for left-skewed significance testing.
#' @noRd
#' @useDynLib CaDrA ks_plot_wrap_
#' 
#' @return need a return value here
ks_plot_wrap <- function(n_x, y, weight, alt=c("less", "greater", "two.sided")){
  
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
#' (\code{"two.sided"} or \code{"greater"} or \code{"less"}). 
#' Default is \code{less} for left-skewed significance testing.
#' @noRd
#' @useDynLib CaDrA ks_genescore_wrap_
#' 
#' @return need a return value here
ks_genescore_wrap <- function(n_x, y, weight, 
                              alt=c("less", "greater", "two.sided")) {
  
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
#' (\code{"two.sided"} or \code{"greater"} or \code{"less"}). 
#' Default is \code{less} for left-skewed significance testing.
#' @param weight a vector of weights to use if performing a weighted-KS test
#' @noRd
#' @useDynLib CaDrA ks_genescore_mat_
#' 
#' @return Two lists: score and p-value
ks_gene_score_mat <- function(
  mat, 
  alt=c("less", "greater", "two.sided"), 
  weight
) {
  
  if(!is.matrix(mat)) 
    stop("Input argument to ks_gene_score_mat function is not a matrix")
  
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


