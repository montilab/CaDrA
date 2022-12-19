

#' Wilcoxon Rank Sum Scoring Method
#'
#' Compute directional Wilcoxon rank sum score for each row of a 
#' given binary matrix
#'
#' @param mat a matrix of binary features (required).
#' @param ranks a vector of sample rankings use to perform the 
#' \code{Wilcoxon test}. Default is \code{NULL}. If NULL, 
#' then samples are assumed to be ordered by increasing rankings. 
#' If not NULL, ranks must include labels or names that associated 
#' with the colnames of the feature matrix.
#' @param alternative a character string specifies an alternative 
#' hypothesis testing (\code{"two.sided"} or \code{"greater"} or 
#' \code{"less"}). Default is \code{less} for left-skewed significance testing. 
#'
#' @return A data frame with two columns: \code{score} and \code{p_value}
#' @examples
#' 
#' # Load R library
#' library(Biobase)
#
#' # Load pre-computed expression set
#' data(sim.ES)
#' 
#' # Define additional parameters and run the function
#' result <- wilcox_rowscore(
#'   mat = exprs(sim.ES), 
#'   ranks = NULL,
#'   alternative = "less"
#' )
#' 
#' @export
#' @importFrom purrr map_dfr
wilcox_rowscore <- function
(
  mat,
  ranks = NULL,
  alternative = c("less", "greater", "two.sided")
)
{
  
  alternative <- match.arg(alternative)
  
  verbose("Using Wilcoxon method for features scoring")
  
  ## Make sure mat variable is a matrix
  mat <- as.matrix(mat)
  
  # If mat has only one column, it must be converted to a 
  # row-wise matrix form as it is needed for 
  # backward_forward_search() computation
  # mat must have rownames to track features and columns to track samples
  # for n = 1 case, it is only in backward_forward_search(), 
  # thus we can assign a random labels to it
  if(ncol(mat) == 1){
    mat <- matrix(t(mat), nrow=1, byrow=TRUE, 
                  dimnames = list("my_label", rownames(mat))) 
  }
  
  # Check if the matrix has only binary 0 or 1 values 
  if(length(mat) == 0 || !is.matrix(mat) || 
     any(!mat %in% c(0,1)) || any(is.na(mat)))
    stop("mat variable must be a matrix of binary values (no empty values).")
  
  # Make sure the mat variable has rownames for features tracking
  if(is.null(rownames(mat)))
    stop("The mat object does not have rownames or ",
         "featureData to track the features by. ",
         "Please provide unique features or rownames for the expression matrix")
  
  # If ranks for samples are not provided, assume it's ordered by 
  # decreasing ranks and assign ranks as 1:N (N: number of samples)
  if(length(ranks) > 0){
    if(length(ranks) != ncol(mat)){
      stop("'The provided ranks must have the same length as the number of ",
           "columns in the expression matrix.")
    }else{
      # check if ranks has any labels or names
      if(length(names(ranks)) == 0){
        stop("The ranks object must have names or labels that match ",
             "the colnames of the expression matrix.\n")
      }
      
      # check if ranks has labels or names that matches the 
      # colnames of the expression matrix
      if(any(!names(ranks) %in% colnames(mat))){
        stop("The ranks object have names or labels that do not match the ",
             "colnames of the expression matrix.\n")
      }
      
      # match colnames of expression matrix with names of provided ranks values
      # if nrow = 1, if it is, convert to matrix form as it is needed for 
      # backward_forward_search with one dimension matrix computation
      if(nrow(mat) == 1){
        mat <- matrix(t(mat[,names(ranks)]), 
                      nrow=1, byrow=TRUE, 
                      dimnames = list(rownames(mat), colnames(mat)))
      }else{
        mat <- mat[,names(ranks)]
      }
    }
  }else{
    ranks <- seq(1, ncol(mat))
  }
  
  # Check if the dataset has any all 0 or 1 features 
  # (these are to be removed since they are not informative)
  if(any(rowSums(mat) == 0) || any(rowSums(mat) == ncol(mat))){
    warning("Provided dataset has features that are either all 0 or 1. ",
            "These features will be removed from the computation.")
    mat <- mat[!(rowSums(mat) == 0 | rowSums(mat) == ncol(mat)),]
  }
  
  # Make sure matrix is not empty after removing uninformative features
  if(nrow(mat) == 0){
    stop("After removing features that are either all 0 or 1. ",
         "There are no more features remained for downsteam computation.\n")
  }
  
  #Compute the wilcox rank sum statitic and p-value per row in the matrix
  wilcox <- apply(X=mat, MARGIN=1, function(x, r=ranks){
    wilcox_score(x=r[which(x==1)], 
                     y=r[which(x==0)], alternative=alternative) 
  })
  
  # Convert list to data.frame
  wilcox <- data.frame(score=wilcox[1,], p_value=wilcox[2,])
  rownames(wilcox) <- rownames(mat)
  
  return(wilcox)
  
}



#' Compute rank sum scores for a given binary feature
#'
#' @param x an integer ranked values for group 1
#' @param y an integer ranked values for group 2
#' @param mu a number uses as an optional parameter to form a null hypothesis. 
#' Default is \code{0}.
#' @param alternative alternative hypothesis for p-value calculation 
#' (\code{"two.sided"} or \code{"greater"} or \code{"less"}). 
#' Default is \code{less} for left-skewed significance testing.
#' @param paired whether to perform paired test. Default is \code{FALSE}.
#' @param exact whether to compute exact p-value. Default is \code{FALSE}.
#' @param correct whether to consider continuity correction for p-value. 
#' Default is \code{TRUE}.
#' 
#' @examples 
#' 
#' # Load R library
#' library(Biobase)
#
#' # Load pre-computed expression set
#' data(sim.ES)
#' 
#' # set seed
#' set.seed(123)
#' 
#' # Provide a vector of continuous scores for a target profile 
#' # The scores must have labels or names 
#' # that match the colnames of expression matrix
#' input_score = rnorm(n = ncol(sim.ES))
#' names(input_score) <- colnames(sim.ES)
#' 
#' # Sort input_score from highest to lowest values
#' input_score <- sort(input_score, decreasing=TRUE)
#' 
#' # Extract the expression matrix
#' mat <- exprs(sim.ES)
#' 
#' # Re-order the samples by input_score sorted from highest to lowest values
#' mat <- mat[,names(input_score)]
#' 
#' # Define ranks of input_score
#' ranks <- seq_along(input_score)
#' 
#' # Define alternative
#' alternative <- "less"
#' 
#' # Compute the wilcox rank sum statitic and p-value per row in the matrix
#' wilcox <- apply(X=mat, MARGIN=1, function(x, r=ranks){
#'  wilcox_score(x=r[which(x==1)], y=r[which(x==0)], 
#'  alternative=alternative) 
#' })
#'
#' @return a data frame with two columns: \code{score} and \code{p_value}
#' @export
#'
#' @importFrom stats pnorm pwilcox
wilcox_score <- function 
(
  x,                                                  
  y,                                                  
  mu = 0,                                             
  alternative = c("less", "greater", "two.sided"),                          
  paired = FALSE,                                     
  exact = FALSE,                                      
  correct = TRUE
) 
{
  
  alternative <- match.arg(alternative)
  
  if (!missing(mu) && ((length(mu) > 1L) || !is.finite(mu))) 
    stop("'mu' must be a single number")
  if (!is.numeric(x)) 
    stop("'x' must be numeric")
  if (!is.numeric(y)) 
    stop("'y' must be numeric")
  
  DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  
  if (length(x) < 1L) {
    #print(x)
    stop("not enough (finite) 'x' observations")}
  
  if (length(y) < 1L){
    #print(y)
    stop("not enough 'y' observations")}
  
  METHOD <- "Wilcoxon rank sum test"
  
  ##### Modification ######
  # Take input as ranks instead of continuous measures 
  # (normally internally ranked: see below)
  r <- c(x,y)
  
  #r <- rank(c(x - mu, y))
  
  n.x <- as.double(length(x))
  n.y <- as.double(length(y))
  
  if (is.null(exact)) 
    exact <- (n.x < 50) && (n.y < 50)
  
  STATISTIC <- c(W = sum(r[seq_along(x)]) - n.x * (n.x + 1)/2)
  TIES <- (length(r) != length(unique(r)))
  
  if (exact && !TIES) {
    
    PVAL <- switch(alternative, two.sided = {
      p <- if (STATISTIC > (n.x * n.y/2)) 
        pwilcox(STATISTIC - 1, n.x, n.y, 
                lower.tail = FALSE) else pwilcox(STATISTIC, n.x, n.y)
      min(2 * p, 1)
    }, greater = {
      pwilcox(STATISTIC - 1, n.x, n.y, lower.tail = FALSE)
    }, less = pwilcox(STATISTIC, n.x, n.y))
    
  } else {
    
    NTIES <- table(r)
    z <- STATISTIC - n.x * n.y/2
    SIGMA <- sqrt((n.x * n.y/12) * ((n.x + n.y + 1) - 
                                      sum(NTIES^3 - NTIES)/
                                      ((n.x + n.y) * (n.x + n.y - 1))))
    
    if (correct) {
      CORRECTION <- switch(alternative, two.sided = sign(z) * 
                             0.5, greater = 0.5, less = -0.5)
      METHOD <- paste(METHOD, "with continuity correction")
    }
    
    z <- (z - CORRECTION)/SIGMA
    
    PVAL <- switch(alternative, less = pnorm(z), 
                   greater = pnorm(z, lower.tail = FALSE), 
                   two.sided = 2 * min(pnorm(z), 
                                       pnorm(z, lower.tail = FALSE)))
    
    if (exact && TIES) 
      warning("cannot compute exact p-value with ties")
    
  }
  
  names(mu) <- ifelse (paired || !is.null(y), "location shift", "location") 
  
  RVAL <- list(statistic = STATISTIC, 
               parameter = NULL, p.value = as.numeric(PVAL), 
               null.value = mu, alternative = alternative, method = METHOD, 
               data.name = DNAME)
  
  class(RVAL) <- "htest"
  
  return(c(score=RVAL$statistic, p_value=RVAL$p.value))         
  
}

