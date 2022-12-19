#' 
#' Kolmogorov-Smirnov Scoring Method
#' 
#' Compute directional KS scores for each row of a given binary matrix
#' 
#' @param mat a matrix of binary features to compute row-wise scores based on 
#' the \code{Kolmogorov-Smirnov} test.
#' @param weights a vector of weights to perform a \code{weighted-KS} test. 
#' Default is \code{NULL}. If not NULL, weights must include labels or names 
#' that associated with the colnames of the feature matrix.
#' @param alternative a character string specifies an alternative hypothesis 
#' testing (\code{"two.sided"} or \code{"greater"} or \code{"less"}). 
#' Default is \code{less} for left-skewed significance testing.
#' 
#' @return A data frame with two columns: \code{score} and \code{p_value}
#' @examples
#' 
#' #' Load R library
#' library(Biobase)
#' 
#' #' Load pre-computed expression set
#' data(sim.ES)
#' 
#' #' Define additional parameters and start the function
#' result <- ks_rowscore(
#'   mat = exprs(sim.ES),
#'   weights = NULL,
#'   alternative = "less"
#' )
#' 
#' @export
ks_rowscore <- function
(
  mat, 
  weights = NULL,
  alternative = c("less", "greater", "two.sided")
)
{
  
  alternative <- match.arg(alternative)
  
  verbose("Using Kolmogorov-Smirnov method for features scoring")
  
  ## Make sure mat variable is a matrix
  mat <- as.matrix(mat)
  
  # If mat has only one column, it must be converted to a row-wise matrix 
  # form as it is needed for backward_forward_search() computation
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
    stop("The mat object does not have rownames or featureData to ",
         "track the features by. Please provide unique features or rownames ",
         "for the expression matrix.\n")
  
  # If weights for samples are not provided, assume it's ordered by 
  # decreasing weights and assign weights as 1:N (N: number of samples)
  if(length(weights) > 0){
    if(length(weights) != ncol(mat)){
      stop("'The provided weights must have the same length as the number of ",
           "columns in the expression matrix.\n")
    }else{
      # check if weights has any labels or names
      if(length(names(weights)) == 0){
        stop("The weights object must have names or labels that ",
             "match the colnames of the expression matrix.\n")
      }
      
      # check if weights has labels or names that matches the 
      # colnames of the expression matrix
      if(any(!names(weights) %in% colnames(mat))){
        stop("The weights object have names or labels that do not match ",
             "the colnames of the expression matrix.\n")
      }
      
      # match colnames of expression matrix with 
      # names of provided weights values
      # iif nrow = 1, if it is, convert to matrix form as it is needed 
      # for backward_forward_search with one dimension matrix computation
      if(nrow(mat) == 1){
        mat <- matrix(t(mat[,names(weights)]), nrow=1, byrow=TRUE, 
                      dimnames = list(rownames(mat), colnames(mat))) 
      }else{
        mat <- mat[,names(weights)]
      }
    }
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
         "There are no more features remained for downsteam computation.")
  }
  
  # Compute the ks statistic and p-value per row in the matrix
  ks <- ks_rowscore_calc(mat=mat, alt=alternative, weight=weights)
  
  # Convert list to data.frame
  ks <- data.frame(score=ks[1,], p_value=ks[2,])
  rownames(ks) <- rownames(mat)
  
  return(ks)
  
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
ks_rowscore_calc <- function(
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
