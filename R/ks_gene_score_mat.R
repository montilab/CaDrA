
#' Kolmogorov-Smirnov Scoring Method
#'
#' Compute directional KS scores for each row of a given binary matrix
#' 
#' @param mat a matrix of binary features to compute row-wise scores based on the \code{Kolmogorov-Smirnov} test.
#' @param weights a vector of weights to perform a \code{weighted-KS} test. Default is \code{NULL}. If not NULL, weights must include labels or names that associated with the colnames of the feature matrix.
#' @param alternative a character string specifies an alternative hypothesis testing (\code{"two.sided"} or \code{"greater"} or \code{"less"}). Default is \code{less} for left-skewed significance testing. 
#' @param verbose a logical value indicates whether or not to print the diagnostic messages. Default is \code{FALSE}. 
#'
#' @return A data frame with two columns: \code{score} and \code{p_value}
#' @examples
#' 
#' # Load R library
#' library(Biobase)
#' 
#' # Load pre-computed expression set
#' data(sim.ES)
#' 
#' # Define additional parameters and start the function
#' ks_gene_score_mat_result <- ks_gene_score_mat(
#'   mat = exprs(sim.ES), 
#'   weights = NULL,
#'   alternative = "less"
#' )
#' 
#' @export 
#' @importFrom purrr map_dfr
ks_gene_score_mat <- function
(
  mat, 
  weights = NULL,
  alternative = c("two.sided", "greater", "less"),
  verbose = FALSE
)
{
  
  # Set up verbose option
  options(verbose=FALSE)
  
  # Check if the matrix has only binary 0 or 1 values 
  if(length(mat) == 0 || !is.matrix(mat) || any(!mat %in% c(0,1)) || any(is.na(mat)))
    stop("mat variable must be a matrix of binary values (no empty values).\n")
  
  # Make sure the mat variable has rownames for features tracking
  if(is.null(rownames(mat)))
    stop("The mat object does not have rownames or featureData to track the features by. Please provide unique features or rownames for the expression matrix.\n")
  
  # If weights for samples are not provided, assume it's ordered by decreasing weights and assign weights as 1:N (N: number of samples)
  if(length(weights) > 0){
    if(length(weights) != ncol(mat)){
      stop("'The provided weights must have the same length as the number of columns in the expression matrix.\n")
    }else{
      # check if weights has any labels or names
      if(length(names(weights)) == 0){
        stop("The weights object must have names or labels that match the colnames of the expression matrix.\n")
      }
      
      if(any(!names(weights) %in% colnames(mat))){
        stop("The weights object have names or labels that do not match the colnames of the expression matrix.\n")
      }
      
      # match colnames of expression matrix with names of provided weights values
      # iif nrow = 1, if it is, convert to matrix form as it is needed for backward_forward_search with one dimension matrix computation
      if(nrow(mat) == 1){
        mat <- matrix(t(mat[,names(weights)]), nrow=1, byrow=T, dimnames = list(rownames(mat), colnames(mat))) 
      }else{
        mat <- mat[,names(weights)]
      }
    }
  }
  
  # Check if the dataset has any all 0 or 1 features (these are to be removed since they are not informative)
  if(any(rowSums(mat) == 0) || any(rowSums(mat) == ncol(mat))){
    warning("Provided dataset has features that are either all 0 or 1. These features will be removed from the computation.\n")
    mat <- mat[!(rowSums(mat) == 0 | rowSums(mat) == ncol(mat)),]
  }
  
  # Make sure matrix is not empty after removing uninformative features
  if(nrow(mat) == 0){
    stop("After removing features that are either all 0 or 1. There are no more features remained for downsteam computation.\n")
  }
  
  # Give a warning if matrix has nrow < 2
  if(nrow(mat) < 2)
    warning("You are computing a row-wise statistic over a matrix with nrow < 2.\n")
  
  # If no alternative is specified, we use "less" as default.
  if(length(alternative) == 0 || nchar(alternative) == 0){
    warning("No alternative hypothesis specified. Using 'less' by default ..\n")
    alternative <- "less"
  }else if(length(alternative) == 1 && !alternative %in% c("two.sided", "greater", "less")){
    stop(paste0(alternative, collapse=", "), " is not a valid alternative hypothesis. Alternative hypothesis must be 'two.sided', 'greater', or 'less'.\n")
  }else if(length(alternative) > 1 && all(!alternative %in% c("two.sided", "greater", "less"))){
    stop(paste0(alternative, collapse=", "), " is not a valid alternative hypothesis. Alternative hypothesis must be 'two.sided', 'greater', or 'less'.\n")
  }else if(length(alternative) > 1 && any(alternative %in% c("two.sided", "greater", "less"))){
    alternative <- alternative[which(alternative %in% c("two.sided", "greater", "less"))][1]
    warning("More than one alternative hypothesis were specified. Only the first valid alternative hypothesis, '", alternative, "', is used.\n")
  } 
  
  #Compute the ks statitic and p-value per row in the matrix
  ks <- 1:nrow(mat) %>% 
    purrr::map_dfr(
      function(r){
        #r=1;
        x = mat[r,]; n.x = length(x); y = which(x==1);     
        ks_gene_score(n.x = n.x, y = y, weights = weights, alternative=alternative)
      }
    )
  
  colnames(ks) <- c("score", "p_value")
  rownames(ks) <- rownames(mat)
  
  return(ks)
  
}

