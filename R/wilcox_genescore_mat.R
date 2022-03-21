

#' Wilcoxon Rank Sum Scoring Method
#'
#' Compute directional Wilcoxon rank sum score for each row of a given binary matrix
#'
#' @param mat a matrix of binary features (required).
#' @param ranks a vector of sample rankings use to perform the \code{Wilcoxon test}. Default is \code{NULL}. If NULL, then samples are assumed to be ordered by increasing rankings. If not NULL, ranks must include labels or names that associated with the colnames of the feature matrix.
#' @param alternative a character string specifies an alternative hypothesis testing (\code{"two.sided"} or \code{"greater"} or \code{"less"}). Default is \code{less} for left-skewed significance testing. 
#' @param verbose a logical value indicates whether or not to print the diagnostic messages. Default is \code{FALSE}. 
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
#' wilcox_genescore_mat_result <- wilcox_genescore_mat(
#'   mat = exprs(sim.ES), 
#'   ranks = NULL,
#'   alternative = "less"
#' )
#' 
#' @export
#' @importFrom purrr map_dfr
wilcox_genescore_mat <- function
(
  mat,
  ranks = NULL,
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
  
  # If ranks for samples are not provided, assume it's ordered by decreasing ranks and assign ranks as 1:N (N: number of samples)
  if(length(ranks) > 0){
    if(length(ranks) != ncol(mat)){
      stop("'The provided ranks must have the same length as the number of columns in the expression matrix.\n")
    }else{
      # check if ranks has any labels or names
      if(length(names(ranks)) == 0){
        stop("The ranks object must have names or labels that match the colnames of the expression matrix.\n")
      }
      
      # check if ranks has labels or names that matches the colnames of the expression matrix
      if(any(!names(ranks) %in% colnames(mat))){
        stop("The ranks object have names or labels that do not match the colnames of the expression matrix.\n")
      }
      
      # match colnames of expression matrix with names of provided ranks values
      # if nrow = 1, if it is, convert to matrix form as it is needed for backward_forward_search with one dimension matrix computation
      if(nrow(mat) == 1){
        mat <- matrix(t(mat[,names(ranks)]), nrow=1, byrow=T, dimnames = list(rownames(mat), colnames(mat)))
      }else{
        mat <- mat[,names(ranks)]
      }
    }
  }else{
    ranks <- seq(1, ncol(mat))
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
    alternative <- "less"
    warning("No alternative hypothesis specified. Using 'less' by default.\n")
  }else if(length(alternative) == 1 && !alternative %in% c("two.sided", "greater", "less")){
    stop(paste0(alternative, collapse=", "), " is not a valid alternative hypothesis. Alternative hypothesis must be 'two.sided', 'greater', or 'less'.\n")
  }else if(length(alternative) > 1 && all(!alternative %in% c("two.sided", "greater", "less"))){
    stop(paste0(alternative, collapse=", "), " is not a valid alternative hypothesis. Alternative hypothesis must be 'two.sided', 'greater', or 'less'.\n")
  }else if(length(alternative) > 1 && any(alternative %in% c("two.sided", "greater", "less"))){
    alternative <- alternative[which(alternative %in% c("two.sided", "greater", "less"))][1]
    warning("More than one alternative hypothesis were specified. Only the first valid alternative hypothesis, '", alternative, "', is used.\n")
  }
  
  #Compute the wilcox rank sum statitic and p-value per row in the matrix
  wilcox <- 1:nrow(mat) %>% 
    purrr::map_dfr(
      function(r){
        #r=1;
        feature = mat[r,]; x = ranks[which(feature==1)]; y = ranks[which(feature==0)]    
        wilcox_genescore(x = x, y = y, alternative=alternative) 
      }
    )
  
  colnames(wilcox) <- c("score", "p_value")
  rownames(wilcox) <- rownames(mat)
  
  return(wilcox)
  
}



