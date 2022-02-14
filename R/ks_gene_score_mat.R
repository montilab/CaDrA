

#' Row-wise matrix Kolmogorov-Smirnov scoring
#'
#' Compute directional KS scores for each row of a given binary matrix
#' @param mat A matrix of binary features to compute row-wise scores based on the Kolmogorov-Smirnov test
#' @param weights a vector of weights to perform a weighted-KS test. Default is NULL. Value is passed to ks_gene_score() function.
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided","less" or "greater". Value passed to ks_gene_score() function.
#' @param verbose a logical indicating whether or not to verbose diagnostic messages. Default is TRUE. 
#'
#' @return A data frame with two columns: \code{score} and \code{p_value}
#' @export 
#' @importFrom purrr map_dfr
ks_gene_score_mat <- function
(
  mat, 
  weights = NULL,
  alternative = c("two.sided", "greater", "less"),
  verbose = TRUE
)
{
  
  # Set up verbose option
  options(verbose=verbose)
  
  # Check if the ES is provided
  if(length(mat) == 0 || !is.matrix(mat) || any(!mat %in% c(0,1)))
    stop("mat variable must be a matrix with binary values only.\n")
  
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
  
  # Check if weight variable is provided
  if(length(weights) > 0){
    if(length(weights) != ncol(mat))
      stop("'The provided weights must have the same length as the number of columns in the expression matrix.\n")
    verbose("Using weighted method for Kolmogorov-Smirnov testing with the provided weights...\n")
  }
  
  # Check if the dataset has any all 0 or 1 features (these are to be removed since they are not informative)
  if(any(rowSums(mat) == 0) || any(rowSums(mat) == ncol(mat))){
    warning("Provided dataset has features that are either all 0 or 1. These features will be removed from the computation.\n")
    mat <- mat[!(rowSums(mat) == 0 | rowSums(mat) == ncol(mat)),]
  }
  
  if(nrow(mat) < 2)
    warning("You are computing a row-wise statistic over a matrix with nrow < 2.\n")
  
  #Compute the ks statitic and p-value per row in the matrix
  ks <- 1:nrow(mat) %>% 
    purrr::map_dfr(
      function(r){
        #r=1;
        x=mat[r,]; y=which(x==1)
        ks_gene_score(x=x, y=y, alternative=alternative, weights=weights) 
      }
    )
  
  return(ks)
  
}

