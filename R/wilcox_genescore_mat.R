
#' Row-wise matrix Wilcoxon rank sum scoring
#'
#' @param mat matrix of binary features to compute row-wise scores for based on the Wilcoxon rank sum test
#' @param rank a vector of ranks to use when performing the Wilcoxon test. Default is NULL. If NULL, then samples are assumed to be ordered by increasing ranking. Value passed to wilcox_genescore() function 
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided","less" or "greater". Value passed to wilcox_genescore() function
#'
#' @return A data frame
#' @export
#' @importFrom purrr map_dfr
wilcox_genescore_mat <- function
(
  mat,
  rank = NULL,
  alternative = "less"
)
{
  
  # If no alternative is specified, we use "less" as default.
  if(length(alternative) == 0 || nchar(alternative) == 0){
    warning("No alternative hypothesis specified. Using 'less' by default ..\n")
    alternative <- "less"
  }else if(length(alternative) == 1 && !alternative %in% c("two.sided", "greater", "less")){
    warning(paste0(alternative, collapse=", "), " is not a valid alternative hypothesis. Alternative hypothesis specified must be 'two.sided', 'greater', or 'less'. Using 'less' by default.\n")
    alternative <- "less"    
  }else if(length(alternative) > 1 && all(!alternative %in% c("two.sided", "greater", "less"))){
    warning(paste0(alternative, collapse=", "), " is not a valid alternative hypothesis. Alternative hypothesis specified must be 'two.sided', 'greater', or 'less'. Using 'less' by default.\n")
    alternative <- "less"
  }else if(length(alternative) > 1 && any(alternative %in% c("two.sided", "greater", "less"))){
    alternative <- alternative[which(alternative %in% c("two.sided", "greater", "less"))][1]
    warning("More than one alternative hypothesis were specified. Only the first valid alternative hypothesis, '", alternative, "', is used.\n")
  }
  
  # If ranks for samples are not provided, assume it's ordered by decreasing rank and assign rank 1:N (N: number of samples)
  if(length(rank) > 0){
    if(length(rank) < ncol(mat))
      stop("'The provided rank must have the same length as the number of samples in the feature matrix.\n")
    verbose("Using provided rank for Wilcoxon rank sum testing.\n")
  }else{
    rank <- seq(1, ncol(mat))
  }
  
  ###### REMOVE INVALID FEATURE #####
  #####################################
  
  # Check if the dataset has any all 0 or 1 features (these are to be removed since they are not informative)
  if(any(rowSums(mat) == 0) || any(rowSums(mat) == ncol(mat))){
    warning("Provided dataset has features that are either all 0 or 1. These features will be removed from the computation.\n")
    mat <- mat[!(rowSums(mat) == 0 | rowSums(mat) == ncol(mat)),]
  }
  
  if(nrow(mat) < 2)
    warning("You are computing a row-wise statistic over a matrix with nrow < 2.\n")
  
  #Compute the wilcox rank sum statitic and p-value per row in the matrix
  wilcox <- 1:nrow(mat) %>% 
    purrr::map_dfr(
      function(r){
        #r=1;
        feature = mat[r,]; x = rank[which(feature==1)]; y = rank[which(feature==0)]    
        wilcox_genescore(x=x, y=y, alternative=alternative)
      }
    )
  
  return(wilcox)
  
}
