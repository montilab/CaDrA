

#' Row-wise matrix Kolmogorov-Smirnov scoring
#' 
#' Compute directional KS scores for each row of a given binary matrix
#' @param mat matrix of binary features to compute row-wise scores for based on the  Kolmogorov-Smirnov test
#' @param weight a vector of weights to use if performing a weighted-KS test. Default is NULL. Value passed to ks.genescore() function 
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided","less" or "greater". Value passed to ks.genescore() function
#' @export
custom_genescore_mat <- function
(
  mat, 
  weight = NULL,
  alternative = "less"                          
)
{
  
  # ## Required R packages
  # library(CaDrA)
  # library(Biobase)
  # library(tidyverse)
  # 
  # ## Source the ks genescore method
  # source("/Users/reinachau/Documents/CaDrA/CaDrA/R/ks_genescore.R")
  # 
  # ## Read in the simulated dataset from CaDrA
  # data("sim.ES")
  # 
  # ## Define the parameters for testing
  # mat <- exprs(sim.ES)
  # weight = NULL                                 # weights for weighted score (see Subramanian et al.) (usually sort(score))
  # alternative = "less"                          # alternative hypothesis for p-value calculation
  
  # Check if the ES is provided
  if(length(mat) == 0 || !is.matrix(mat) || any(!mat %in% c(0,1)))
    stop("mat variable must be a matrix with binary values only.\n")
    
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
  
  # Check if weight variable is provided
  if(length(weight) > 0){
    
    if(length(weight) < ncol(mat))
      stop("'The provided weight must have the same length as the number of samples in feature matrix.\n")
    
    print("Using weighted method for KS testing with the provided weight..\n")
    
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

  #Compute the ks statitic and p-value per row in the matrix
  ks <- 1:nrow(mat) %>% 
    map_dfr(
      function(r){
        #r=1;
        x=mat[r,]; y=which(x==1)
        ks_genescore(x=x, y=y, alternative=alternative, weight=weight) 
      }
    )

  return(ks)
  
}

