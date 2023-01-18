#'
#' Kolmogorov-Smirnov Scoring Method
#'
#' Compute directional KS scores for each row of a given binary feature matrix
#'
#' @param FS a feature set of binary features. It can be a matrix or
#' a \code{SummarizedExperiment} class object from SummarizedExperiment package.
#' @param input_score a vector of continuous scores representing a phenotypic
#' readout of interest such as protein expression, pathway activity, etc.
#' The \code{input_score} object must have names or labels that match the column
#' names of FS object.
#' @param weight a vector of weights to perform a \code{weighted-KS} test.
#' Default is \code{NULL}. If not NULL, \code{weight} must have labels or names
#' that match labels of \code{input_score}.
#' @param alternative a character string specifies an alternative hypothesis
#' testing (\code{"two.sided"} or \code{"greater"} or \code{"less"}).
#' Default is \code{less} for left-skewed significance testing.
#' @param warning a logical value indicates whether or not to print the 
#' diagnostic messages. Default is \code{TRUE}
#' 
#' @noRd
#' @useDynLib CaDrA ks_genescore_mat_
#'
#' @return A matrix with two columns: \code{score} and \code{p_value}
#' @import SummarizedExperiment
ks_rowscore <- function
(
  FS,
  input_score,
  weight = NULL,
  alternative = c("less", "greater", "two.sided"),
  warning = TRUE
)
{

  alternative <- match.arg(alternative)
  
  # Check of FS and input_score are valid inputs
  if(warning == TRUE) check_data_input(FS = FS, input_score = input_score, warning=warning)
  
  # Get the feature names
  feature_names <- rownames(FS)
  
  # KS is a ranked-based method
  # So we need to sort input_score from highest to lowest values
  input_score <- sort(input_score, decreasing=TRUE)

  # Re-order the matrix based on the order of input_score
  FS <- FS[, names(input_score)]
  
  # Extract the feature binary matrix
  if(class(FS)[1] == "SummarizedExperiment"){
    mat <- as.matrix(SummarizedExperiment::assay(FS))
  }else if(class(FS)[1] == "matrix"){
    mat <- as.matrix(FS)
  }else{
    mat <- matrix(t(FS), nrow=1, byrow=TRUE,
                  dimnames=list(feature_names, names(FS)))
  }

  # Check if weight is provided
  if(length(weight) > 0){
    # Check if weight has any labels or names
    if(is.null(names(weight)))
      stop("The weight object must have names or labels that ",
           "match the labels of input_score\n")

    # Make sure its labels or names match the
    # the labels of input_score
    weight <- as.numeric(weight[names(input_score)])
  }

  # Get the alternative hypothesis testing method
  alt_int <- switch(alternative, two.sided=0L, less=1L, greater=-1L, 1L)

  # Compute the ks statistic and p-value per row in the matrix
  ks <- .Call(ks_genescore_mat_, mat, weight, alt_int)

  # Convert results as matrix
  # KS method returns both score and p-values
  ks_mat <- matrix(NA, nrow=nrow(mat), 
                   ncol=2, byrow=TRUE, 
                   dimnames=list(rownames(mat), c("score", "p_value")))
  
  ks_mat[,1] <- ks[1,]
  ks_mat[,2] <- ks[2,]
  
  return(ks_mat)

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
