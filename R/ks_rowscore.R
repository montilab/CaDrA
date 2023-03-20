#'
#' Kolmogorov-Smirnov Scoring Method
#'
#' Compute directional KS scores for each row of a given binary feature matrix
#'
#' @param FS_mat a matrix of binary features where 
#' rows represent features of interest (e.g. genes, transcripts, exons, etc...)
#' and columns represent the samples.
#' @param input_score a vector of continuous scores representing a phenotypic
#' readout of interest such as protein expression, pathway activity, etc.
#' The \code{input_score} object must have names or labels that match the column
#' names of FS_mat object.
#' @param weight a vector of weights to perform a \code{weighted-KS} test.
#' Default is \code{NULL}. If not NULL, \code{weight} must have labels or names
#' that match labels of \code{input_score}.
#' @param alternative a character string specifies an alternative hypothesis
#' testing (\code{"two.sided"} or \code{"greater"} or \code{"less"}).
#' Default is \code{less} for left-skewed significance testing.
#' @param metric a character string specifies a metric to search for 
#' best features. \code{"pval"} or \code{"stat"} may be used which is 
#' corresponding to p-value or score statistic. Default is \code{pval}. 
#' 
#' @noRd
#' @useDynLib CaDrA ks_genescore_mat_
#' 
#' @examples 
#' 
#' mat <- matrix(c(1,0,1,0,0,0,0,0,1,0, 
#'                 0,0,1,0,1,0,1,0,0,0,
#'                 0,0,0,0,1,0,1,0,1,0), nrow=3)
#'
#' colnames(mat) <- 1:10
#' row.names(mat) <- c("TP_1", "TP_2", "TP_3")
#'
#' set.seed(42)
#' input_score = rnorm(n = ncol(mat))
#' names(input_score) <- colnames(mat)
#' 
#' ks_rs <- ks_rowscore(
#'    FS_mat = mat,
#'    input_score = input_score,
#'    weight = NULL,
#'    alternative = "less",
#'    metric = "pval"
#' )
#'
#' @return return a vector of row-wise scores where its labels or names 
#' must match the row names of \code{FS_mat} object
#' 
ks_rowscore <- function
(
  FS_mat,
  input_score,
  weight = NULL,
  alternative = c("less", "greater", "two.sided"),
  metric = c("stat", "pval")
)
{

  metric <- match.arg(metric)
  alternative <- match.arg(alternative)
  
  # KS is a ranked-based method
  # So we need to sort input_score from highest to lowest values
  input_score <- sort(input_score, decreasing=TRUE)
  
  # Re-order the matrix based on the order of input_score
  FS_mat <- FS_mat[, names(input_score), drop=FALSE]  
  
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
  ks <- .Call(ks_genescore_mat_, FS_mat, weight, alt_int)

  # Obtain score statistics and p-values from KS method
  stat <- ks[1,]
  pval <- ks[2,]
  
  # Compute the scores according to the provided metric
  scores <- ifelse(rep(metric, nrow(FS_mat)) %in% "pval", -log(pval), stat)
  names(scores) <- rownames(FS_mat)

  return(scores)
  
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
#' 
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



#' ks_genescore wrapper
#'
#' Compute directional Kolmogorov-Smirnov scores for each row of a given vector
#' @param n_x length of ranked list
#' @param y positions of geneset items in ranked list (ranks)
#' @param weight a vector of weights
#' @param alt alternative hypothesis for p-value calculation
#' (\code{"two.sided"} or \code{"greater"} or \code{"less"}).
#' Default is \code{less} for left-skewed significance testing.
#' 
#' @noRd
#' @useDynLib CaDrA ks_genescore_wrap_
#'
#' @return a numeric vector of lenght 2 with 2 values: score and p-value
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
