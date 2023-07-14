
#' Calculate row-wise scores of a given binary feature set based on
#' a given scoring method
#'
#' @param FS a matrix of binary features or a SummarizedExperiment class object 
#' from SummarizedExperiment package where rows represent features of interest 
#' (e.g. genes, transcripts, exons, etc...) and columns represent the samples. 
#' The assay of FS contains binary (1/0) values indicating the presence/absence 
#' of omics features.
#' @param input_score a vector of continuous scores representing a phenotypic
#' readout of interest such as protein expression, pathway activity, etc.
#' 
#' NOTE: \code{input_score} object must have names or labels that match the 
#' column names of \code{FS} object.
#' @param seed_names a vector of one or more features representing known 
#' “causes” of activation or features associated with a response of interest.
#' It is applied for \code{method = "revealer"} only.
#' @param method a character string specifies a scoring method that is
#' used in the search. There are 6 options: (\code{"ks_pval"} or \code{ks_score}
#' or \code{"wilcox_pval"} or \code{wilcox_score} or 
#' \code{"revealer"} (conditional mutual information from REVEALER) or
#' \code{"custom"} (a customized scoring method)). 
#' Default is \code{ks_pval}.
#' @param custom_function if method is \code{"custom"}, specifies
#' the name of the customized function here. Default is \code{NULL}.
#' 
#' NOTE: \code{custom_function} must take \code{FS} and \code{input_score} 
#' as its input arguments, and its final result must return a vector of row-wise 
#' scores ordered from most significant to least significant where its labels or 
#' names matched the row names of \code{FS} object.
#' @param custom_parameters if method is \code{"custom"}, specifies a list of
#' additional arguments (excluding \code{FS} and \code{input_score}) to be 
#' passed to \code{custom_function}. Default is \code{NULL}.
#' @param alternative a character string specifies an alternative hypothesis
#' testing (\code{"two.sided"} or \code{"greater"} or \code{"less"}).
#' Default is \code{less} for left-skewed significance testing.
#' 
#' NOTE: This argument is applied to KS and Wilcoxon method
#' @param weights if method is \code{ks_score} or \code{ks_pval}, specifying a 
#' vector of weights will perform a weighted-KS testing. Default is \code{NULL}.
#' @param do_check a logical value indicates whether or not to validate if the  
#' given parameters (\code{FS} and \code{input_score}) are valid inputs. 
#' Default is \code{TRUE}.
#' @param verbose a logical value indicates whether or not to print the
#' diagnostic messages. Default is \code{FALSE}.
#' @param ... additional parameters to be passed to \code{custom_function}
#' 
#' @return return a vector of row-wise positive scores where it is ordered from 
#' most significant to least significant (e.g. from highest to lowest values) 
#' and its labels or names must match the row names of \code{FS} object
#' 
#' @examples
#' 
#' # Create a feature matrix
#' mat <- matrix(c(1,0,1,0,0,0,0,0,1,0, 
#'                 0,0,1,0,1,0,1,0,0,0,
#'                 0,0,0,0,1,0,1,0,1,0), nrow=3)
#' 
#' colnames(mat) <- 1:10
#' row.names(mat) <- c("TP_1", "TP_2", "TP_3")
#' 
#' # Create a vector of observed input scores
#' set.seed(42)
#' input_score = rnorm(n = ncol(mat))
#' names(input_score) <- colnames(mat)
#'
#' # Run the ks method
#' ks_rowscore_result <- calc_rowscore(
#'   FS = mat,
#'   input_score = input_score,
#'   seed_names = NULL,
#'   method = "ks_pval",
#'   weights = NULL,
#'   alternative = "less"
#' )
#'
#' # Run the wilcoxon method
#' wilcox_rowscore_result <- calc_rowscore(
#'   FS = mat,
#'   input_score = input_score,
#'   seed_names = NULL,
#'   method = "wilcox_pval",
#'   alternative = "less"
#' )
#'
#' # Run the revealer method
#' revealer_rowscore_result <- calc_rowscore(
#'   FS = mat,
#'   input_score = input_score,
#'   method = "revealer",
#'   seed_names = NULL
#' )
#' 
#' # A customized function using ks-test function
#' customized_rowscore <- function(FS, input_score, seed_names = NULL, alternative="less"){
#'   
#'   ks <- apply(FS, 1, function(r){ 
#'     x = input_score[which(r==1)]; 
#'     y = input_score[which(r==0)];
#'     res <- ks.test(x, y, alternative=alternative)
#'     return(c(res$statistic, res$p.value))
#'   })
#'   
#'   # Obtain score statistics and p-values from KS method
#'   stat <- ks[1,]
#'   
#'   # Change values of 0 to the machine lowest value to avoid taking -log(0)
#'   pval <- ks[2,]
#'   pval[which(pval == 0)] <- .Machine$double.xmin
#'   
#'   # Compute the -log scores for pval
#'   # Make sure scores has names that match the row names of FS object
#'   scores <- -log(pval)
#'   names(scores) <- rownames(FS)
#'   
#'   return(scores)
#'   
#' }
#' 
#' # Search for best features using a custom-defined function
#' custom_rowscore_result <- calc_rowscore(
#'   FS = mat,
#'   input_score = input_score,
#'   seed_names = NULL,
#'   method = "custom",
#'   custom_function = customized_rowscore,            
#'   custom_parameters = NULL  
#' )
#'
#' @export
#' @import SummarizedExperiment
calc_rowscore <- function(
    FS,
    input_score,
    seed_names = NULL,
    method = c("ks_pval", "ks_score", "wilcox_pval", "wilcox_score", 
               "revealer", "custom"),
    custom_function = NULL,
    custom_parameters = NULL,   
    alternative = c("less", "greater", "two.sided"),
    weights = NULL,
    do_check = TRUE,
    verbose = FALSE,
    ...
){

  # Set up verbose option
  options(verbose = verbose)
  
  # Match arguments
  method <- match.arg(method)
  alternative <- match.arg(alternative)
  
  # Check if FS is a matrix or a SummarizedExperiment class object
  if(!is(FS, "SummarizedExperiment") && !is(FS, "matrix"))
    stop("'FS' must be a matrix or a SummarizedExperiment class object
         from SummarizedExperiment package")
  
  # Retrieve the feature binary matrix
  if(is(FS, "SummarizedExperiment")){
    FS_mat <- SummarizedExperiment::assay(FS)
  }else{
    FS_mat <- FS
  }
  
  # Check if FS and input_score are valid inputs
  if(do_check == TRUE) 
    check_data_input(FS_mat = FS_mat, input_score = input_score, 
                     seed_names = seed_names, do_check = do_check)
  
  # Define metric value based on a given scoring method
  if(length(grep("score", method)) > 0){
    metric <- "stat"
  }else{
    metric <- "pval"   
  }
  
  # Extract only the method value (no metric info)
  # based on a given method string
  method <- gsub("_score|_pval", "", method)
  
  # Select the appropriate method to compute row-wise directional scores
  rscores <- switch(
    method,
    ks = ks_rowscore(
      FS = FS_mat,
      input_score = input_score,
      seed_names = seed_names,
      weights = weights,
      alternative = alternative,
      metric = metric
    ),
    wilcox = wilcox_rowscore(
      FS = FS_mat,
      input_score = input_score,
      seed_names = seed_names,
      alternative = alternative,
      metric = metric
    ),
    revealer = revealer_rowscore(
      FS = FS_mat,
      input_score = input_score,
      seed_names = seed_names,
      assoc_metric = "IC"
    ),
    custom = custom_rowscore(
      FS = FS,
      input_score = input_score,
      seed_names = seed_names,
      custom_function = custom_function,
      custom_parameters = custom_parameters,
      method = method, 
      alternative = alternative, 
      weights = weights, 
      do_check = do_check, 
      verbose = verbose
    )
  )
  
  # Re-order FS in a decreasing order (from most to least significant)
  # This comes in handy when doing the top-N evaluation of
  # the top N 'best' features
  rscores <- rscores[order(rscores, decreasing=TRUE)]

  return(rscores)

}



