
#' Calculate row-wise scores of a given binary feature set based on
#' a given scoring method
#'
#' @param FS_mat a matrix of binary features where 
#' rows represent features of interest (e.g. genes, transcripts, exons, etc...)
#' and columns represent the samples.
#' @param input_score a vector of continuous scores representing a phenotypic
#' readout of interest such as protein expression, pathway activity, etc.
#' The \code{input_score} object must have names or labels that match the column
#' names of FS object.
#' @param method a character string specifies a scoring method that is
#' used in the search. There are 4 options: (\code{"ks"} or \code{"wilcox"} or
#' \code{"revealer"} (conditional mutual information from REVEALER) or
#' or \code{custom} (a user-customized scoring method). 
#' Default is \code{ks}.
#' @param alternative a character string specifies an alternative hypothesis
#' testing (\code{"two.sided"} or \code{"greater"} or \code{"less"}).
#' Default is \code{less} for left-skewed significance testing.
#' @param weight if method is \code{ks}, specifies a vector of weights
#' to perform a weighted-KS testing. Default is \code{NULL}.
#' @param seed_names a vector of one or more features representing known “causes”
#' of activation or features associated with a response of interest.
#' It is used for \code{method = "revealer"} only.
#' @param custom_function if method is \code{"custom"}, just pass the name of the
#' customized function. Default is \code{NULL}.
#' @param custom_parameters if method is \code{"custom"}, specifies a list of
#' additional arguments (excluding \code{FS_mat} and \code{input_score}) to be passed
#' to \code{custom_function()}. Default is \code{NULL}.
#' @param do_check a logical value indicates whether or not to validate if the  
#' given parameters (FS and input_score) are valid inputs. 
#' Default is \code{TRUE}.
#' @param ... additional parameters to be passed to the custom_function() 
#' 
#' @return return a vector of scores ordered from most significant to least
#' significant where its labels or names must match the row names of FS_mat object
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
#'   FS_mat = mat,
#'   input_score = input_score,
#'   method = "ks_pval",
#'   weight = NULL,
#'   alternative = "less"
#' )
#'
#' # Run the wilcoxon method
#' wilcox_rowscore_result <- calc_rowscore(
#'   FS_mat = mat,
#'   input_score = input_score,
#'   method = "wilcox_pval",
#'   alternative = "less"
#' )
#'
#' # Run the revealer method
#' revealer_rowscore_result <- calc_rowscore(
#'   FS_mat = mat,
#'   input_score = input_score,
#'   method = "revealer",
#'   seed_names = NULL
#' )
#' 
#' # A customized function using ks-test function
#' customized_rowscore <- function(FS_mat, input_score, alternative="less"){
#'   
#'   ks <- apply(FS_mat, 1, function(r){ 
#'     x = input_score[which(r==1)]; 
#'     y = input_score[which(r==0)];
#'     res <- ks.test(x, y, alternative=alternative)
#'     return(c(res$statistic, res$p.value))
#'   })
#'   
#'   # Obtain score statistics and p-values from KS method
#'   stat <- ks[1,]
#'   pval <- ks[2,]
#'   
#'   # Compute the -log scores for pval
#'   # Make sure scores has names that match the row names of FS_mat object
#'   scores <- -log(pval)
#'   names(scores) <- rownames(FS_mat)
#'   
#'   # Remove scores that are Inf as it is resulted from
#'   # taking the -log(0). They are uninformative.
#'   scores <- scores[scores != Inf]  
#'   
#'   # Re-order FS_mat in a decreasing order (from most to least significant)
#'   # This comes in handy when doing the top-N evaluation of
#'   # the top N 'best' features
#'   scores <- scores[order(scores, decreasing=TRUE)]
#'   
#'   return(scores)
#'   
#' }
#' 
#' # Search for best features using a custom-defined function
#' custom_rowscore_result <- calc_rowscore(
#'   FS_mat = mat,
#'   input_score = input_score,
#'   method = "custom",
#'   custom_function = customized_rowscore,            
#'   custom_parameters = NULL  
#' )
#'
#' @export
calc_rowscore <- function(
    FS_mat,
    input_score,
    method = c("ks_pval", "ks_score", "wilcox_pval", "wilcox_score", "revealer", "custom"),
    alternative = c("less", "greater", "two.sided"),
    weight = NULL,
    seed_names = NULL,
    custom_function = NULL,
    custom_parameters = NULL,
    do_check = TRUE,
    ...
){

  # Match arguments
  method <- match.arg(method)
  alternative <- match.arg(alternative)
  
  # Check if FS and input_score are valid inputs
  if(do_check == TRUE) 
    check_data_input(FS = FS_mat, input_score = input_score, do_check=do_check)

  # Define metric value based on a given scoring method
  if(length(grep("score", method)) > 0){
    metric <- "stat"
  }else{
    metric <- "pval"
  }
  
  # Extract only the method value (no metric info)
  # based on a given method string
  method <- gsub("_score|_pval", "", method)
  
  # Create a list of known arguments (excluding FS_mat and input_score)
  # that can be passed to custom_function()
  known_parameters <- list(method = method, alternative = alternative, 
                           weight = weight, seed_names = seed_names, 
                           ...)
  
  # Select the appropriate method to compute row-wise directional scores
  rscore <- switch(
    method,
    ks = ks_rowscore(
      FS_mat = FS_mat,
      input_score = input_score,
      weight = weight,
      alternative = alternative,
      metric = metric
    ),
    wilcox = wilcox_rowscore(
      FS_mat = FS_mat,
      input_score = input_score,
      alternative = alternative,
      metric = metric
    ),
    revealer = revealer_rowscore(
      FS_mat = FS_mat,
      input_score = input_score,
      seed_names = seed_names,
      assoc_metric = "IC"
    ),
    custom = custom_rowscore(
      FS_mat = FS_mat,
      input_score = input_score,
      custom_function = custom_function,
      custom_parameters = custom_parameters,
      known_parameters = known_parameters
    )
  )

  return(rscore)

}



