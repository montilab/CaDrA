
#' Customized Scoring Method
#' 
#' Compute row-wise scoring for each row of a given binary feature matrix
#' @param custom_function a customized function to perform row-wise scoring for each row of a given binary feature matrix (required). \code{NOTE:} This function must return a data frame with one or two columns: \code{score} or \code{p_value} or \code{both}.
#' @param custom_parameters a list of arguments to be passed to the custom_function() (required).
#' @param verbose a logical value indicates whether or not to print the diagnostic messages. Default is \code{FALSE}. 
#'
#' @return a data frame with one or two columns: \code{score} and \code{p_value}
#' @examples
#' # Load R library
#' library(Biobase)
#' 
#' # Examples of customized function using ks-test function
#' customized_genescore_mat <- function(mat, target, alternative){
#'   result <- 1:nrow(mat) %>% 
#'     purrr::map_dfr(
#'       function(r){ 
#'         feature = mat[r,];
#'         x = target[which(feature==1)]; y = target[which(feature==0)];
#'         res <- ks.test(x, y, alternative=alternative)
#'         return(data.frame(score=res$statistic, p_value=res$p.value))
#'    })
#' }
#' 
#' # Load pre-computed expression set
#' data(sim.ES)
#' 
#' # Extract the binary feature matrix
#' mat = exprs(sim.ES)
#' 
#' # Check if the dataset has any all 0 or 1 features 
#' # These are to be removed since they are not informative
#' if(any(rowSums(mat) == 0) || any(rowSums(mat) == ncol(mat))){
#'   mat <- mat[!(rowSums(mat) == 0 | rowSums(mat) == ncol(mat)),]
#' }
#' 
#' # Provide a vector of continuous scores for a target profile
#' target = rnorm(n = ncol(sim.ES))
#'  
#' # Define additional parameters and start the function
#' custom_genescore_result <- custom_genescore_mat(
#'   custom_function = customized_genescore_mat, 
#'   custom_parameters = list(mat = mat, target = target, alternative = "less")
#' )
#'  
#' @export
#' @importFrom purrr map_dfr
custom_genescore_mat <- function
(
  custom_function,
  custom_parameters,
  verbose = FALSE
)
{
  
  # Setup verbose option definition
  options(verbose=FALSE)
  
  # check if the custom_function is indeed a function
  if(!is.function(custom_function)){
    stop("custom_function parameter must be a function.")
  }
  
  # check if the custom_parameters contains a list with arguments
  if(!is.list(custom_parameters)){
    stop("custom_parameters must be a list that contains arguments to pass to custom_function().")
  }
  
  custom <- tryCatch({
    base::do.call(custom_function, custom_parameters)
  }, error = function(err){
    stop("There is an error while trying to perform custom_function(). Check your function and the provided paramaters again!")
  })
  
  if(!is.data.frame(custom)){
    stop("The custom function must be a data frame with one or two columns: 'score' or 'p_value' or 'both'.")
  }else{
    if(any(!c("score", "p_value") %in% colnames(custom))){
      stop("The custom function must return a data frame with one or two columns: 'score' or 'p_value' or 'both'.")
    }
    custom <- custom[,c("score", "p_value") ]
  }
  
  return(custom)
  
}

