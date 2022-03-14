
#' Customized Scoring Method
#' 
#' Compute row-wise scoring for each row of a given binary feature matrix
#' @param mat a matrix of binary features (required). The provided \code{mat} along with \code{custom_parameters} will be passed as arguments to custom_function() which is used to compute row-wise scoring for each row of a given binary feature.
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
#' # Provide a vector of continuous scores for a target profile
#' target = rnorm(n = ncol(sim.ES))
#'  
#' # Define additional parameters and start the function
#' custom_genescore_result <- custom_genescore_mat(
#'   mat = mat,
#'   custom_function = customized_genescore_mat, 
#'   custom_parameters = list(target = target, alternative = "less")
#' )
#'  
#' @export
#' @importFrom purrr map_dfr
custom_genescore_mat <- function
(
  mat,
  custom_function,
  custom_parameters,
  verbose = FALSE
)
{
  
  # Setup verbose option definition
  options(verbose=FALSE)
  
  # Check if the matrix has only binary values and no empty values
  if(length(mat) == 0 || !is.matrix(mat) || any(!mat %in% c(0,1)) || any(is.na(mat)))
    stop("mat variable must be a matrix of binary values (no empty values).\n")
  
  # Make sure the mat variable has rownames for features tracking
  if(is.null(rownames(mat)))
    stop("The mat object does not have rownames or featureData to track the features by. Please provide unique features or rownames for the expression matrix.\n")
  
  # Check if the dataset has any all 0 or 1 features (these are to be removed since they are not informative)
  if(any(rowSums(mat) == 0) || any(rowSums(mat) == ncol(mat))){
    warning("Provided dataset has features that are either all 0 or 1. These features will be removed from the computation.\n")
    mat <- mat[!(rowSums(mat) == 0 | rowSums(mat) == ncol(mat)),]
  }
  
  if(nrow(mat) == 0){
    stop("After removing features that are either all 0 or 1. There are no more features remained for downsteam computation.\n")
  }
  
  if(nrow(mat) < 2)
    warning("You are computing a row-wise statistic over a matrix with nrow < 2.\n")
  
  # check if the custom_function is indeed a function
  if(!is.function(custom_function)){
    stop("custom_function parameter must be a function.")
  }
  
  # check if the custom_parameters contains a list with arguments
  if(!is.list(custom_parameters)){
    stop("custom_parameters must be a list that contains arguments to pass to custom_function().")
  }else{
    if(is.null(names(custom_parameters))){
      stop("custom_parameters must be a list with labels or names to each value.")
    }else{
      if(length(names(custom_parameters)) != length(unique(names(custom_parameters)))){
        stop("custom_parameters must be a list with unique labels or names.")
      }
    }
  }
  
  if("mat" %in% names(custom_parameters)){
    warning("Removing 'mat' from the custom_parameter list. Using the required 'mat' as argument instead.")
    custom_parameters <- within(custom_parameters, rm(mat))
  }
  
  ## create a list for the data matrix 
  mat_list = list(mat = mat)
  
  ## append a list of provided parameters to matrix list
  custom_parameters <- append(mat_list, custom_parameters)

  ## check if the function runs with no errors
  custom <- tryCatch({
    base::do.call(custom_function, custom_parameters)
  }, error = function(err){
    stop("There is an error while trying to perform custom_function(). Check your function and the provided paramaters again!")
  })
  
  ## check if the custom is a data frame
  if(!is.data.frame(custom)){
    stop("The custom function must be a data frame with one or two columns: 'score' or 'p_value' or 'both'.")
  }else{
    if(all(!c("score", "p_value") %in% colnames(custom))){
      stop("The custom function must return a data frame with one or two columns: 'score' or 'p_value' or 'both'.")
    }
  }
  
  return(custom)
  
}

