
#' Customized Scoring Method
#' 
#' Compute row-wise scoring for each row of a given binary feature matrix
#' @param mat a matrix of binary features (required). \code{NOTE:} The provided \code{mat} along with \code{target} and \code{custom_parameters} will be passed as arguments to custom_function() which is later used to compute row-wise scoring for each row of a given binary feature.
#' @param target a vector of continuous values for a target profile (required). \code{target} must include labels or names that associated with the colnames of the binary feature matrix. \code{NOTE:} \code{target} will be passed as one of the arguments to custom_function().
#' @param custom_function a customized function to perform row-wise scoring for each row of a given binary feature matrix (required). \code{NOTE:} This function must return a data frame with one or two columns: \code{score} or \code{p_value} or \code{both}.
#' @param custom_parameters a list of additional arguments to be passed to the custom_function() (exluding \code{mat} and \code{target} parameters).
#' @param verbose a logical value indicates whether or not to print the diagnostic messages. Default is \code{FALSE}. 
#'
#' @return a data frame with one or two columns: \code{score} or \code{p_value} or \code{both}
#' @examples
#' 
#' # Examples of a customized function using ks-test function
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
#' # Load R library
#' library(Biobase)
#' 
#' # Load pre-computed expression set
#' data(sim.ES)
#' 
#' # Extract the binary feature matrix
#' mat = exprs(sim.ES)
#' 
#' # set seed
#' set.seed(123)
#' 
#' # Provide a vector of continuous scores for a target profile
#' target = rnorm(n = ncol(sim.ES))
#' names(target) <- colnames(sim.ES)
#' 
#' # Define additional parameters and start the function
#' custom_genescore_result <- custom_genescore_mat(
#'   mat = mat,
#'   target = target,
#'   custom_function = customized_genescore_mat, 
#'   custom_parameters = list(alternative = "less")
#' )
#'  
#' @export
#' @importFrom purrr map_dfr
custom_genescore_mat <- function
(
  mat,
  target,
  custom_function,
  custom_parameters = NULL,
  verbose = FALSE
)
{
  
  # Setup verbose option definition
  options(verbose=FALSE)
  
  # Check if the matrix has only binary values and no empty values
  if(length(mat) == 0 || !is.matrix(mat) || any(!mat %in% c(0,1)) || any(is.na(mat)))
    stop("mat variable must be a matrix of binary values (no empty values).\n")
  
  # Check if target is provided and no empty values
  if(length(target) == 0 || any(!is.numeric(target)) || any(is.na(target)))
    stop("target variable must be provided and are numeric with no empty values.\n")
  
  # Make sure the mat variable has rownames for features tracking
  if(is.null(rownames(mat)))
    stop("The mat object does not have rownames or featureData to track the features by. Please provide unique features or rownames for the expression matrix.\n")
  
  # Make sure the target variable has names as the colnames in mat
  if(is.null(names(target)))
    stop("The target object must have names or labels to track the samples by. Please provide the sample names or labels that matches the colnames of the expression matrix.\n")
  
  # Make sure the target has the same length as number of samples in mat
  if(length(target) != ncol(mat)){
    stop("The target variable must have the same length as the number of columns in mat.\n")
  }else{
    # check if target has any labels or names
    if(length(names(target)) == 0){
      stop("The target object must have names or labels that match the colnames of the expression matrix.\n")
    }
    
    # check if target has labels or names that matches the colnames of the expression matrix
    if(any(!names(target) %in% colnames(mat))){
      stop("The target object have names or labels that do not match the colnames of the expression matrix.\n")
    }
    # match colnames of expression matrix with names of provided target values
    # if nrow = 1, if it is, convert to matrix form as it is needed for backward_forward_search with one dimension matrix computation
    if(nrow(mat) == 1){
      mat <- matrix(t(mat[,names(target)]), nrow=1, byrow=T, dimnames = list(rownames(mat), colnames(mat)))
    }else{
      mat <- mat[,names(target)]
    }
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
  
  # check if the custom_function is indeed a function
  if(!is.function(custom_function)){
    stop("custom_function parameter must be a function.")
  }
  
  # if custom_parameters is not NULL, check if the custom_parameters contains a list with arguments
  if(!is.null(custom_parameters)){
    if(!is.list(custom_parameters)){
      stop("custom_parameters must be a list that contains arguments to pass to custom_function().")
    }else{
      if(is.null(names(custom_parameters))){
        stop("custom_parameters must be a list with labels or names that attach to each of its values.")
      }else{
        if(length(names(custom_parameters)) != length(unique(names(custom_parameters)))){
          stop("custom_parameters must be a list with unique labels or names that attach to each of its values.")
        }
      }
    }
  }
  
  # check if custom_function() required mat and target as arguments
  custom_args <- names(formals(custom_function))
  
  # if custom_args does not contain 'mat' as parameter, gives a warning
  if(!"mat" %in% custom_args){
    stop("custom_function() must contain 'mat' as one of its arguments (required).")
  }
  
  # if custom_args does not contain 'target' as parameter, gives a warning
  if(!"target" %in% custom_args){
    stop("custom_function() must contain 'target' as one of its arguments (required).")
  }
  
  # using the required mat argument as parameter. removing any mat variables from the custom_parameter list since it is redundant. 
  if("mat" %in% names(custom_parameters)){
    warning("Removing 'mat' from the custom_parameter list. Using the required 'mat' as argument instead.")
    custom_parameters <- within(custom_parameters, rm(mat))
  }
  
  # using the required target argument as parameter. removing any target variables from the custom_parameter list since it is redundant.
  if("target" %in% names(custom_parameters)){
    warning("Removing 'target' from the custom_parameter list. Using the required 'target' as argument instead.")
    custom_parameters <- within(custom_parameters, rm(target))
  }
  
  ## create a list for the required binary data matrix and target variables
  mat_list = list(mat = mat, target = target)
  
  ## combine a mat and target variables with additional parameters provided in the custom_parameters list
  custom_parameters <- base::append(mat_list, custom_parameters)
  
  # check if there are any missing arguments needed to run custom_function()
  if(!all(custom_args %in% names(custom_parameters))){
    stop(paste0("Missing parameters: ", paste0(custom_args[which(!custom_args %in% names(custom_parameters))], collapse = ", "), ", needed in custom_function()."))
  }
  
  # only extract the required arguments to passed to the custom_function()
  custom_parameters <- custom_parameters[which(names(custom_parameters) %in% custom_args)]
  
  ## check if the function runs with no errors
  custom <- tryCatch({
    base::do.call(custom_function, custom_parameters)
  }, error = function(err){
    stop(err)
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

