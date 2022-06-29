
#' Customized Scoring Method
#' 
#' Compute row-wise scoring for each row of a given binary feature matrix using a custom-defined function
#' 
#' @param mat a matrix of binary features (required). \code{NOTE:} The provided \code{mat} along with \code{input_score} and \code{custom_parameters} will be passed as arguments to custom_function() which is later used to compute row-wise scoring for each row of a given binary feature.
#' @param input_score a vector of continuous values for a response of interest (required). \code{input_score} must include labels or names that associated with the colnames of the binary feature matrix. \code{NOTE:} \code{input_score} will be passed as one of the arguments to custom_function().
#' @param custom_function a customized function to perform row-wise scoring for each row of a given binary feature matrix (required). \code{NOTE:} This function must return a data frame with one or two columns: \code{score} or \code{p_value} or \code{both}.
#' @param custom_parameters a list of additional arguments to be passed to the custom_function() (exluding \code{mat} and \code{input_score} parameters).
#' @param verbose a logical value indicates whether or not to print the diagnostic messages. Default is \code{FALSE}. 
#'
#' @return a data frame with one or two columns: \code{score} or \code{p_value} or \code{both}
#' @examples
#'
#' # Load R library
#' library(Biobase)
#'  
#' # Examples of a customized function using ks-test function
#' customized_genescore_mat <- function(mat, input_score, alternative){
#'   result <- 1:nrow(mat) %>% 
#'     purrr::map_dfr(
#'       function(r){ 
#'         feature = mat[r,];
#'         x = input_score[which(feature==1)]; y = input_score[which(feature==0)];
#'         res <- ks.test(x, y, alternative=alternative)
#'         return(data.frame(score=res$statistic, p_value=res$p.value))
#'    })
#' }
#' 
#' # Load pre-computed expression set
#' data(sim.ES)
#' 
#' mat = exprs(sim.ES)
#' 
#' # set seed
#' set.seed(123)
#' 
#' # Provide a vector of continuous scores for a target profile
#' input_score = rnorm(n = ncol(sim.ES))
#' names(input_score) <- colnames(sim.ES)
#' 
#' # Define additional parameters and start the function
#' custom_genescore_result <- custom_genescore_mat(
#'   mat = mat,
#'   input_score = input_score,
#'   custom_function = customized_genescore_mat, 
#'   custom_parameters = list(alternative = "less")
#' )
#'  
#' @export
#' @importFrom purrr map_dfr
custom_genescore_mat <- function
(
  mat,
  input_score,
  custom_function,
  custom_parameters = NULL,
  verbose = FALSE
)
{
  
  # Setup verbose option definition
  options(verbose=verbose)
  
  ## Make sure mat variable is a matrix
  mat <- as.matrix(mat)
  
  # If mat has only one column, it must be converted to a row-wise matrix form as it is needed for backward_forward_search() computation
  # mat must have rownames to track features and columns to track samples
  # for n = 1 case, it is only in backward_forward_search(), thus we can assign a random labels to it
  if(ncol(mat) == 1){
    mat <- matrix(t(mat), nrow=1, byrow=TRUE, dimnames = list("my_label", rownames(mat))) 
  }
  
  # Check if the matrix has only binary 0 or 1 values 
  if(length(mat) == 0 || !is.matrix(mat) || any(!mat %in% c(0,1)) || any(is.na(mat)))
    stop("mat variable must be a matrix of binary values (no empty values).\n")
  
  # Make sure the input mat has rownames for features tracking
  if(is.null(rownames(mat)))
    stop("The mat object does not have rownames or featureData to track the features by. Please provide unique features or rownames for the expression matrix.\n")
  
  # Check input_score is provided and is a continuous values with no NAs
  if(length(input_score) == 0 || any(!is.numeric(input_score)) || any(is.na(input_score)))
    stop("input_score must be a vector of continous values (with no NAs) where the vector names match the colnames of the expression matrix (required).\n")
  
  # Make sure the input_score has names or labels that are the same as the colnames of mat
  if(is.null(names(input_score)))
    stop("The input_score object must have names or labels to track the samples by. Please provide unique sample names or labels that matches the colnames of the expression matrix.\n")
  
  # Make sure the input_score has the same length as number of samples in mat
  if(length(input_score) != ncol(mat)){
    stop("The input_score variable must have the same length as the number of columns in mat.\n")
  }else{
    # check if input_score has any labels or names
    if(length(names(input_score)) == 0){
      stop("The input_score object must have names or labels that match the colnames of the expression matrix.\n")
    }
    
    # check if input_score has labels or names that matches the colnames of the expression matrix
    if(any(!names(input_score) %in% colnames(mat))){
      stop("The input_score object have names or labels that do not match the colnames of the expression matrix.\n")
    }
    
    # match colnames of expression matrix with names of provided input_score values
    # iif nrow = 1, if it is, convert to matrix form as it is needed for backward_forward_search with one dimension matrix computation
    if(nrow(mat) == 1){
      mat <- matrix(t(mat[,names(input_score)]), nrow=1, byrow=TRUE, dimnames = list(rownames(mat), colnames(mat))) 
    }else{
      mat <- mat[,names(input_score)]
    }
  }
  
  # Check if the dataset has any all 0 or 1 features (these are to be removed since they are not informative)
  if(any(rowSums(mat) == 0) || any(rowSums(mat) == ncol(mat))){
    warning("The provided matrix has some features that are either all 0 or 1. These features will be removed from downsteam computation.\n")
    mat <- mat[!(rowSums(mat) == 0 | rowSums(mat) == ncol(mat)),]
  }
  
  # Make sure matrix is not empty after removing uninformative features
  if(nrow(mat) == 0){
    stop("After removing features that are either all 0 or 1. There are no more features remained for downsteam computation.\n")
  }
  
  # Give a warning if matrix has nrow < 2
  if(nrow(mat) < 2)
    verbose("You are computing a row-wise statistic over a matrix with nrow < 2.\n")
  
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
  
  # check if custom_function() required mat and input_score as arguments
  custom_args <- names(formals(custom_function))
  
  # if custom_args does not contain 'mat' as parameter, gives a warning
  if(!"mat" %in% custom_args){
    stop("custom_function() must contain 'mat' as one of its arguments (required).")
  }
  
  # if custom_args does not contain 'input_score' as parameter, gives a warning
  if(!"input_score" %in% custom_args){
    stop("custom_function() must contain 'input_score' as one of its arguments (required).")
  }
  
  # using the required mat argument as parameter. removing any mat variables from the custom_parameter list since it is redundant. 
  if("mat" %in% names(custom_parameters)){
    warning("Removing 'mat' from the custom_parameter list. Using the required 'mat' as argument instead.")
    custom_parameters <- within(custom_parameters, rm(mat))
  }
  
  # using the required input_score argument as parameter. removing any input_score variables from the custom_parameter list since it is redundant.
  if("input_score" %in% names(custom_parameters)){
    warning("Removing 'input_score' from the custom_parameter list. Using the required 'input_score' as argument instead.")
    custom_parameters <- within(custom_parameters, rm(input_score))
  }
  
  ## create a list for the required binary data matrix and input_score variables
  mat_list = list(mat = mat, input_score = input_score)
  
  ## combine a mat and input_score variables with additional parameters provided in the custom_parameters list
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

