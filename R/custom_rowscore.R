
#' Customized Scoring Method
#'
#' Compute row-wise scoring for each row of a given binary feature matrix
#' using a custom-defined function
#'
#' @param FS a feature set of binary features. It can be a matrix or
#' a \code{SummarizedExperiment} class object from SummarizedExperiment package.
#' @param input_score a vector of continuous scores representing a phenotypic
#' readout of interest such as protein expression, pathway activity, etc.
#' The \code{input_score} object must have names or labels that match the column
#' names of FS object.
#' @param custom_function a customized function to perform row-wise scoring.
#' NOTE: this function must take FS and input_score as input arguments, and its final
#' result must return a matrix with one or two columns:
#' \code{score} or \code{p_value} or \code{both (score and p_value)}.
#' @param custom_parameters a list of additional arguments to be passed to
#' the custom_function (excluding \code{FS} and \code{input_score}).
#' @param warning a logical value indicates whether or not to print the 
#' diagnostic messages. Default is \code{TRUE}
#' 
#' @noRd
#'
#' @return a matrix with one or two columns: \code{score} or
#' \code{p_value} or \code{both} (score and p_value)
#' @import SummarizedExperiment
custom_rowscore <- function
(
  FS,
  input_score,
  custom_function,
  custom_parameters = NULL,
  warning = TRUE
)
{
  
  # Check of FS and input_score are valid inputs
  if(warning == TRUE) check_data_input(FS = FS, input_score = input_score, max_size = 1)
  
  # Get the feature names
  feature_names <- rownames(FS)
  
  # Extract the feature binary matrix
  if(is(FS, "SummarizedExperiment")){
    mat <- as.matrix(SummarizedExperiment::assay(FS))
  }else if(is(FS, "matrix")){
    mat <- as.matrix(FS)
  }else{
    mat <- matrix(t(FS), nrow=1, byrow=TRUE,
                  dimnames=list(feature_names, names(FS)))
  }

  # check if the custom_function is indeed a function
  if(!is.function(custom_function)){
    stop("custom_function must be a function.")
  }

  # if custom_parameters is not NULL, check if the custom_parameters
  # contains a list with arguments
  if(!is.null(custom_parameters)){
    if(!is.list(custom_parameters)){
      stop("custom_parameters must be a list that contains arguments to pass ",
           "to custom_function().")
    }else{
      if(is.null(names(custom_parameters))){
        stop("custom_parameters must be a list with labels or names ",
             "that attach to each of its values.")
      }else{
        if(length(names(custom_parameters)) !=
           length(unique(names(custom_parameters)))){
          stop("custom_parameters must be a list with unique labels or ",
               "names that attach to each of its values.")
        }
      }
    }
  }

  # check if custom_function() required FS and input_score as arguments
  custom_args <- names(formals(custom_function))

  # if custom_args does not contain 'FS' as parameter, gives a warning
  if(!"FS" %in% custom_args){
    stop("custom_function() must contain 'FS' as ",
         "one of its arguments (required).")
  }

  # if custom_args does not contain 'input_score' as parameter, gives a warning
  if(!"input_score" %in% custom_args){
    stop("custom_function() must contain 'input_score' ",
         "as one of its arguments (required).")
  }

  # using the required FS argument as parameter. removing any FS variables
  # from the custom_parameter list since it is redundant.
  if("FS" %in% names(custom_parameters)){
    warning("Removing 'FS' from the custom_parameter list. ",
            "Using the provided 'FS' argument instead.")
    custom_parameters <- within(custom_parameters, rm(FS))
  }

  # using the required input_score argument as parameter.
  # removing any input_score variables from the custom_parameter
  # list since it is redundant.
  if("input_score" %in% names(custom_parameters)){
    warning("Removing 'input_score' from the custom_parameter list. ",
            "Using the provided 'input_score' argument instead.")
    custom_parameters <- within(custom_parameters, rm(input_score))
  }

  ## create a list for the required binary data matrix and input_score variables
  mat_list <- list(FS = FS, input_score = input_score)

  ## combine a FS and input_score variables with additional parameters
  # provided in the custom_parameters list
  custom_parameters <- base::append(mat_list, custom_parameters)

  # check if there are any missing arguments needed to run custom_function()
  if(!all(custom_args %in% names(custom_parameters))){
    stop(paste0("Missing parameters: ",
                paste0(custom_args[which(!custom_args %in%
                                           names(custom_parameters))],
                       collapse = ", "), ", needed in custom_function()."))
  }

  # only extract the required arguments to passed to the custom_function()
  custom_parameters <-
    custom_parameters[which(names(custom_parameters) %in% custom_args)]

  ## check if the function runs with no errors
  custom <- tryCatch({
    base::do.call(custom_function, custom_parameters)
  }, error = function(err){
    stop(err)
  })

  ## check if the custom is a matrix
  if(!is.matrix(custom)){
    stop("The custom function must return a matrix with one or ",
         "two columns: 'score' or 'p_value' or 'both'.")
  }else if(all(!c("score", "p_value") %in% colnames(custom))){
    stop("The custom function must return a matrix with one or ",
         "two columns: 'score' or 'p_value' or 'both'.")
  }

  return(custom)

}

