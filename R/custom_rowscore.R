
#' Customized Scoring Method
#'
#' Compute row-wise scores for each row of a given binary feature matrix
#' using a custom-defined function
#'
#' @param FS_mat a matrix of binary features where 
#' rows represent features of interest (e.g. genes, transcripts, exons, etc...)
#' and columns represent the samples.
#' @param input_score a vector of continuous scores representing a phenotypic
#' readout of interest such as protein expression, pathway activity, etc.
#' The \code{input_score} object must have names or labels that match the column
#' names of FS_mat object.
#' @param custom_function a customized function to perform row-wise scoring.
#' NOTE: this function must take FS_mat and input_score as input arguments, 
#' and its final result must return a matrix with one or two columns:
#' \code{score} or \code{score} and \code{p_value}.
#' @param custom_parameters a list of additional arguments to be passed to
#' the custom_function (excluding \code{FS_mat} and \code{input_score}).
#' 
#' @noRd
#'
#' @return return a vector of scores ordered from most significant to least
#' significant where its labels or names match the row names of FS_mat object
#' 
custom_rowscore <- function
(
  FS_mat,
  input_score,
  custom_function,
  custom_parameters = NULL
)
{
  
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

  # check if custom_function() required FS_mat and input_score as arguments
  custom_args <- names(formals(custom_function))

  # if custom_args does not contain 'FS_mat' as parameter, gives a warning
  if(!"FS_mat" %in% custom_args){
    stop("custom_function() must contain 'FS_mat' as ",
         "one of its arguments (required).")
  }

  # if custom_args does not contain 'input_score' as parameter, gives a warning
  if(!"input_score" %in% custom_args){
    stop("custom_function() must contain 'input_score' ",
         "as one of its arguments (required).")
  }

  # using the required FS_mat argument as parameter. removing any FS_mat variables
  # from the custom_parameter list since it is redundant.
  if("FS_mat" %in% names(custom_parameters)){
    warning("Removing 'FS_mat' from the custom_parameter list. ",
            "Using the provided 'FS_mat' argument instead.")
    custom_parameters <- within(custom_parameters, rm(FS_mat))
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
  mat_list <- list(FS_mat = FS_mat, input_score = input_score)

  ## combine a FS_mat and input_score variables with additional parameters
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

  ## Check if the function runs with no errors
  custom <- tryCatch({
    base::do.call(custom_function, custom_parameters)
  }, error = function(err){
    stop(err)
  })

  ## Make sure custom function returns a vector of scores with no NAs
  if(length(custom) == 0 || any(!is.numeric(custom)) || any(is.na(custom)))
    stop("The custom function must return a vector of continuous scores ",
         "(with no NAs) where it has names or labels that match the row names ",
         "(or feature names) of the FS_mat object.\n")
  
  # Make sure the custom has names or labels that are the
  # same as the feature names as the FS_mat object
  if(is.null(names(custom)))
    stop("The custom function must return a vector of continuous scores ",
         "(with no NAs) where it has names or labels that match the row names ",
         "(or feature names) of the FS_mat object.\n")
  
  # Make sure the custom has the same length as number of features in FS_mat
  if(length(custom) != nrow(FS_mat))
    stop("The custom function must return a vector of continuous scores that has ",
         "the same length as the number of rows in FS_mat object.\n")
  
  # Make sure the custom has names or labels that are the
  # same as the feature names as the FS_mat object  
  if(any(!names(custom) %in% rownames(FS_mat)))
    stop("The custom object must have names or labels that match the feature names ",
         "(or row names) of FS_mat object.")
  
  return(custom)
  
}

