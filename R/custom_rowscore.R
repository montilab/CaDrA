
#' Customized Scoring Method
#'
#' Compute a row-wise score for each row of a given binary feature matrix
#' using a custom-defined function
#'
#' @param FS_mat a matrix of binary features where rows represent features of 
#' interest (e.g. genes, transcripts, exons, etc...) and columns represent 
#' the samples.
#' @param input_score a vector of continuous scores representing a phenotypic
#' readout of interest such as protein expression, pathway activity, etc.
#' The \code{input_score} object must have names or labels that match the column
#' names of \code{FS_mat} object.
#' @param custom_function a customized function which computes a row-wise 
#' score for each row of a given binary feature matrix (FS_mat).
#' 
#' NOTE: custom_function() must take FS_mat (or FS) and input_score as its
#' input arguments, and its final result must return a vector of row-wise scores 
#' ordered from most significant to least significant where its labels or names 
#' matched the row names of FS_mat object.
#' @param custom_parameters a list of additional arguments to be passed to  
#' custom_function() (excluding \code{FS_mat} (or FS) and \code{input_score}).
#' Default is NULL.
#' @param known_parameters a list of known parameters that existed in
#' the previous function that can be passed to custom_function() if and only 
#' if they were not provided by custom_parameters. Default is NULL.
#' 
#' @noRd
#' 
#' @examples 
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
#' # Search for best features using a custom-defined function
#' custom_rs <- custom_rowscore(
#'   FS_mat = mat,
#'   input_score = input_score,
#'   custom_function = customized_rowscore,            
#'   custom_parameters = NULL  
#' )
#' 
#' @return return a vector of row-wise scores where its labels or names 
#' must match the row names of \code{FS_mat} object
#' 
custom_rowscore <- function
(
  FS_mat,
  input_score,
  custom_function,
  custom_parameters = NULL,
  known_parameters = NULL
)
{
  
  # Check if the custom_function is indeed a function
  if(!is.function(custom_function))
    stop("custom_function must be a function.")
  
  # If custom_parameters is provided, check if it is a list and 
  # has labels or names that associated with each of its value
  # e.g. custom_parameters = list(alternative = 'less')
  if((!is.null(custom_parameters) && !is.list(custom_parameters)) || 
     (!is.null(custom_parameters) && is.list(custom_parameters) && 
      is.null(names(custom_parameters))))    
    stop("custom_parameters must be a list with labels or names ",
         "that attach to each of its values. \nFor example: ",
         "custom_parameters = list(alternative = 'less')")
  
  # Extract all formal arguments required by custom_function()
  custom_args <- formals(custom_function)
  
  # Check if custom_function() requires 'FS_mat' or "FS" as its argument
  if(all(!c("FS_mat", "FS") %in% names(custom_args)))
    stop("custom_function() must take 'FS_mat' (or 'FS') as ",
         "one of its arguments (required).")
  
  # Check if custom_function() requires 'input_score' as its argument
  if(!"input_score" %in% names(custom_args))
    stop("custom_function() must take 'input_score' ",
         "as one of its arguments (required).")
  
  ## Create a list with only the required variables 
  req_args <- list(FS_mat=FS_mat, input_score=input_score)
  
  # Combine custom_parameters, required variables, and a list of 
  # known parameters together 
  # However, we must exclude FS_mat, input_score from custom_parameters and 
  # excluding FS_mat, input_score, and custom_parameters from known parameters 
  # as they would be redundant
  combined_parameters <- c(
    req_args, 
    custom_parameters[
      which(!names(custom_parameters) %in% names(req_args))],
    known_parameters[
      which(!names(known_parameters) %in% c(names(req_args), names(custom_parameters)))]
  )
  
  # Extract a list of custom_function() parameters that existed in combined variables
  included_parameters <- combined_parameters[
    which(names(combined_parameters) %in% names(custom_args))] 
  
  # Check if some parameters not existed in a list of combined variables
  excluded_parameters <- custom_args[
    which(!names(custom_args) %in% names(combined_parameters))]
  
  # If some parameters are excluded, check to see if that argument has a default value
  # Finally, return all necessarily arguments to be passed to custom_function()
  all_parameters <- c(
    included_parameters, 
    lapply(seq_along(excluded_parameters), function(s){
      #s=2;
      if(!is.null(excluded_parameters[[s]]) && excluded_parameters[[s]] == "")
        stop("argument '", names(excluded_parameters)[s], 
             "' is missing from custom_function() with no default value") 
      else
        return(excluded_parameters[s])
    }) |> unlist())
  
  ## Check if the function runs with no errors
  custom <- tryCatch({
    base::do.call(custom_function, all_parameters)
  }, error = function(err){
    stop(err)
  })

  # Make sure custom_function() returns a vector of scores with no NAs
  # where it has labels or feature names that match the row names of FS_mat 
  # (or FS) object
  if(length(custom) == 0 || any(!is.numeric(custom)) || any(is.na(custom)) || 
     is.null(names(custom)) || any(!names(custom) %in% rownames(FS_mat)))
    stop("The custom_function must return a vector of continuous scores ",
         "(with no NAs) where its feature names or labels must match ",
         "the row names of FS_mat object.\n")
  
  return(custom)
  
}

