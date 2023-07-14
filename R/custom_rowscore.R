
#' Customized Scoring Method
#'
#' Compute a row-wise score for each row of a given binary feature matrix
#' using a custom-defined function
#'
#' @param FS a matrix of binary features or a SummarizedExperiment class object 
#' from SummarizedExperiment package where rows represent features of interest 
#' (e.g. genes, transcripts, exons, etc...) and columns represent the samples. 
#' The assay of FS contains binary (1/0) values indicating the presence/absence 
#' of omics features.
#' @param input_score a vector of continuous scores representing a phenotypic
#' readout of interest such as protein expression, pathway activity, etc.
#' The \code{input_score} object must have names or labels that match the column
#' names of \code{FS} object.
#' @param seed_names a vector of one or more features representing known 
#' causes of activation or features associated with a response of interest, 
#' \code{e.g. input_score}. Default is NULL.
#' @param custom_function a customized function which computes a row-wise 
#' score for each row of a given binary feature set (FS).
#' 
#' NOTE: \code{custom_function()} must take \code{FS} and \code{input_score} as 
#' its input arguments, and its final result must return a vector of row-wise 
#' scores ordered from most significant to least significant where its labels 
#' or names matched the row names of FS object.
#' @param custom_parameters a list of additional arguments to be passed to  
#' \code{custom_function()} (excluding \code{FS} and \code{input_score}).
#' Default is NULL.
#' @param ... additional parameters to be passed to \code{custom_function}
#' 
#' @noRd
#' 
#' @examples 
#' 
#' # A customized function using ks-test function
#' customized_rowscore <- function(FS, input_score, seed_names=NULL, alternative="less"){
#'   
#'   # Check if seed_names is provided
#'   if(!is.null(seed_names)){
#'     # Taking the union across the known seed features
#'     if(length(seed_names) > 1) {
#'       seed_vector <- as.numeric(ifelse(colSums(FS[seed_names,]) == 0, 0, 1))
#'     }else{
#'       seed_vector <- as.numeric(FS[seed_names,])
#'     }
#'      
#'     # Remove the seeds from the binary feature matrix
#'     # and taking logical OR btw the remaining features with the seed vector
#'     locs <- match(seed_names, row.names(FS))
#'     FS <- base::sweep(FS[-locs,], 2, seed_vector, `|`)*1
#'      
#'     # Check if there are any features that are all 1s generated from
#'     # taking the union between the matrix
#'     # We cannot compute statistics for such features and thus they need
#'     # to be filtered out
#'     if(any(rowSums(FS) == ncol(FS))){
#'       warning("Features with all 1s generated from taking the matrix union ",
#'               "will be removed before progressing...\n")
#'       FS <- FS[rowSums(FS) != ncol(FS),]
#'     }
#'   }
#'    
#'   # KS is a ranked-based method
#'   # So we need to sort input_score from highest to lowest values
#'   input_score <- sort(input_score, decreasing=TRUE)
#'    
#'   # Re-order the matrix based on the order of input_score
#'   FS <- FS[, names(input_score), drop=FALSE]  
#'   
#'   # Compute the scores using the KS method
#'   ks <- apply(FS, 1, function(r){ 
#'     x = input_score[which(r==1)]; 
#'     y = input_score[which(r==0)];
#'     res <- ks.test(x, y, alternative=alternative)
#'     return(c(res$statistic, res$p.value))
#'   })
#'   
#'   # Obtain score statistics
#'   stat <- ks[1,]
#'   
#'   # Obtain p-values and change values of 0 to the machine lowest value 
#'   # to avoid taking -log(0)
#'   pval <- ks[2,]
#'   pval[which(pval == 0)] <- .Machine$double.xmin
#'   
#'   # Compute the -log(pval)
#'   # Make sure scores has names that match the row names of FS object
#'   scores <- -log(pval)
#'   names(scores) <- rownames(FS)
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
#'   FS = mat,
#'   input_score = input_score,
#'   seed_names = NULL,
#'   custom_function = customized_rowscore,            
#'   custom_parameters = NULL  
#' )
#' 
#' @return return a vector of row-wise positive scores where its labels or names 
#' must match the row names of \code{FS} object
#' 
custom_rowscore <- function
(
  FS,
  input_score,
  seed_names = NULL,
  custom_function,
  custom_parameters = NULL,
  ...
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
  
  # Check if custom_function() requires 'FS' as its argument
  if(!"FS" %in% names(custom_args))
    stop("custom_function() must take 'FS'as ",
         "one of its arguments (required).")
  
  # Check if custom_function() requires 'input_score' as its argument
  if(!"input_score" %in% names(custom_args))
    stop("custom_function() must take 'input_score' ",
         "as one of its arguments (required).")
  
  # Check if custom_function() requires 'seed_names' as its argument
  if(!"seed_names" %in% names(custom_args))
    stop("custom_function() must take 'seed_names' ",
         "as one of its arguments (required).")
  
  ## Create a list with only the required variables 
  req_parameters <- list(FS=FS, input_score=input_score, seed_names=seed_names)
  
  # Obtain additional parameters
  additional_parameters <- list(...)
  
  # Combine custom_parameters, required variables, and a list of 
  # known parameters together 
  # However, we must exclude FS, input_score from custom_parameters 
  # excluding FS, input_score, and custom_parameters from known parameters 
  # as they would be redundant
  combined_parameters <- c(
    req_parameters, 
    custom_parameters[
      which(!names(custom_parameters) %in% names(req_parameters))],
    additional_parameters[
      which(!names(additional_parameters) %in% c(names(req_parameters), 
                                            names(custom_parameters)))]
  )
  
  # Extract a list of custom_function() parameters that existed in 
  # combined variables
  included_parameters <- combined_parameters[
    which(names(combined_parameters) %in% names(custom_args))] 
  
  # Check if some parameters not existed in a list of combined variables
  excluded_parameters <- custom_args[
    which(!names(custom_args) %in% names(combined_parameters))]
  
  # If some parameters are excluded, 
  # check to see if that argument has a default value
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
  # where it has labels or feature names that match the row names of FS 
  # (or FS) object
  if(length(custom) == 0 || any(!is.numeric(custom)) || any(is.na(custom)) || 
     is.null(names(custom)) || any(!names(custom) %in% rownames(FS)))
    stop("The custom_function must return a vector of continuous scores ",
         "(with no NAs) where its feature names or labels must match ",
         "the row names of FS object.\n")
  
  return(custom)
  
}

