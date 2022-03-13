
#' \code{REVEALER} Scoring Method
#' 
#' Compute conditional mutual information of \code{x} and \code{y} given \code{z} for each row of a given binary feature matrix
#' @param mat a matrix of binary features (required)
#' @param target a vector of continuous values of a target profile (required). target must include labels or names that associated with the colnames of the feature matrix. 
#' @param target_match a direction of target matching (\code{"negative"} or \code{"positive"}). Use \code{"positive"} to match the higher values of the target, \code{"negative"} to match the lower values. Default is \code{positive}. 
#' @param seed_names one or more features(s) that associated with the activation of a given target profile
#' @param seed_combination_op An operation to consolidate and summarize vectors of seed_names (\code{"max"} or \code{"min"} or \code{"mean"} or \code{"median"}). Default is \code{max}.
#' @param assoc_metric an association metric: \code{"IC"} information coefficient or \code{"COR"} correlation. Default is \code{IC}.
#' @param verbose a logical value indicates whether or not to print the diagnostic messages. Default is \code{FALSE}. 
#'
#' @return a data frame with two columns: \code{score} and \code{p_value}
#' @examples
#' # Load R library
#' library(Biobase)
#'
#' # Load pre-computed expression set
#' data(sim.ES)
#' 
#' # Provide a vector of continuous scores for a target profile
#' target = rnorm(n = ncol(sim.ES))
#' names(target) <- colnames(sim.ES)
#' 
#' # Define additional parameters and start the function
#' revealer_genescore_result <- revealer_genescore_mat(
#'   mat = exprs(sim.ES), target = target, target_match = "positive", assoc_metric = "IC"
#' )
#'  
#' @export
#' @importFrom purrr map_dfr
revealer_genescore_mat <- function
(
  mat,                                   
  target, 
  target_match = "positive",             
  seed_names = NULL,
  seed_combination_op = "max", 
  assoc_metric = "IC",
  verbose = FALSE
)
{
  
  # Setup verbose option definition
  options(verbose=FALSE)
  
  # Check if the matrix has only binary values and no empty values
  if(length(mat) == 0 || !is.matrix(mat) || any(!mat %in% c(0,1)) || any(is.na(mat)))
    stop("mat variable must be a matrix with binary values (no empty values).\n")
  
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
    if(any(names(target) != colnames(mat))){
      stop("The target object must have names or labels that matches the colnames of the expression matrix.\n")
    }
    mat <- mat[,names(target)]
  }
  
  # Check if seed_names is provided
  if(length(seed_names) > 0){
    if(any(!seed_names %in% rownames(mat))){
      stop("The provided seed_names must be part of the features or rownames of the expression matrix.\n")
    }   
  }else{
    seed_names <- "NULLSEED"   
  }
  
  # If target_match variable is not specified, use "positive" as default.
  if(length(target_match) == 0 || nchar(target_match) == 0){
    warning("The target_match variable is not specified. Using 'positive' by default ..\n")
    target_match <- "positive"
  }else if(length(target_match) == 1 && !target_match %in% c("positive", "negative")){
    stop(paste0(target_match, collapse=", "), " is not a valid target_match value. The target_match variable must be 'positive' or 'positive'.")
  }else if(length(target_match) > 1 && all(!target_match %in% c("positive", "negative"))){
    stop(paste0(target_match, collapse=", "), " is not a valid target_match value. The target_match variable must be 'positive' or 'positive'.")
  }else if(length(target_match) > 1 && any(target_match %in% c("positive", "negative"))){
    target_match <- target_match[which(target_match %in% c("positive", "negative"))][1]
    warning("More than one target_match values were specified. Only the first valid target_match value, '", target_match, "', is used.\n")
  }
  
  # If seed_combination_op variable is not specified, use "max" as default.
  if(length(seed_combination_op) == 0 || nchar(seed_combination_op) == 0){
    warning("The seed_combination_op variable is not specified. Using 'max' by default ..\n")
    seed_combination_op <- "max"
  }else if(length(seed_combination_op) == 1 && !seed_combination_op %in% c("max", "min", "mean", "median")){
    stop(paste0(seed_combination_op, collapse=", "), " is not a valid seed_combination_op value. The seed_combination_op variable must be 'max' or 'min' or 'mean' or 'median'.")
  }else if(length(seed_combination_op) > 1 && all(!seed_combination_op %in% c("max", "min", "mean", "median"))){
    stop(paste0(seed_combination_op, collapse=", "), " is not a valid seed_combination_op value. The seed_combination_op variable must be 'max' or 'min' or 'mean' or 'median'.")
  }else if(length(seed_combination_op) > 1 && any(seed_combination_op %in% c("max", "min", "mean", "median"))){
    seed_combination_op <- seed_combination_op[which(seed_combination_op %in% c("max", "min", "mean", "median"))][1]
    warning("More than one seed_combination_op values were specified. Only the first valid seed_combination_op value, '", seed_combination_op, "', is used.")
  }
  
  # Check if the dataset has any all 0 or 1 features (these are to be removed since they are not informative)
  if(any(rowSums(mat) == 0) || any(rowSums(mat) == ncol(mat))){
    warning("The provided matrix has some features that are either all 0 or 1. These features will be removed from the computation.\n")
    mat <- mat[!(rowSums(mat) == 0 | rowSums(mat) == ncol(mat)),]
  }
  
  if(nrow(mat) < 2)
    warning("You are computing a row-wise statistic over a matrix with nrow < 2.\n")
  
  # Define seed from given seed_names
  if (seed_names == "NULLSEED") {
    seed <- as.vector(rep(0, ncol(mat)))      
  } else {
    if (length(seed_names) > 1) {
      seed <- apply(mat[seed_names,], MARGIN=2, FUN=seed_combination_op)
    } else {
      seed <- mat[seed_names,]
    }
    locs <- match(seed_names, row.names(mat))
    mat <- mat[-locs,]
  }
  
  # Compute MI and % explained with original seed(s)
  cmi <- 1:nrow(mat) %>% 
    purrr::map_dfr(
      function(r){
        #r=1;
        revealer_genescore(x=target, y=mat[r,], z=seed, assoc_metric=assoc_metric, target_match=target_match) 
      }
    )
  
  # Only score value return from revealer
  colnames(cmi) <- c("score")
  rownames(cmi) <- rownames(mat)
  
  return(cmi)
  
}

