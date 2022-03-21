
#' \code{REVEALER} Scoring Method
#' 
#' Compute conditional mutual information of \code{x} and \code{y} given \code{z} for each row of a given binary feature matrix
#' @param mat a matrix of binary features (required)
#' @param target a vector of continuous values of a target profile (required). target must include labels or names that associated with the colnames of the feature matrix. 
#' @param target_match a direction of target matching (\code{"negative"} or \code{"positive"}). Use \code{"positive"} to match the higher values of the target, \code{"negative"} to match the lower values. Default is \code{positive}. 
#' @param seed_names one or more features(s) that associated with the activation of a given target profile
#' @param assoc_metric an association metric: \code{"IC"} information coefficient or \code{"COR"} correlation. Default is \code{IC}.
#' @param verbose a logical value indicates whether or not to print the diagnostic messages. Default is \code{FALSE}. 
#'
#' @return a data frame with one column: \code{score}
#' @examples
#' 
#' # Load R library
#' library(Biobase)
#'
#' # Load pre-computed expression set
#' data(sim.ES)
#' 
#' # set seed
#' set.seed(123)
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
  assoc_metric = c("IC", "COR"),
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
  
  # Check if seed_names is provided
  if(length(seed_names) > 0){
    if(any(!seed_names %in% rownames(mat))){
      stop(paste0("The provided seed_names, ", paste0(seed_names, collapse = ","), ", does not exist among the rownames of expression matrix.\n"))
    } 
  }else{
    seed_names <- "NULLSEED"   
  }
  
  # If target_match variable is not specified, use "positive" as default.
  if(length(target_match) == 0 || nchar(target_match) == 0){
    warning("The target_match variable is not specified. Using 'positive' by default ..\n")
    target_match <- "positive"
  }else if(length(target_match) == 1 && !target_match %in% c("positive", "negative")){
    stop(paste0(target_match, collapse=", "), " is not a valid target_match value. The target_match variable must be 'positive' or 'negative'.")
  }else if(length(target_match) > 1 && all(!target_match %in% c("positive", "negative"))){
    stop(paste0(target_match, collapse=", "), " is not a valid target_match value. The target_match variable must be 'positive' or 'negative'.")
  }else if(length(target_match) > 1 && any(target_match %in% c("positive", "negative"))){
    target_match <- target_match[which(target_match %in% c("positive", "negative"))][1]
    warning("More than one target_match values were specified. Only the first valid target_match value, '", target_match, "', is used.\n")
  }
  
  # If assoc_metric variable is not specified, use "IC" as default.
  if(length(assoc_metric) == 0 || nchar(assoc_metric) == 0){
    warning("The assoc_metric variable is not specified. Using 'IC' by default ..\n")
    assoc_metric <- "IC"
  }else if(length(assoc_metric) == 1 && !assoc_metric %in% c("IC", "COR")){
    stop(paste0(assoc_metric, collapse=", "), " is not a valid assoc_metric value. The assoc_metric variable must be 'IC' or 'COR'.")
  }else if(length(assoc_metric) > 1 && all(!assoc_metric %in% c("IC", "COR"))){
    stop(paste0(assoc_metric, collapse=", "), " is not a valid assoc_metric value. The assoc_metric variable must be 'IC' or 'COR'.")
  }else if(length(assoc_metric) > 1 && any(assoc_metric %in% c("IC", "COR"))){
    assoc_metric <- assoc_metric[which(assoc_metric %in% c("IC", "COR"))][1]
    warning("More than one assoc_metric values were specified. Only the first valid assoc_metric value, '", assoc_metric, "', is used.\n")
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
  if(nrow(mat) < 2){
    warning("You are computing a row-wise statistic over a matrix with nrow < 2.\n")
  }
  
  # Define seed from given seed_names
  if (seed_names == "NULLSEED") {
    seed <- as.vector(rep(0, ncol(mat)))      
  } else {
    if (length(seed_names) > 1) {
      seed <- apply(mat[seed_names,], MARGIN=2, FUN=max)
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
  
  # Only score value from revealer is returned
  colnames(cmi) <- c("score")
  rownames(cmi) <- rownames(mat)
  
  return(cmi)

}

