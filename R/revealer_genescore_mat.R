

#' Row-wise matrix conditional mutual information I(x, y | z) from REVEALER
#' 
#' Compute conditional mutual information scores of x and y given z for each row of a given binary matrix
#' @param mat row matrix of binary features to compute row-wise scores for based on the REVEALER MUTUALLY EXCLUSIVE test
#' @param target a matrix of continuous functional response of interest 
#' @param target_match Direction of the match (negative or positive). Use "positive" to match the higher values of the target, "negative" to match the lower values. Default is positive. 
#' @param seed_names Starting seed, one or more binary features(s) representing known causes of activation or features associated with the target
#' @param seed_combination_op Operation to consolidate and summarize seeds to one vector of values. "max" is default
#' @param assoc_metric Assocication Metric: "IC" information coeff. (default) or "COR" correlation.
#' @param verbose a logical indicating whether or not to verbose diagnostic messages. Default is TRUE. 
#'
#' @return A data frame with two columns: \code{score} and \code{p_value}
#' @export
#' @importFrom purrr map_dfr
revealer_genescore_mat <- function
(
  mat,                                   
  target,      
  target_match = c("positive", "negative"),             
  seed_names = NULL,
  seed_combination_op = c("min", "max"), 
  assoc_metric = "IC",
  verbose = TRUE
)
{
  
  # Setup verbose option definition
  options(verbose=verbose)
  
  # Check if the matrix has only binary 0 or 1 values 
  if(length(mat) == 0 || !is.matrix(mat) || any(!mat %in% c(0,1)))
    stop("mat variable must be a matrix with binary values only.\n")
  
  # Check if the target_score has the same length as the number of samples in the feature matrix
  if(length(target) != ncol(mat))
    stop("'The provided target must have the same length as the number of columns in the expression matrix.\n")
  
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
  
  # Exclude samples with target == NA
  if(length(which(is.na(target))) > 0){
    verbose(paste0("Excluding NAs from target variable..."))
    locs <- which(!is.na(target))
    target <- target[locs]
    mat <- mat[,locs]
    verbose(paste0("Length of target after removing NAs: ", length(target)))   
  }
  
  # Check if the dataset has any all 0 or 1 features (these are to be removed since they are not informative)
  if(any(rowSums(mat) == 0) || any(rowSums(mat) == ncol(mat))){
    warning("Provided dataset has features that are either all 0 or 1. These features will be removed from the computation.\n")
    mat <- mat[!(rowSums(mat) == 0 | rowSums(mat) == ncol(mat)),]
  }
  
  if(nrow(mat) < 2)
    warning("You are computing a row-wise statistic over a matrix with nrow < 2.\n")
  
  # Check if seed_names is provided
  if (is.null(seed_names)) seed_names <- "NULLSEED"
  
  ## Reordering by the direction of matched target
  if (target_match == "negative") {
    ind <- order(target, decreasing=F)
  } else if (target_match == "positive"){
    ind <- order(target, decreasing=T)
  }
  
  target <- target[ind]
  mat <- mat[,ind]
  
  # Define seeds
  if (seed_names == "NULLSEED") {
    
    seed <- as.vector(rep(0, ncol(mat)))      
    seed_vectors <- as.matrix(t(seed))
    
  } else {
    
    verbose("Location(s) of seed(s):")
    verbose(match(seed_names, row.names(mat)))
    
    if (length(seed_names) > 1) {
      seed <- apply(mat[seed_names,], MARGIN=2, FUN=seed_combination_op)
      seed_vectors <- as.matrix(mat[seed_names,])
    } else {
      seed <- mat[seed_names,]
      seed_vectors <- as.matrix(t(mat[seed_names,]))
    }
    
    locs <- match(seed_names, row.names(mat))
    mat <- mat[-locs,]
    
  }
  
  # Compute MI and % explained with original seed(s)
  cmi <- 1:nrow(mat) %>% 
    purrr::map_dfr(
      function(r){
        #r=1;
        revealer_genescore(x=target, y=mat[r,], z=seed_vectors, assoc_metric=assoc_metric, target_match=target_match) 
      }
    )
  
  return(cmi)
  
}

