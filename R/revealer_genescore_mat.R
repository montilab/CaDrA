
#' Row-wise matrix conditional mutual information I(x, y | z) from REVEALER
#' 
#' Compute conditional mutual information scores of x and y given z for each row of a given binary matrix
#' @param mat row matrix of binary features to compute row-wise scores for based on the REVEALER MUTUALLY EXCLUSIVE test
#' @param target a matrix of continuous functional response of interest 
#' @param target_match Direction of the match (negative or positive). Use "positive" to match the higher values of the target, "negative" to match the lower values. Default is positive. 
#' @param seed_names Starting seed, one or more binary features(s) representing known causes of activation or features associated with the target
#' @param seed_combination_op Operation to consolidate and summarize seeds to one vector of values. "max" is default
#' @param exclude_features Features to exclude in the search iterations
#' @param normalize_features Feature row normalization: FALSE (default) or "standardize" or "0.1.rescaling"
#' @param count_thres_low Filter out features with less than count.thres.low. NULL by default.
#' @param count_thres_high Filter out features with more than count.thres.high. NULL by default.
#' @param assoc_metric Assocication Metric: "IC" information coeff. (default) or "COR" correlation.
#' @param verbose a logical indicating whether or not to verbose diagnostic messages. Default is FALSE. 
#'
#' @return A data frame
#' @export
#' @importFrom purrr map_dfr
revealer_genescore_mat <- function
(
  mat,                                   
  target,      
  target_match = "positive",             
  seed_names = NULL,
  seed_combination_op = "max", 
  exclude_features = NULL,
  normalize_features = FALSE,
  count_thres_low = NULL,                   
  count_thres_high = NULL,
  assoc_metric = "IC",
  verbose = FALSE
)
{
  
  # Setup verbose option definition
  options(verbose=verbose)
  
  # Check if the target_score has the same length as the number of samples in the feature matrix
  if(length(target) < ncol(mat))
    stop("'The provided target must have the same length as the number of samples in the feature data matrix.\n")
  
  # If target_match variable is not specified, use "positive" as default.
  if(length(target_match) == 0 || nchar(target_match) == 0){
    warning("The target_match variable is not specified. Using 'positive' by default ..\n")
    target_match <- "positive"
  }else if(length(target_match) == 1 && !target_match %in% c("positive", "negative")){
    warning(paste0(target_match, collapse=", "), " is not a valid target_match value. The target_match variable must be 'positive' or 'positive'. Using 'positive' by default.\n")
    target_match <- "positive"    
  }else if(length(target_match) > 1 && all(!target_match %in% c("positive", "negative"))){
    warning(paste0(target_match, collapse=", "), " is not a valid target_match value. The target_match variable must be 'positive' or 'positive'. Using 'positive' by default.\n")
    target_match <- "positive"
  }else if(length(target_match) > 1 && any(target_match %in% c("positive", "negative"))){
    target_match <- target_match[which(target_match %in% c("positive", "negative"))][1]
    warning("More than one target_match values were specified. Only the first valid target_match value, '", target_match, "', is used.\n")
  }
  
  # Check if seed_names is provided
  if (is.null(seed_names)) seed_names <- "NULLSEED"
  
  # Exclude samples with target == NA
  verbose(paste("Initial target length:", length(target)))      
  locs <- seq(1, length(target))[!is.na(target)]
  target <- target[locs]
  mat <- mat[,locs]
  verbose(paste("Target length after excluding NAs:", length(target)))    
  
  ## Reordering by the direction of matched target
  if (target_match == "negative") {
    ind <- order(target, decreasing=F)
  } else {
    ind <- order(target, decreasing=T)
  }
  
  target <- target[ind]
  mat <- mat[,ind]
  
  # Eliminate flat, sparse or features that are too dense
  if (!is.null(count_thres_low) && !is.null(count_thres_high)) {
    
    sum_rows <- rowSums(mat)
    seed_flag <- rep(0, nrow(mat))
    
    if (seed_names != "NULLSEED") {
      locs <- match(seed_names, row.names(mat))
      locs <- locs[!is.na(locs)]
      seed_flag[locs] <- 1
    }
    
    retain <- rep(0, nrow(mat))
    
    for (i in 1:nrow(mat)) {
      #i=1
      if ((sum_rows[i] >= count_thres_low) && (sum_rows[i] <= count_thres_high)) retain[i] <- 1
      if (seed_flag[i] == 1) retain[i] <- 1
    }
    
    mat <- mat[retain == 1,]
    verbose(paste("Number of features kept:", sum(retain), "(", signif(100*sum(retain)/length(retain), 3), " percent)"))
    
  }
  
  # Normalize features FALSE (DEFAULT), standardized, or 0.1.rescaling
  if (normalize_features == "standardized") {
    for (i in 1:nrow(mat)) {
      mean_row <- mean(mat[i,])
      sd_row <- ifelse(sd(mat[i,]) == 0, 0.1*mean_row, sd(mat[i,]))
      mat[i,] <- (mat[i,] - mean_row)/sd_row
    }
  } else if (normalize_features == "0.1.rescaling") {
    for (i in 1:nrow(mat)) {
      max_row <- max(mat[i,])
      min_row <- min(mat[i,])
      range_row <- ifelse(max_row == min_row, 1, max_row - min_row)
      mat[i,] <- (mat[i,] - min_row)/range_row
    }
  }
  
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
  
  # Exclude user-specified features 
  if (!is.null(exclude_features)) {
    locs <- match(exclude_features, row.names(mat))
    locs <- locs[!is.na(locs)]
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

