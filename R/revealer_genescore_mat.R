
#' Row-wise matrix conditional mutual information I(x, y | z) from REVEALER
#' 
#' Compute conditional mutual information scores of x and y given z for each row of a given binary matrix
#' @param mat matrix of binary features to compute row-wise scores for based on the REVEALER MUTUALLY EXCLUSIVE test
#' @param target a matrix of continuous functional response of interest 
#' @param target_match Direction of the match (negative or positive). Use "positive" to match the higher values of the target, "negative" to match the lower values. Default is positive. 
#' @param seed_names Starting seed, one or more binary features(s) representing known “causes” of activation or features associated with the target
#' @param seed_combination_op Operation to consolidate and summarize seeds to one vector of values. "max" is default
#' @param exclude_features Features to exclude in the search iterations
#' @param normalize.features Feature row normalization: FALSE (default) or "standardize" or "0.1.rescaling"
#' @param consolidate_identical_features Consolidate identical features: FALSE (default) or "identical" or "similar" 
#' @param cons_features_hamming_thres If consolidate.identical.features = "similar" then consolidate features within this Hamming dist. thres.
#' @param save_preprocessed_features_dataset save preprocessed features dataset. NULL is default    
#' @param max_n_iter Maximum number of iterations to perform. max_n_iter = 5 by default.
#' @param n_markers Top n hits to display at each iteration. n_markers = 30 by default.
#' @param count_thres_low Filter out features with less than count.thres.low 1's. count_thres_low = 5 by default.
#' @param count_thres_high Filter out features with more than count.thres.low 1's. count_thres_high = 30 by default.
#' @param n_perm Number of permutations (x number of genes) for computing p-vals and FRDs. n_perm = 2 by default.
#' @param assoc_metric Assoc. Metric: "IC" information coeff. (default) or "COR" correlation.
#' @param r_seed Random number generation seed (34578 is default)
#' @return 
#' @export
revealer_genescore_mat <- function
(
  mat,                                   
  target,      
  target_match = "positive",             
  seed_names = NULL,
  seed_combination_op = "max", 
  exclude_features = NULL,
  normalize_features = FALSE,
  consolidate_identical_features = FALSE,
  cons_features_hamming_thres = NULL,
  save_preprocessed_features_dataset = NULL,
  max_n_iter = 5,            
  n_markers = 30,
  count_thres_low = NULL,                   
  count_thres_high = NULL,
  n_perm = 10,
  assoc_metric = "IC",
  r_seed = 34578,
  verbose = TRUE
)
{
  
  # ## Required R packages
  # library(NMF)
  # library(MASS)
  # library(ppcor)
  # library(misc3d)
  # library(smacof)
  # library(maptools)
  # library(RColorBrewer)
  # library(CaDrA)
  # library(Biobase)
  # library(tidyverse)
  # 
  # ## Source the revealer genescore method
  # source("/Users/reinachau/Documents/CaDrA/CaDrA/R/revealer_genescore.R")
  # 
  # ## Read in the simulated dataset from CaDrA
  # data("sim.ES")
  # mat <- exprs(sim.ES)
  #
  # ## Define the parameters for testing
  # target = rnorm(n=ncol(mat), mean=0, sd=2)    
  # target_match = "positive"             
  # seed_names = NULL
  # seed_combination_op = "max"
  # exclude_features = NULL
  # normalize_features = FALSE
  # consolidate_identical_features = FALSE
  # cons_features_hamming_thres = NULL
  # save_preprocessed_features_dataset = "/Users/reinachau/Documents/CaDrA/CaDrA/R/save_preprocessed_features_dataset.csv"
  # max_n_iter = 2            
  # n_markers = 30
  # count_thres_low = 3                   
  # count_thres_high = 50
  # n_perm = 10
  # assoc_metric = "IC"
  # r_seed = 34578
  # verbose = TRUE
  
  # Check if the target_score has the same length as the number of samples in the feature matrix
  if(length(target) < ncol(mat))
    stop("'The provided target must have the same length as the number of samples in the feature data matrix.\n")
  
  # Check if feature matrix contains only binary values
  feature_values <- lapply(mat, 1, unique) %>% unlist()
  
  if(any(!feature_values %in% c(0,1))
     stop("Feature matrix must contain only binary values (0 or 1)")
     
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
  
  ## Set R seeds
  set.seed(r_seed)
  
  # Exclude samples with target == NA
  print(paste("Initial target length:", length(target)))      
  locs <- seq(1, length(target))[!is.na(target)]
  target <- target[locs]
  mat <- mat[,locs]
  print(paste("Target length after excluding NAs:", length(target)))    
  
  ## reordering by the direction of matched target
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
    print(paste("Number of features kept:", sum(retain), "(", signif(100*sum(retain)/length(retain), 3), " percent)"))
    
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
    
    print("Location(s) of seed(s):")
    print(match(seed_names, row.names(mat)))
    
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
  
  # Consolidate identical features
  # This is a very fast way to eliminate perfectly identical features compared with what we do below in "similar"
  if (consolidate_identical_features == "identical") {  
    
    summary_vectors <- apply(mat, MARGIN=1, FUN=paste, collapse="")
    ind <- order(summary_vectors)
    summary_vectors <- summary_vectors[ind]
    mat <- mat[ind,]
    taken <- i_count <- rep(0, length(summary_vectors))
    
    i <- 1
    while (i <= length(summary_vectors)) {
      j <- i + 1
      while ((summary_vectors[i] == summary_vectors[j]) & (j <= length(summary_vectors))) {
        j <- j + 1
      }
      i_count[i] <- j - i
      if (i_count[i] > 1) taken[seq(i + 1, j - 1)] <- 1
      i <- j
    }
    
    if (sum(i_count) != length(summary_vectors)) stop("ERROR")     # Add counts in parenthesis
    row.names(mat) <- paste(row.names(mat), " (", i_count, ")", sep="")
    mat <- mat[taken == 0,]

    # This uses the hamming distance to consolidate similar features up to the Hamming dist. threshold 
    
  } else if (consolidate_identical_features == "similar") { 
    
    hamming_matrix <- hamming.distance(mat)
    taken <- rep(0, nrow(mat))
    
    for (i in 1:nrow(mat)) {
      if (taken[i] == 0) { 
        similar_features <- row.names(mat)[hamming_matrix[i,] <= cons_features_hamming_thres]
        if (length(similar_features) > 1) {
          row.names(mat)[i]  <- paste(row.names(mat)[i], " [", length(similar_features), "]", sep="") # Add counts in brackets
          locs <- match(similar_features, row.names(mat))
          taken[locs] <- 1
          taken[i] <- 0
        }
      }
    }
    
    mat <- mat[taken == 0,]
    
  }
  
  print(paste("Number of features (after filtering and consolidation):", nrow(mat)))
 
  # Save filtered and consolidated file
  if (!is.null(save_preprocessed_features_dataset)) {
    write.csv(data.frame(mat), row.names = row.names(mat), file = save_preprocessed_features_dataset)
  } 
  
  # Compute MI and % explained with original seed(s)
  cmi <- 1:nrow(mat) %>% 
    map_dfr(
      function(r){
        #r=1;
        revealer_genescore(x=target, y=mat[r,], z=seed_vectors, assoc_metric=assoc_metric, target_match=target_match) 
      }
    )
  
  return(cmi)
  
}

