
#' \code{REVEALER} Scoring Method
#' 
#' Compute conditional mutual information of \code{x} and \code{y} 
#' given \code{z} for each row of a given binary feature matrix
#' @param mat a matrix of binary features (required)
#' @param input_score a vector of continuous scores of a response of interest 
#' (required). \code{input_score} must have labels or names that associated 
#' with the colnames of the feature matrix. 
#' @param target_match a direction of target matching (\code{"negative"} or 
#' \code{"positive"}). Use \code{"positive"} to match higher values of 
#' \code{input_score}, \code{"negative"} to match lower values of 
#' \code{input_score}. Default is \code{positive}. 
#' @param seed_names one or more features representing known “causes” 
#' of activation or features associated with a response of interest
#' @param assoc_metric an association metric: \code{"IC"} for information 
#' coefficient or \code{"COR"} for correlation. Default is \code{IC}. 
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
#' input_score <- rnorm(n = ncol(sim.ES))
#' names(input_score) <- colnames(sim.ES)
#' 
#' # Define additional parameters and start the function
#' revealer_genescore_result <- revealer_genescore_mat(
#'   mat = exprs(sim.ES), input_score = input_score, 
#'   target_match = "positive", assoc_metric = "IC"
#' )
#'  
#' @export
revealer_genescore_mat <- function
(
  mat,                                   
  input_score, 
  target_match = c("positive", "negative"),             
  seed_names = NULL,
  assoc_metric = c("IC", "COR")
)
{

  assoc_metric <- match.arg(assoc_metric)
  target_match <- match.arg(target_match)
  
  ## Make sure mat variable is a matrix
  mat <- as.matrix(mat)
  
  # If mat has only one column, it must be converted to a row-wise matrix 
  # form as it is needed for backward_forward_search() computation
  # mat must have rownames to track features and columns to track samples
  # for n = 1 case, it is only in backward_forward_search(), thus we can 
  # assign a random labels to it
  if(ncol(mat) == 1){
    mat <- matrix(t(mat), nrow=1, byrow=TRUE, 
                  dimnames = list("my_label", rownames(mat))) 
  }
  
  # Check if the matrix has only binary values and no empty values
  if(length(mat) == 0 || !is.matrix(mat) || 
     any(!mat %in% c(0,1)) || any(is.na(mat)))
    stop("mat variable must be a matrix of binary values (no empty values).")
  
  # Check if input_score is provided and no empty values
  if(length(input_score) == 0 || 
     any(!is.numeric(input_score)) || any(is.na(input_score)))
    stop("input_score variable must be provided and are numeric ",
         "with no empty values.\n")
  
  # Make sure the mat variable has rownames for features tracking
  if(is.null(rownames(mat)))
    stop("The mat object does not have rownames or featureData to ",
         "track the features by. Please provide unique features or ",
         "rownames for the expression matrix.\n")
  
  # Make sure the input_score variable has names as the colnames in mat
  if(is.null(names(input_score)))
    stop("The input_score object must have names or labels to track ",
         "the samples by. Please provide the sample names or labels that ",
         "matches the colnames of the expression matrix.\n")
  
  # Make sure the input_score has the same length as number of samples in mat
  if(length(input_score) != ncol(mat)){
    stop("The input_score variable must have the same length as ",
         "the number of columns in mat.\n")
  }else{
    # check if input_score has any labels or names
    if(length(names(input_score)) == 0){
      stop("The input_score object must have names or labels that ",
           "match the colnames of the expression matrix.\n")
    }
    
    # check if input_score has labels or names that matches the 
    # colnames of the expression matrix
    if(any(!names(input_score) %in% colnames(mat))){
      stop("The input_score object have names or labels that do not match ",
           "the colnames of the expression matrix.\n")
    }
    
    # match colnames of expression matrix with names of 
    # provided input_score values
    # if nrow = 1, if it is, convert to matrix form as it is needed for 
    # backward_forward_search with one dimension matrix computation
    if(nrow(mat) == 1){
      mat <- matrix(t(mat[,names(input_score)]), nrow=1, byrow=TRUE, 
                    dimnames = list(rownames(mat), colnames(mat))) 
    }else{
      mat <- mat[,names(input_score)]
    }
  }
  
  # Check if seed_names is provided
  if(length(seed_names) > 0){
    if(any(!seed_names %in% rownames(mat))){
      stop(paste0("The provided seed_names, ", 
                  paste0(seed_names, collapse = ","), ", 
                  does not exist among the rownames of expression matrix."))
    } 
  }

  # Check if the dataset has any all 0 or 1 features 
  # (these are to be removed since they are not informative)
  if(any(rowSums(mat) == 0) || any(rowSums(mat) == ncol(mat))){
    warning("The provided matrix has some features that are either all 0s or 1s.",
            "These features will be removed from downsteam computation.\n")
    mat <- mat[!(rowSums(mat) == 0 | rowSums(mat) == ncol(mat)),]
  }
  
  # Make sure matrix is not empty after removing uninformative features
  if(nrow(mat) == 0){
    stop("After removing features that are either all 0s or 1s. ",
         "There are no more features remained for downsteam computation.\n")
  }
  
  # Define seed from given seed_names
  if(length(seed_names) == 0){
    seed <- as.vector(rep(0, ncol(mat)))      
  } else {
    if (length(seed_names) > 1) {
      # Consolidate or summarize one or more seeds to one vector of values
      seed <- as.numeric(ifelse(colSums(mat[seed_names,]) == 0, 0, 1))
    } else {
      seed <- mat[seed_names,]
    }
    locs <- match(seed_names, row.names(mat))
    mat <- mat[-locs,]
  }
  
  # Compute MI and % explained with original seed(s)
  cmi <- seq_len(nrow(mat)) %>% 
    purrr::map_dbl(
      function(r){
        revealer_genescore(x=input_score, y=mat[r,], 
                           z=seed, assoc_metric=assoc_metric, 
                           target_match=target_match) 
      }
    )
  
  # Convert list to data.frame
  # Only score value from revealer is returned
  cmi <- data.frame(score=cmi)
  rownames(cmi) <- rownames(mat)

  return(cmi)
  
}

