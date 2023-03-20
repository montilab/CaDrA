
# Function to whether print diagnostic messages or not
verbose <- function(...){

  # Fetch verbose option set in candidate_search() function
  opt <- getOption("verbose", FALSE)
  if(!opt) return(invisible(NULL))
  msgs <- list(...)
  message(msgs,"\n")

}

#' Pre-filter features
#'
#' Pre-filter a dataset prior running \code{candidate_search()} to avoid
#' testing features that are too prevalent or too sparse across samples in
#' the dataset
#' @param FS a SummarizedExperiment object containing binary features where
#' rows represent features of interest (e.g. genes, transcripts, exons, etc.)
#' and columns represent the samples.
#' @param max_cutoff a numeric value between 0 and 1 describing the absolute
#' prevalence of a feature across all samples in the FS object which the
#' feature will be filtered out. Default is 0.6 (feature that occur in
#' 60 percent or more of the samples will be removed)
#' @param min_cutoff a numeric value between 0 and 1 describing the absolute
#' prevalence of a feature across all samples in the FS object which the
#' feature will be filtered out. Default is 0.03 (feature that occur in
#' 3 percent or less of the samples will be removed)
#' @return A SummarizedExperiment object with only the filtered-in features
#' given the filtered thresholds
#' @param verbose a logical value indicates whether or not to print the
#' diagnostic messages. Default is \code{FALSE}.
#' @examples
#'
#' # Load pre-computed feature set
#' data(sim_FS)
#'
#' # Filter out features having < 3% and > 60% prevalence across all samples
#' # by (default)
#' sim_FS_filt1 <- prefilter_data(FS = sim_FS)
#'
#' # Change the min cut-off to 1% prevalence, instead of the default of 3%
#' sim_FS_filt2 <- prefilter_data(FS = sim_FS, min_cutoff  = 0.01)
#'
#' # Change the max cut-off to 65% prevalence, instead of the default of 60%
#' sim_FS_filt3 <- prefilter_data(FS = sim_FS, max_cutoff = 0.65)
#'
#' @export
#' @import SummarizedExperiment
prefilter_data <- function(
    FS,
    max_cutoff = 0.6,
    min_cutoff = 0.03,
    verbose = FALSE
){

  # Set up verbose option
  options(verbose = verbose)
  
  # Check if FS is a SummarizedExperiment class object
  if(!is(FS, "SummarizedExperiment"))
    stop("'FS' must be SummarizedExperiment class object
         from SummarizedExperiment package")
  
  # Compute the frequency of feature occurrence across all samples
  # (i.e. fraction of samples having the feature)
  frac <- round(rowSums(SummarizedExperiment::assay(FS))/ncol(FS), 2)

  verbose("Pre-filtering features...\n")
  verbose("Removing features having < ", min_cutoff*100, "% and > ",
          max_cutoff*100, "% occurence in sample set\n")

  FS <- FS[ (frac >= min_cutoff) & (frac <= max_cutoff) , ]

  verbose(nrow(FS), " features retained out of ", length(frac),
          " supplied features in the FS object\n")

  return(FS)

}


#' Checks if feature set and input scores are valid dataset
#'
#' @param FS a matrix of binary features or a SummarizedExperiment class object 
#' from SummarizedExperiment package where rows represent features of interest 
#' (e.g. genes, transcripts, exons, etc...) and columns represent the samples. 
#' The assay of FS contains binary (1/0) values indicating the presence/absence 
#' of omics features.
#' @param input_score a vector of continuous scores of a molecular phenotype of
#' interest such as protein expression, pathway activity, etc.
#' NOTE: The \code{input_score} object must have names or labels that match the column
#' names of FS object.
#' @param do_check a logical value indicates whether or not to validate if the  
#' given parameters (FS and input_score) are valid inputs. Default is \code{TRUE}
#' 
#' @noRd
#' 
#' @examples 
#' 
#' # Load pre-computed feature set
#' data(sim_FS)
#'
#' # Load pre-computed input-score
#' data(sim_Scores)
#'
#' check_data_input(
#'  FS = sim_FS,
#'  input_score = sim_Scores
#' )
#' 
#' @return a filtered feature set and input scores with overlapping samples
#' @import SummarizedExperiment
check_data_input <- function(
    FS,
    input_score,
    do_check = TRUE
){
  
  if(do_check == FALSE) return(NULL)
  
  # Check if FS is a matrix or a SummarizedExperiment class object
  if(!is(FS, "SummarizedExperiment") && !is(FS, "matrix"))
    stop("'FS' must be a matrix or a SummarizedExperiment class object
         from SummarizedExperiment package")
  
  # Retrieve the feature binary matrix
  if(is(FS, "SummarizedExperiment")){
    mat <- as.matrix(SummarizedExperiment::assay(FS))
  }else{
    mat <- FS
  }
  
  # Check if the matrix has only binary 0 or 1 values
  if(length(mat) == 0 || any(!mat %in% c(0,1)) || any(is.na(mat)))
    stop("FS object must contain binary values 0s or 1s (no empty values).")
  
  # Make sure the FS object has row names for features tracking
  if(is.null(rownames(mat)))
    stop("The FS object does not have row names to ",
         "track features by. Please provide unique features or row names ",
         "for the FS object.\n")
  
  # Make sure the FS object has row names for features tracking
  if(is.null(colnames(mat)))
    stop("The FS object does not have column names to ",
         "track samples by. Please provide unique sample names ",
         "for the FS object.\n")
  
  # Check if input_score is provided and contains a vector of continuous values 
  # (with no NAs). Additionally, check if it has names or labels that match 
  # the column names of the FS object
  if(length(input_score) == 0 || any(!is.numeric(input_score)) || 
     any(is.na(input_score)) || is.null(names(input_score)) ||
     any(!names(input_score) %in% colnames(mat)))
    stop("input_score must contain a vector of continuous scores ",
         "(with no NAs), and its vector names or labels must match the column ",
         "names of the FS object.\n")
  
  # Check if the features have either all 0s or 1s values
  if(any(rowSums(mat) %in% c(0, ncol(mat)) ))
    stop("The FS object has features that are either all 0s or 1s. ",
         "These features must be removed from the FS object as ",
         "they are uninformative.")
  
}


#' Check top_N value for candidate_search()
#'
#' Checks top_N value is valid and return a list of indices of top N features 
#' with the best scores to start candidate_search()
#'
#' @param rowscore a vector of directional scores computed based on a
#' given scoring method, feature set (FS), and input score (input_score) using
#' calc_rowscore() function. The scores are ordered from most significant 
#' to least significant.
#' @param feature_names a list of feature names (or row names) in FS object
#' @param top_N an integer specifies the number of features to start the
#' search over (e.g. starting from the top 'N' features with the best scores). 
#' Default is 1.
#' If \code{top_N} is provided, then \code{search_start} parameter will be ignored. 
#' @param search_start a list of character strings (separated by commas)
#' which specify a list of feature names to start the search with. 
#' If \code{search_start} is provided, then \code{top_N} parameter will be ignored. 
#' Default is \code{NULL}.
#' 
#' @noRd
#' 
#' @examples 
#' 
#' # Load library
#' library(SummarizedExperiment)
#' 
#' # Load pre-computed feature set
#' data(sim_FS)
#'
#' # Load pre-computed input-score
#' data(sim_Scores)
#'
#' rowscore <- calc_rowscore(
#'  FS_mat = SummarizedExperiment::assay(sim_FS),
#'  input_score = sim_Scores,
#'  method = "ks_pval",
#'  alternative = "less",
#'  weight = NULL
#' )
#' 
#' top_N_index <- check_top_N(
#'  rowscore = rowscore,
#'  top_N = 7,
#'  search_start = NULL,
#'  feature_names = rownames(sim_FS)
#' )
#'
#' @return If top_N is given, a vector of indices of top N features with the best
#' scores will be returned and used to start the candidate_search() with. 
#' Otherwise, the candidate_search() will start the search with a list of indices 
#' of features defined in search_start.
check_top_N <- function(
    rowscore, 
    feature_names,
    top_N = 1, 
    search_start = NULL
){
  
  # Check if search_start is given
  if(is.null(search_start)){
    
    if(is.na(top_N) || length(top_N) == 0 || top_N <= 0)
      stop("Please specify a NUMERIC top_N value to evaluate over top N ",
           "features (top_N must be >= 1).\n")
    
    if(top_N > length(feature_names))
      stop("Please specify a top_N value that is less than the number of ",
           "features in the FS object.\n")
    
    if(top_N > 10)
      warning("top_N value specified is greater than 10. ",
              "This may result in a longer search time.\n")
    
    # Getting the top N features with the biggest scores
    top_features <- names(rowscore[seq_len(top_N)])
    
    # Retrieve the indexes of features with the best scores in FS object
    search_feature_index <- lapply(seq_along(top_features), function(f){
      #f=1;
      which(feature_names == top_features[f])
    }) |> unlist()
    
    verbose("Evaluating search over top ", length(search_feature_index), " features\n")
    
  }else{
    
    search_start <- strsplit(as.character(search_start), ",", fixed=TRUE) |>
      unlist() |>
      trimws()
    
    if(!is.na(top_N) && length(top_N) > 0)
      warning("Since search_start variable is given, ",
              "evaluating over top_N value will be ignored.\n")
    
    # User-specified feature names
    verbose("Starting with specified feature names...\n")
    
    if(length(search_start) == 0 || any(!search_start %in% feature_names))
      stop("The provided starting features: ", 
           paste0(search_start[which(!search_start %in% feature_names)], collapse=", "),
           " does not exist among the row names of FS object.\n")
    
    # Get the index of the search_start strings
    search_feature_index <- lapply(seq_along(search_start), function(f){
      #f=1;
      which(feature_names == search_start[f])
    }) |> unlist()
    
  }
  
  return(search_feature_index)
  
}



#' ks_test_d_wrap_ wrapper
#'
#' Compute directional Kolmogorov-Smirnov scores
#' @param n_x length of ranked list
#' @param y positions of geneset items in ranked list (ranks)
#' @param alt alternative hypothesis for p-value calculation
#' (\code{"two.sided"} or \code{"greater"} or \code{"less"}).
#' Default is \code{less} for left-skewed significance testing.
#' 
#' @noRd
#' @useDynLib CaDrA ks_test_d_wrap_
#'
#' @return need a return value here
ks_test_double_wrap <- function(n_x, y, alt=c("less", "greater", "two.sided")) {


  if(length(alt) > 0){
    alt_int<- switch(alt, two.sided=0L, less=1L, greater=-1L, 1L)
  } else {
    alt_int <- 1L
  }

  # If input is an empty vector
  if( n_x != length(y) | length(y) < 1) return (NULL)

  res <- .Call(ks_test_d_wrap_,  as.integer(n_x), as.numeric(y), alt_int)
  res

}


#' Random permutation matrix generator
#'
#' Produces a random permutation score matrix given a vector of sample-specific
#' scores representing a phenotypic readout of interest such as protein expression,
#' pathway activity, etc.
#' @param input_score a vector of continuous scores of a molecular phenotype of 
#' interest such as protein expression, pathway activity, etc.
#' NOTE: The \code{input_score} object must have names or labels to track 
#' samples by.
#' @param n_perm a number of permutations to generate. This determines
#' the number of rows in the permutation matrix.
#' @return a matrix of values where each row contains scores of a single 
#' permuted \code{input_score}.
#'
#' @examples
#'
#' # Load pre-simulated scores
#' data(sim_Scores)
#'
#' # Set seed for permutation
#' set.seed(123)
#'
#' # Define number of permutations
#' n_perm = 1000
#'
#' # Generate permuted scores
#' perm_matrix <- generate_permutations(
#'   input_score = sim_Scores,
#'   n_perm = n_perm
#' )
#'
#' @export
generate_permutations <- function(
    input_score,
    n_perm
){
  
  # Check input_score is provided and is a continuous values with no NAs
  stopifnot("input_score must contain a vector of continuous values (no NAs)"= 
              length(input_score) > 0 && all(is.numeric(input_score)) && 
              all(!is.na(input_score)))
  
  # Make sure the input_score has names or labels to track samples by
  stopifnot("input_score object must have names or labels to track samples by"=
              !is.null(names(input_score)))
  
  # Get number of samples
  n <- length(input_score)
  
  # Create permutation matrix
  perm <- matrix(NA, nrow=n_perm, ncol=n)
  colnames(perm) <- names(input_score)

  # Sample the input scores
  for(i in seq_len(n_perm)){
    perm[i,] <- sample(input_score, n, replace=FALSE)
  }

  # Make sure permutations are unique
  stopifnot("Permutations are not unique. Try a different n_perm value." =
              (nrow(perm) == nrow(unique.matrix(perm))))
  
  return(perm)

}

