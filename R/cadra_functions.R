
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
#' prevalence of a feature across all samples in the dataset above which the
#' feature will be filtered out. Default is 0.6 (feature that occur in
#' 60 percent or more of the samples will be removed)
#' @param min_cutoff a numeric value between 0 and 1 describing the absolute
#' prevalence of a feature across all samples in the dataset below which the
#' feature will be filtered out. Default is 0.03 (feature that occur in
#' 3 percent or less of the samples will be removed)
#' @return An SummarizedExperiment object with only the filtered-in features
#' given the filter thresholds specified
#' @param warning a logical value indicates whether or not to print the
#' diagnostic messages. Default is \code{FALSE}.
#' @examples
#'
#' # Load pre-computed feature set
#' data(sim_FS)
#'
#' # Filter out features having < 3 and > 60% prevalence across all samples
#' # (default)
#' sim_FS_filt1 <- prefilter_data(FS = sim_FS)
#'
#' # Change the min cut-off to 1% prevalence, instead of the default 3%
#' sim_FS_filt2 <- prefilter_data(FS = sim_FS, min_cutoff  = 0.01)
#'
#' # Change the max cut-off to 65% prevalence, instead of the default 60%
#' sim_FS_filt3 <- prefilter_data(FS = sim_FS, max_cutoff = 0.65)
#'
#' @export
#' @import SummarizedExperiment
prefilter_data <- function(
    FS,
    max_cutoff = 0.6,
    min_cutoff = 0.03,
    warning = FALSE
){

  # Set up verbose option
  options(verbose = warning)
  
  # Check if FS is a SummarizedExperiment class object
  if(!is(FS, "SummarizedExperiment"))
    stop("'FS' must be SummarizedExperiment class object
         from SummarizedExperiment package")
  
  # Compute the frequency of feature occurrence across all samples
  # (i.e. fraction of samples having the feature)
  frac <- round(rowSums(assay(FS))/ncol(FS), 2)

  verbose("Pre-filtering features...\n\n")
  verbose("Removing features having < ", min_cutoff*100, "and > ",
          max_cutoff*100, " % occurence in sample set...\n")

  FS <- FS[ (frac >= min_cutoff) & (frac <= max_cutoff) , ]

  verbose(nrow(FS)," features retained out of ", length(frac),
          " supplied features in dataset\n\n")

  return(FS)

}


#' Checks if feature set and input scores are valid dataset
#'
#' @param FS a matrix of binary features or a SummarizedExperiment class object 
#' from SummarizedExperiment package where rows represent features of interest 
#' (e.g. genes, transcripts, exons, etc...) and columns represent the samples. 
#' The assay of FS contains binary (1/0) values indicating the presence/absence 
#' of ‘omics’ features.
#' @param input_score a vector of continuous scores of a molecular phenotype of
#' interest such as protein expression, pathway activity, etc.
#' The \code{input_score} object must have names or labels that match the column
#' names of FS object.
#' @param do_check a logical value indicates whether or not to validate if the  
#' given parameters (FS and input_score) are valid inputs. 
#' Default is \code{TRUE}
#' 
#' @noRd
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
  if(!is(FS, "SummarizedExperiment") & !is(FS, "matrix"))
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
  
  # Check input_score is provided and is a continuous values with no NAs
  if(length(input_score) == 0 ||
     any(!is.numeric(input_score)) || any(is.na(input_score)))
    stop("input_score must contain a vector of continuous values ",
         "(with no NAs) where its vector names match the colnames of ",
         "the FS object.\n")

  # Make sure the input_score has names or labels that are the
  # same as the column names of FS object
  # No need for this check - this will be checked when we compare names with
  # the column names of mat (see below)
  # if(is.null(names(input_score)))
  #   stop("The input_score object must have names or labels to track ",
  #        "the samples by. Please provide unique sample names or labels ",
  #        "that match the column names of the FS object.\n")
  
  # Make sure the input_score has the same length as number of samples in FS
  # We do not need this check either as it will fail if the following if is not TRUE
  # if(length(input_score) != ncol(mat))
  #   stop("The input_score must have the same length ",
  #        "as the number of columns in FS.\n")
  
  if(any(!names(input_score) == colnames(mat)))
      stop("The input_score object must have names or labels that ",
      "match the exact order of the column names of the FS object.")

  # Check if the features have either all 0s or 1s values
  if(any(rowSums(mat) %in% c(0, ncol(mat)) ))
    stop("The FS object has features that are either all 0s or 1s. ",
         "These features must be removed from the FS object as ",
         "they are uninformative.")
  
  if(is.na(max_size) || max_size <= 0 || max_size > nrow(FS)){
    stop("Please specify an integer value specifies a maximum size ",
         "that a meta-feature can extend to do for a given search ",
         "(max_size must be >= 1)",
         "and max_size must be lesser than the number of features in",
         "FS\n")
  }
  
}


#' Check top_N value for candidate_search
#'
#' Checks top_N value is valid
#'
#' @param top_N an integer specifies the number of features to start the
#' search over, starting from the top 'N' features in each case. If \code{top_N}
#' is provided, then \code{search_start} parameter will be ignored. 
#' @param search_start a list of character strings (separated by commas)
#' which specifies feature names within the expression set object to start
#' the search with. If \code{search_start} is provided, then \code{top_N}
#' parameter will be ignored. Default is \code{NULL}.
#' @param rownames row names of input SummarizedExperiment object
#' 
#' @noRd
#'
#' @return an integer value - start search index

check_top_N <- function(
    top_N, search_start, rownames
){
  
  # Check if search_start is given
  if(is.null(search_start)){
    
    if(is.na(top_N) || length(top_N)==0 || top_N <= 0){
      stop("Please specify a NUMERIC top_N value to evaluate over top N ",
           "features (top_N must be >= 1).\n")
    }
    
    if(top_N > length(rownames))
      stop("Please specify a top_N value that is less than the number of ",
           "features in the FS.\n")
    
    if(top_N > 10)
      warning("top_N value specified is greater than 10. ",
              "This may result in a longer search time.\n")
    
    # Start the search with top N features based on their sorted indexes
    search_feature_index <- seq_len(top_N)
    
    verbose("Evaluating search over top ", length(search_feature_index),
            " features\n")
    
  } else {
    
    search_start <- strsplit(as.character(search_start), ",", fixed=TRUE) |>
      unlist() |>
      trimws()
    
    if(!is.na(top_N) && length(top_N) > 0){
      warning("Since start_search variable is given, ",
              "evaluating over top_N value will be ignored.\n")
    }
    
    # User-specified feature name
    # (has to be a character from rownames(1:nrow(FS)))
    verbose("Starting with specified feature names...\n")
    
    if(length(search_start) == 0 || any(!search_start %in% rownames))
      stop("Provided starting feature does not exist among FS's rownames.\n\n")
    
    # Get the index of the search_start strings and start
    # the search with the defined indexes
    search_feature_index <- which(rownames %in% search_start)
    
  } # end else (!is.null)
  
  
  search_feature_index
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
#' Produces a random permutation score matrix given a vector of sample-specific scores
#' representing a phenotypic readout of interest such as protein expression,
#' pathway activity, etc.
#' @param input_score a vector of continuous scores of a molecular phenotype of interest
#' such as protein expression, pathway activity, etc.
#' @param n_perm a number of permutations to generate. This determines
#' the number of rows in the permutation matrix.
#' @return a matrix of ordered values where each row contains the order of
#' a single permuted \code{input_score}.
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
              length(input_score) > 0 & all(is.numeric(input_score)) & 
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

  stopifnot("Permutations are not unique. Try a different n_perm value." =
              (nrow(perm) == nrow(unique.matrix(perm))))
  
  return(perm)

}

