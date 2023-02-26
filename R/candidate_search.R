#'
#' Candidate Search
#'
#' Performs heuristic search on a set of binary features to determine whether
#' there are features whose union is more skewed (enriched at the extremes)
#' than either features alone. This is the main functionality of the \code{CaDrA}
#' package. 
#' 
#' NOTE: The legacy function \code{topn_eval()} is equivalent to the recommended
#' \code{candidate_search()} function
#' @param FS a SummarizedExperiment class object from SummarizedExperiment package
#' where rows represent features of interest (e.g. genes, transcripts, exons, etc.) 
#' and columns represent the samples. The assay of FS contains binary (1/0) values 
#' indicating the presence/absence of omics features.
#' @param input_score a vector of continuous scores representing a phenotypic
#' readout of interest such as protein expression, pathway activity, etc.
#' 
#' NOTE: \code{input_score} object must have names or labels that match the column
#' names of \code{FS} object.
#' @param method a character string specifies a scoring method that is
#' used in the search. There are 6 options: (\code{"ks_pval"} or \code{ks_score}
#' or \code{"wilcox_pval"} or \code{wilcox_score} or 
#' \code{"revealer"} (conditional mutual information from REVEALER) or
#' \code{"custom"} (a customized scoring method)). 
#' Default is \code{ks_pval}.
#' @param custom_function if method is \code{"custom"}, specifies
#' the name of the customized function here. Default is \code{NULL}.
#' 
#' NOTE: custom_function() must take FS_mat (or FS) and input_score as its
#' input arguments, and its final result must return a vector of row-wise scores 
#' ordered from most significant to least significant where its labels or names 
#' matched the row names of FS_mat (or FS) object.
#' @param custom_parameters if method is \code{"custom"}, specifies a list of
#' additional arguments (excluding \code{FS_mat} (or FS) and \code{input_score}) 
#' to be passed to the custom_function(). Default is \code{NULL}.
#' @param alternative a character string specifies an alternative hypothesis
#' testing (\code{"two.sided"} or \code{"greater"} or \code{"less"}).
#' Default is \code{less} for left-skewed significance testing.
#' 
#' NOTE: This argument is applied to KS and Wilcoxon method
#' @param weight if method is \code{ks_score} or \code{ks_pval}, specifying a 
#' vector of weights will perform a weighted-KS testing. Default is \code{NULL}.
#' @param search_start a list of character strings (separated by commas)
#' which specifies feature names within the FS object to start
#' the search with. If \code{search_start} is provided, then \code{top_N}
#' parameter will be ignored. Default is \code{NULL}.
#' @param top_N an integer specifies the number of features to start the
#' search over, starting from the top 'N' features in each case. If \code{top_N}
#' is provided, then \code{search_start} parameter will be ignored. Default is
#' \code{1}.
#' @param search_method a character string specifies an algorithm to filter
#' out the best features (\code{"forward"} or \code{"both"}). Default is
#' \code{both} (i.e. backward and forward).
#' @param max_size an integer specifies a maximum size that a meta-feature
#' can extend to do for a given search. Default is \code{7}.
#' @param best_score_only a logical value indicates whether or not to return 
#' the best score corresponding to each top N searches ONLY.
#' Default is \code{FALSE}.
#' @param do_plot a logical value indicates whether or not to plot the
#' overlapping features of the resulting meta-feature matrix. 
#' 
#' NOTE: plot can only be produced if the resulting meta-feature matrix contains 
#' more than 1 feature (e.g. length(search_start) > 1 or top_N > 1). 
#' Default is \code{FALSE}.
#' @param do_check a logical value indicates whether or not to validate if the  
#' given parameters (FS and input_score) are valid inputs. 
#' Default is \code{TRUE}.
#' @param verbose a logical value indicates whether or not to print the
#' diagnostic messages. Default is \code{FALSE}.
#'
#' @return If \code{best_score_only} is set to \code{TRUE}, the function will
#' return a list of objects containing ONLY the best score of the union 
#' meta-feature matrix for each top N searches. If \code{best_score_only} is set
#' to \code{FALSE}, a list of objects containing the returned meta-feature matrix,
#' as well as its corresponding best score and observed input scores are returned.
#' 
#' @examples
#'
#' # Load pre-computed feature set
#' data(sim_FS)
#'
#' # Load pre-computed input scores
#' data(sim_Scores)
#'
#' # Define additional parameters and run the function
#' candidate_search_result <- candidate_search(
#'   FS = sim_FS, input_score = sim_Scores, 
#'   method = "ks_pval", alternative = "less", weight = NULL, 
#'   search_start = NULL, top_N = 3, search_method = "both",
#'   max_size = 7, best_score_only = FALSE
#' )
#'
#' @export
#' @import SummarizedExperiment
candidate_search <- function(
    FS,
    input_score,
    method = c("ks_pval", "ks_score", "wilcox_pval", "wilcox_score", "revealer", "custom"),
    custom_function = NULL,
    custom_parameters = NULL,
    alternative = c("less", "greater", "two.sided"),
    weight = NULL,
    search_start = NULL,
    top_N = 1,
    search_method = c("both", "forward"),
    max_size = 7,
    best_score_only = FALSE,
    do_plot = FALSE,
    do_check = TRUE,
    verbose = FALSE
){

  # Set up verbose option
  options(verbose = verbose)

  # Match arguments
  method <- match.arg(method)
  alternative <- match.arg(alternative)
  search_method <- match.arg(search_method)

  # Select the appropriate method to compute scores based on
  # skewness of a given binary matrix
  # Return a vector of scores that already been ordered from 
  # most significant to least significant
  # This comes in handy when doing the top-N evaluation of
  # the top N 'best' features
  rowscore <- calc_rowscore(
    FS_mat = SummarizedExperiment::assay(FS),
    input_score = input_score,
    method = method,
    custom_function = custom_function,
    custom_parameters = custom_parameters,   
    alternative = alternative,
    weight = weight,
    seed_names = NULL,
    do_check = do_check,
    verbose = verbose,
    FS = FS,
    search_start = search_start,
    top_N = top_N,
    search_method = search_method,
    max_size = max_size,
    best_score_only = best_score_only,
    do_plot = do_plot
  )
  
  ###### INITIALIZE VARIABLES ###########
  #######################################
  
  # Check if top_N is given and is numeric
  top_N <- as.integer(top_N)
  
  # Get the indices of features to start the search
  search_feature_index <- check_top_N(
    rowscore = rowscore, 
    feature_names = rownames(FS),
    top_N = top_N, 
    search_start = search_start 
  )
  
  ## Check the search_method variable ####
  back_search <- ifelse(search_method == "both", TRUE, FALSE)

  ## Check the max_size variable ####
  max_size <- as.integer(max_size)

  if(is.na(max_size) || length(max_size)==0 || max_size <= 0 || max_size > nrow(FS))
    stop("Please specify a maximum size that a meta-feature can extend to do ",
         "for a given search (max_size must be >= 1)",
         "and max_size must be lesser than the number of features in FS\n")
  
  # Performs the search based on the given indices of the starting best features 
  topn_l <- lapply(seq_along(search_feature_index), function(x){
    # x=1;
    # Extract the index of the feature 
    best_s_index <- search_feature_index[x]
    
    # Get the name of the starting feature and its best score
    start_feature <- rownames(FS)[best_s_index]
    best_feature <- start_feature
    best_s <- rowscore[start_feature]
    
    verbose("Top Feature ", x, ": ", start_feature, "\n")
    verbose("Score: ", best_s, "\n")
    
    # Fetch the vector corresponding to best score
    # Set this as the initial 'meta-feature'
    best_meta <- as.numeric(SummarizedExperiment::assay(FS)[best_s_index,])
    
    # Counter variable for number of iterations
    i <- 0
    
    # Variable to store best score attained over all iterations
    # initialize this to the starting best score
    global_best_s <- best_s

    # Vector of features in the (growing) obtained meta-feature.
    # Begin with just the starting feature
    global_best_s_features <- c()

    ###### BEGIN ITERATIONS ###############
    #######################################

    verbose("\nBeginning candidate search...\n")
    
    while((best_s > global_best_s | i == 0) &&
          (length(global_best_s_features) < max_size)){
      
      verbose("Iteration number: ", (i+1), "\n")
      verbose("Global best score: ", global_best_s, "\n")
      verbose("Previous score: ", best_s, "\n")

      # Update scores and feature set since entry into the
      # loop means there is an improvement (iteration > 0)
      global_best_s <- best_s
      global_best_s_features <- c(global_best_s_features, best_feature)

      verbose("Current feature set: ", global_best_s_features, "\n")

      if(i != 0){

        verbose("Found feature that improves score!\n")

        # Update the new best meta feature (from meta mat)
        best_meta <- meta_mat[hit_best_s_index,]

        # Add that index to the group of indices to be excluded
        # for subsequent checks
        # Here we go off the rownames in the original matrix to
        # find which index to exclude from the FS in subsequent iterations
        best_s_index <- c(best_s_index, which(rownames(FS) == best_feature))

      }

      # Perform a backward check on the list of existing features and
      # update global scores/feature lists accordingly
      if(length(global_best_s_features) > 3 & back_search == TRUE){

        backward_search_results <- forward_backward_check(
          FS = FS,
          input_score = input_score,
          method = method,
          custom_function = custom_function,
          custom_parameters = custom_parameters,
          alternative = alternative,
          weight = weight,
          glob_f = global_best_s_features,     
          glob_f_s = global_best_s,  
          search_start = search_start,
          top_N = top_N,
          search_method = search_method,
          max_size = max_size,
          best_score_only = best_score_only,
          do_plot = do_plot       
        )
        
        # Update globlal features, scores
        global_best_s_features <- backward_search_results[["best_features"]]
        global_best_s <- backward_search_results[["best_scores"]]

        # Update best_meta based on feature set
        best_meta <- as.numeric(ifelse(
          colSums(SummarizedExperiment::assay(FS)[global_best_s_features,]) == 0, 0, 1))

      }

      # Take the OR function between that feature and all
      # other features to see which gives the best score
      # Keep in mind, the number of rows in meta_mat keeps reducing by one
      # each time we find a hit that improves the score
      verbose("Forming meta-feature matrix with all other features in dataset.")

      # Here "*1" is used to convert the boolean back to integer 1's and 0's
      # Notice we remove anything in best_s_index from the original matrix
      # first to form the meta matrix.
      meta_mat <- base::sweep(SummarizedExperiment::assay(FS)[-best_s_index,], 2, best_meta, `|`)*1

      # Check if there are any features that are all 1s generated from
      # taking the union between the matrix
      # We cannot compute statistics for such features and they thus need
      # to be filtered out
      if(any(rowSums(meta_mat) == ncol(meta_mat))){
        warning("Features with all 1s generated from taking the matrix union ",
                "will be removed before progressing...\n")
        meta_mat <- meta_mat[rowSums(meta_mat) != ncol(meta_mat),]
      }

      # With the newly formed 'meta-feature' matrix, 
      # compute row-wise directional scores
      meta_rowscore <- calc_rowscore(
        FS_mat = meta_mat,
        input_score = input_score,
        method = method,
        custom_function = custom_function,
        custom_parameters = custom_parameters,       
        alternative = alternative,
        weight = weight,
        seed_names = NULL,
        do_check = FALSE,
        verbose = FALSE,
        FS = SummarizedExperiment::SummarizedExperiment(assays=meta_mat),
        search_start = search_start,
        top_N = top_N,
        search_method = search_method,
        max_size = max_size,
        best_score_only = best_score_only,
        do_plot = do_plot
      )
      
      # Get the first best score from already ordered meta-scores
      best_s <- meta_rowscore[1]

      # Get feature name of the first best score
      best_feature <- names(best_s)

      # If no improvement (exiting loop)
      if(best_s <= global_best_s){
        verbose("No further improvement in score has been found.")
      }else{
        verbose("Feature that produced best score in combination ",
                "with previous meta-feature: ", best_feature)
        verbose("Score: ", best_s)
        
        # Get the best hit index from meta-feature matrix
        hit_best_s_index <- which(rownames(meta_mat) == best_feature)
      }
      
      # Increment the next loop
      i <- i+1

    } ######### End of while loop

    verbose("\n\nFinished!\n\n")
    verbose("Number of iterations covered: ", i, "\n")
    verbose("Best score attained over ", i , " iterations: ", global_best_s, "\n")

    if(length(global_best_s_features) == 1)
      verbose("No meta-feature that improves the enrichment was found\n")

    verbose("Features returned in FS: ", global_best_s_features, "\n")

    # We don't just want the combination (meta-feature) at the end.
    # We want all the features that make up the meta-feature
    # This can be obtained using the list of indices that were progressively
    # excluded (if at all) in the step-wise procedure
    # If returning only those features that led to the best global score
    FS_best <- FS[global_best_s_features,]

    # Make a list containing two elements.
    #The first will be the FS with the features that gave the best meta-feature
    #The second will be the score corresponding to the meta-feature
    # (named by the starting feature that led to that score)
    # Assign the name of the best meta-feature score to be the
    # starting feature that gave that score
    names(global_best_s) <- start_feature

    return(list("feature_set" = FS_best,
                "input_score" = input_score,
                "score" = global_best_s))
    
  })

  # length of search_feature must be greater 2 in order to generate topn_plot
  if(do_plot) topn_plot(topn_list = topn_l)

  # best_score_only
  if(best_score_only == TRUE){

    # Extract best score from each top N run
    scores_l <- lapply(seq_along(topn_l), function(l){ topn_l[[l]][['score']] })

    # Obtain the best scores with its associated feature names
    best_meta_scores <- unlist(scores_l)

    # Fetch the best score with the highest value
    best_score <- best_meta_scores[order(best_meta_scores, decreasing = TRUE)][1]

    return(best_score)

  }

  # By Default, the function returns the top N candidate search results as a list of lists
  return(topn_l)

}


# Performance backward selection
#' @param FS a SummarizedExperiment class object from SummarizedExperiment package
#' where rows represent features of interest (e.g. genes, transcripts, exons, etc...) 
#' and columns represent the samples. The assay of FS contains binary (1/0) values 
#' indicating the presence/absence of ‘omics’ features.
#' @param input_score a vector of continuous scores representing a phenotypic
#' readout of interest such as protein expression, pathway activity, etc.
#' The \code{input_score} object must have names or labels that match the column
#' names of FS object.
#' @param method a character string specifies a scoring method that is
#' used in the search. There are 4 options: (\code{"ks"} or \code{"wilcox"} or
#' \code{"revealer"} (conditional mutual information from REVEALER) or
#' \code{"custom"} (a user customized scoring method)). Default is \code{ks}.
#' @param custom_function if method is \code{"custom"}, specifies
#' the customized function here. Default is \code{NULL}.
#' @param custom_parameters if method is \code{"custom"}, specifies a list of
#' arguments to be passed to the custom_function(). Default is \code{NULL}.
#' @param alternative a character string specifies an alternative hypothesis
#' testing (\code{"two.sided"} or \code{"greater"} or \code{"less"}).
#' Default is \code{less} for left-skewed significance testing.
#' @param weight if method is \code{ks_pval} or \code{ks_score}, specifying a vector 
#' of weights will perform a weighted-KS testing. Default is \code{NULL}.
#' @param glob_f a vector containing the features (or row names) whose
#' union gives the best score (so far) in the search.
#' Feature names should match those of the provided FS object
#' @param glob_f_s a vector of scores corresponding to the union of the 
#' specified vector of features
#' @param ... additional parameters to be passed to the custom_function() 
#' 
#' @noRd
#' 
#' @return return the set of features with its current best score that is
#' better than the existing meta-feature best score
#' @import SummarizedExperiment
forward_backward_check <- function
(
  FS,
  input_score,
  method,
  custom_function,
  custom_parameters,  
  alternative,
  weight,
  glob_f,
  glob_f_s,
  ...
){

  verbose("Performing backward search...\n")
  verbose("Iterating over ", length(glob_f), " chosen features...\n")

  # Matrix of only global best features so far
  gmat <- SummarizedExperiment::assay(FS[glob_f,])
  rownames(gmat) <- glob_f

  # Here, we make a list that should store the features and their
  # corresponding meta-feature score for each leave-one-out run
  f_names <- list()
  f_scores <- c()

  # We want to see if leaving anyone feature out improves the overall
  # meta-feature score
  # Here we only consider previous features in the meta-feature to remove
  # (i.e. not the last one which was just added)
  for(n in seq_len(length(glob_f)-1)){
    #n=1;
    f_names[[n]] <- glob_f[-n]

    # Take leave-one-out union of features from matrix
    # This will result in a single vector to compute the scores on
    f_union <- ifelse(colSums(gmat[-n,]) == 0, 0, 1)

    # Here we are getting the union of the meta-features
    u_mat <- matrix(f_union, nrow=1, byrow=TRUE,
                    dimnames=list(c("OR"), colnames(FS)))
    
    # Compute scores for this meta feature
    u_rowscore <- calc_rowscore(
      FS_mat = u_mat,
      input_score = input_score,
      method = method,
      custom_function = custom_function,
      custom_parameters = custom_parameters,     
      alternative = alternative,
      weight = weight,
      seed_names = NULL,
      do_check = FALSE,
      verbose = FALSE,
      FS = SummarizedExperiment::SummarizedExperiment(assays=u_mat),
      ...
    )

    # Store score to list
    f_scores <- c(f_scores, u_rowscore[1])

  } # end for loop

  # Obtain index of union meta-feature that gives the best score
  f_best_index <- which.max(f_scores)

  # Obtain best score of the union meta-feature
  f_best_score <- f_scores[f_best_index]

  # Here, we check if the new meta-feature has a computed best score better than 
  # the previous meta-feature's score
  if(f_best_score > glob_f_s){

    verbose("Found improvement on removing existing feature...\n")
    verbose("New feature set: ", f_names[[f_best_index]], "\n")
    verbose("New global best score: ", f_best_score, "\n")

    # Return a set of features that gave a better score than
    # the existing best score and its new best score as well
    return(list(best_features=f_names[[f_best_index]], best_scores=f_best_score))

  }else{

    # Don't change anything. Return the existing
    # best set of features and the corresponding score
    return(list(best_features=glob_f, best_scores=glob_f_s))

  }

}


