#'
#' Candidate Search
#'
#' Performs heuristic search on a set of binary features to determine whether
#' there are features whose union is more skewed (enriched at the extremes)
#' than either features alone. This is the main functionality of
#' the \code{CaDrA} package.
#'
#' NOTE: The legacy function \code{topn_eval} is equivalent to the recommended
#' \code{candidate_search} function
#'
#' @param FS a matrix of binary features or a SummarizedExperiment class object
#' from SummarizedExperiment package where rows represent features of interest
#' (e.g. genes, transcripts, exons, etc...) and columns represent the samples.
#' The assay of FS contains binary (1/0) values indicating the presence/absence
#' of omics features.
#' @param input_score a vector of continuous scores representing a phenotypic
#' readout of interest such as protein expression, pathway activity, etc.
#'
#' NOTE: \code{input_score} object must have names or labels that match the
#' column names of \code{FS} object.
#' @param method a character string specifies a scoring method that is
#' used in the search. There are 7 options: (\code{"ks_pval"} or \code{ks_score}
#' or \code{"wilcox_pval"} or \code{wilcox_score} or
#' \code{"revealer"} (conditional mutual information from REVEALER) or
#' \code{"correlation"} (based on simple correlation - pearson or spearman) or
#' \code{"custom"} (a user-defined scoring method)).
#' Default is \code{ks_pval}.
#' @param method_alternative a character string specifies an alternative
#' hypothesis testing (\code{"two.sided"} or \code{"greater"} or \code{"less"}).
#' Default is \code{less} for left-skewed significance testing.
#'
#' NOTE: This argument only applies to \code{ks_pval} and
#' \code{wilcox_pval} method
#' @param custom_function if method is \code{"custom"}, specifies
#' a user-defined function here. Default is \code{NULL}.
#'
#' NOTE: \code{custom_function} must take FS and input_score as its
#' input arguments and its final result must return a vector of row-wise scores
#' where its labels or names match the row names of \code{FS} object.
#' @param custom_parameters if method is \code{"custom"}, specifies a list of
#' additional arguments (excluding \code{FS} and \code{input_score})
#' to be passed to the \code{custom_function}. For example:
#' custom_parameters = list(alternative = "less"). Default is \code{NULL}.
#' @param weights if method is \code{ks_score} or \code{ks_pval}, specifying a
#' vector of weights will perform a weighted-KS testing. Default is \code{NULL}.
#'
#' NOTE: \code{weights} must have names or labels that match the labels of
#' \code{input_score}.
#' @param search_start a vector of character strings (separated by commas)
#' specifies feature names in the \code{FS} object to start the search with.
#' If \code{search_start} is provided, then \code{top_N} parameter will be
#' ignored and vice versa. Default is \code{NULL}.
#' @param top_N an integer specifies the number of features to start the
#' search over. By default, it starts with the feature that has the highest
#' best score (top_N = 1).
#'
#' NOTE: If \code{top_N} is provided, then \code{search_start} parameter
#' will be ignored and vice versa. If top_N > 10, it may result in a longer
#' search time.
#' @param search_method a character string specifies an algorithm to filter
#' out the best features (\code{"forward"} or \code{"both"}). Default is
#' \code{both} (i.e. backward and forward).
#' @param max_size an integer specifies a maximum size that a meta-feature
#' can extend to do for a given search. Default is \code{7}.
#' @param best_score_only a logical value indicates whether or not to return
#' the best score corresponding to each top N searches only.
#' Default is \code{FALSE}.
#' @param do_plot a logical value indicates whether or not to plot the
#' overlapping features of the resulting meta-feature matrix.
#'
#' NOTE: plot can only be produced if the resulting meta-feature matrix contains
#' more than 1 feature (e.g. length(search_start) > 1 or top_N > 1).
#' Default is \code{FALSE}.
#' @param verbose a logical value indicates whether or not to print the
#' diagnostic messages. Default is \code{FALSE}.
#' @param cmethod correlation method to use - spearman or pearson. Default is "spearman"
#' #' NOTE: This argument only applies to \code{correlation} method only
#'
#' @return If \code{best_score_only = TRUE}, the heuristic search will return
#' the best feature whose its union meta-feature matrix has the highest score
#' among the \code{top_N} feature searches.
#' If \code{best_score_only = FALSE}, a list of objects pertaining to
#' \code{top_N} feature searches will be returned. For each top_N feature search,
#' the candidate search will contain 7 objects: (1) its best meta-feature matrix
#' (\code{feature_set}), (2) its observed input scores (\code{input_score}),
#' (3) its corresponding best score pertaining to the union meta-feature
#' matrix (\code{score}), (4) names of the best meta-features (\code{best_features}),
#' (5) rank of the best meta-features in term of their best scores (\code{best indices}),
#' (6) marginal scores of the best meta-features (\code{marginal_best_scores}),
#' (7) cumulative scores of the best meta-features (\code{cumulative_best_scores}).
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
#'   method = "ks_pval", method_alternative = "less", weights = NULL,
#'   search_start = NULL, top_N = 3, search_method = "both",
#'   max_size = 7, best_score_only = FALSE
#' )
#'
#' @export
#' @import SummarizedExperiment
candidate_search <- function(
    FS,
    input_score,
    method = c("ks_pval", "ks_score", "wilcox_pval", "wilcox_score",
               "revealer", "correlation", "custom"),
    method_alternative = c("less", "greater", "two.sided"),
    cmethod = c("spearman", "pearson"),
    custom_function = NULL,
    custom_parameters = NULL,
    weights = NULL,
    search_start = NULL,
    top_N = 1,
    search_method = c("both", "forward"),
    max_size = 7,
    best_score_only = FALSE,
    do_plot = FALSE,
    verbose = FALSE
){

  # Set up verbose option
  options(verbose = verbose)

  # Match arguments
  method <- match.arg(method)
  method_alternative <- match.arg(method_alternative)
  search_method <- match.arg(search_method)

  if (method == "correlation") {
    cmethod <- match.arg(cmethod)
  }

  # Select the appropriate method to compute scores based on
  # skewness of a given binary matrix
  # Return a vector of scores that already been ordered from
  # most significant to least significant
  # This comes in handy when doing the top-N evaluation of
  # the top N 'best' features
  rowscore <- calc_rowscore(
    FS = FS,
    input_score = input_score,
    meta_feature = NULL,
    method = method,
    method_alternative = method_alternative,
    cmethod = cmethod,
    custom_function = custom_function,
    custom_parameters = custom_parameters,
    weights = weights,
    search_start = search_start,
    top_N = top_N,
    search_method = search_method,
    max_size = max_size,
    best_score_only = best_score_only,
    do_plot = do_plot,
    do_check = TRUE,
    verbose = verbose
  )

  # Re-order rowscore in a decreasing order (from highest to lowest values)
  # This comes in handy when doing the top-N evaluation of
  # top N 'best' features
  rowscore <- rowscore[order(rowscore, decreasing=TRUE)]

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

  if(is.na(max_size) || length(max_size)==0 ||
     max_size <= 0 || max_size > nrow(FS))
    stop("Please specify a maximum size that a meta-feature matrix can extend ",
         "to do for a given search (e.g., max_size must be >= 1) ",
         "and max_size must be lesser than the number of features in FS.\n")

  # Performs the search based on the given indices of the starting best features
  topn_l <- lapply(seq_along(search_feature_index), function(x){
    # x=1;
    # Extract the index of the feature
    feature_best_index <- search_feature_index[x]

    # Get the name of the starting feature and its best score
    start_feature <- rownames(FS)[feature_best_index]
    best_feature <- start_feature
    best_score <- rowscore[best_feature]

    verbose("***Top Feature ", x, ": ", best_feature, "***\n")
    verbose("Score: ", best_score, "\n")

    # Update score and feature set since entry into the
    # loop means there is an improvement (i > 0)
    global_best_s <- best_score

    # Initiate list of global features, indices, and scores
    global_best_s_features <- c()
    global_best_s_scores <- c()

    # Counter variable for number of iterations
    i <- 0

    ###### BEGIN ITERATIONS ###############
    #######################################

    verbose("Beginning candidate search...\n")

    while((best_score > global_best_s | i == 0)
          && (length(global_best_s_features) < max_size)
          && (length(global_best_s_features) < nrow(FS))){

      verbose("- Iteration number: ", (i+1), "\n")
      verbose("Global best score: ", global_best_s, "\n")
      verbose("Previous score: ", best_score, "\n")

      # Update global score with current best score
      global_best_s <- best_score

      # Update meta-feature set, its indices, and best scores
      global_best_s_features <- c(global_best_s_features, best_feature)
      global_best_s_scores <- c(global_best_s_scores, best_score)

      verbose("Current meta-feature set:\n ", global_best_s_features, "\n")

      # Perform a backward check on the list of existing features and
      # update global scores/feature lists accordingly
      if(length(global_best_s_features) > 3 & back_search == TRUE){

        backward_search_results <- forward_backward_check(
          FS = FS,
          input_score = input_score,
          method = method,
          method_alternative = method_alternative,
          cmethod = cmethod,
          custom_function = custom_function,
          custom_parameters = custom_parameters,
          weights = weights,
          search_start = search_start,
          top_N = top_N,
          search_method = search_method,
          max_size = max_size,
          best_score_only = best_score_only,
          do_plot = do_plot,
          do_check = FALSE,  # MAKE SURE DO_CHECK IS SILENCE HERE
          verbose = verbose,
          glob_f = global_best_s_features,
          glob_s = global_best_s,
          glob_f_scores = global_best_s_scores
        )

        # Update global features, scores
        global_best_s_features <- backward_search_results[["best_features"]]
        global_best_s_scores <- backward_search_results[["best_scores"]]
        global_best_s <- backward_search_results[["score"]]

        verbose("Current meta-feature set after performing backward search:\n ",
                global_best_s_features, "\n")

      }

      # Compute row-wise directional scores given known meta features
      meta_rowscore <- calc_rowscore(
        FS = FS,
        input_score = input_score,
        meta_feature = global_best_s_features,
        method = method,
        method_alternative = method_alternative,
        cmethod = cmethod,
        custom_function = custom_function,
        custom_parameters = custom_parameters,
        weights = weights,
        search_start = search_start,
        top_N = top_N,
        search_method = search_method,
        max_size = max_size,
        best_score_only = best_score_only,
        do_plot = do_plot,
        do_check = FALSE,  # MAKE SURE DO_CHECK IS SILENCE HERE
        verbose = FALSE    # MAKE SURE VERBOSE IS SILENCE HERE
      )

      # Set up verbose option
      options(verbose = verbose)

      # Get index of highest best score
      best_index <- which.max(meta_rowscore)

      # Get the highest best score
      best_score <- meta_rowscore[best_index]

      # Get feature name of highest best score
      best_feature <- names(best_score)

      # If no improvement (exiting loop)
      if(best_score > global_best_s){
        verbose("Feature that produced best score in combination ",
                "with current meta-feature set: ", best_feature, "\n")
        verbose("Score: ", best_score, "\n")
      }else{
        verbose("No further improvement in score has been found.\n")
      }

      # Increment the next loop
      i <- i+1

    } ######### End of while loop

    verbose("\n\nFinished!\n")
    verbose("Number of iterations covered: ", i, "\n")
    verbose("Best score attained over ", i ,
            " iterations: ", global_best_s, "\n")

    if(length(global_best_s_features) == 1)
      verbose("No meta-feature that improves the enrichment was found.\n")

    verbose("Features returned in FS: ", global_best_s_features, "\n")

    # We don't just want the combination (meta-feature) at the end.
    # We want all the features that make up the meta-feature
    # This can be obtained using the list of indices that were progressively
    # excluded (if at all) in the step-wise procedure
    # If returning only those features that led to the best global score
    FS_best <- FS[global_best_s_features, , drop=FALSE]

    # Make a list containing two elements.
    #The first will be the FS with the features that gave the best meta-feature
    #The second will be the score corresponding to the meta-feature
    # (named by the starting feature that led to that score)
    # Assign the name of the best meta-feature score to be the
    # starting feature that gave that score
    names(global_best_s) <- start_feature

    # Get indices of best features
    global_best_s_indices <- which(names(rowscore) %in% global_best_s_features)
    names(global_best_s_indices) <- global_best_s_features

    # Get marginal scores of best features
    marginal_best_scores <- rowscore[which(names(rowscore) %in% global_best_s_features)]
    names(marginal_best_scores) <- global_best_s_features

    return(
      list("feature_set" = FS_best,
           "input_score" = input_score,
           "score" = global_best_s,
           "best_features" = global_best_s_features,
           "best_indices" = global_best_s_indices,
           "marginal_best_scores" = marginal_best_scores,
           "cumulative_best_scores" = global_best_s_scores
      )
    )

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

  # By Default, the function returns the top N candidate search results as
  # a list of lists
  return(topn_l)

}


# Performance backward selection
#' @param FS a matrix of binary features or a SummarizedExperiment class object
#' from SummarizedExperiment package where rows represent features of interest
#' (e.g. genes, transcripts, exons, etc...) and columns represent the samples.
#' The assay of FS contains binary (1/0) values indicating the presence/absence
#' of ‘omics’ features.
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
#' arguments to be passed to the \code{custom_function}. Default is \code{NULL}.
#' @param method_alternative a character string specifies an alternative hypothesis
#' testing (\code{"two.sided"} or \code{"greater"} or \code{"less"}).
#' Default is \code{less} for left-skewed significance testing.
#' @param weights if method is \code{ks_pval} or \code{ks_score}, specifying
#' a vector of weights will perform a weighted-KS testing. Default is
#' \code{NULL}.
#' @param glob_f a vector containing the features (or row names) whose
#' union gives the best score (so far) in the search.
#' Feature names should match those of the provided FS object
#' @param glob_f_s a vector of scores corresponding to the union of the
#' specified vector of features
#' @param verbose a logical value indicates whether or not to print the
#' diagnostic messages. Default is \code{FALSE}.
#' @param ... additional parameters to be passed to the \code{custom_function}
#'
#' @noRd
#'
#' @return return the set of features with its current best score that is
#' better than the existing meta-feature best score
#'
#' @import SummarizedExperiment
forward_backward_check <- function
(
  FS,
  input_score,
  method,
  method_alternative,
  cmethod,
  custom_function,
  custom_parameters,
  weights,
  glob_f,
  glob_s,
  glob_f_scores,
  verbose = FALSE,
  ...
){

  # Set up verbose option
  options(verbose = verbose)

  verbose("Performing backward search...\n")
  verbose("Iterating over ", length(glob_f), " chosen features...\n")

  # Retrieve binary feature matrix of only global best features so far
  if(is(FS, "SummarizedExperiment")){
    gmat <- SummarizedExperiment::assay(FS)[glob_f,]
  }else{
    gmat <- FS[glob_f,]
  }

  # Here, we make a list that should store the features and their
  # corresponding meta-feature score for each leave-one-out run
  f_names <- list()
  f_indices <- list()
  f_scores <- list()
  scores <- c()

  # We want to see if leaving anyone feature out improves the overall
  # meta-feature score
  # Here we only consider previous features in the meta-feature matrix
  # (i.e. not the last one which was just added)
  for(n in seq_len(length(glob_f)-1)){
    #n=1;
    # Store features to list
    f_names[[n]] <- glob_f[-n]
    f_scores[[n]] <- glob_f_scores[-n]

    # Take the leave-one-out union of features will result in a single
    # vector to compute the scores on
    f_union <- ifelse(colSums(gmat[-n,]) == 0, 0, 1)

    # Here we are getting the union of the meta-features
    u_mat <- matrix(f_union, nrow=1, byrow=TRUE,
                    dimnames=list(c("UNION"), colnames(FS)))

    # If the class of FS is originally a SummarizedExperiment object
    # Convert the matrix to SummarizedExperiment
    # This makes sure the original class of FS object does not change
    if(is(FS, "SummarizedExperiment")){
      u_FS <- SummarizedExperiment::SummarizedExperiment(assays=u_mat)
    }else{
      u_FS <- u_mat
    }

    # Compute row scores for this meta feature
    u_rowscore <- calc_rowscore(
      FS = u_FS,
      input_score = input_score,
      meta_feature = NULL,
      method = method,
      method_alternative = method_alternative,
      cmethod = cmethod,
      custom_function = custom_function,
      custom_parameters = custom_parameters,
      weights = weights,
      verbose = FALSE,
      ...
    )

    # Get the highest best score
    scores <- c(scores, u_rowscore[which.max(u_rowscore)])

  } # end for loop

  # Set up verbose option
  options(verbose = verbose)

  # Obtain index of union meta-feature that gives the best score
  f_best_index <- which.max(scores)

  # Obtain best score of the union meta-feature
  f_best_score <- scores[f_best_index]

  # Here, we check if the new meta-feature has a computed best score better than
  # the previous meta-feature's score
  if(f_best_score > glob_s){

    verbose("Found improvement on removing existing features...\n")
    verbose("New feature set: ", f_names[[f_best_index]], "\n")
    verbose("New global best score: ", f_best_score, "\n")

    # Return a set of features that gave a better score than
    # the existing best score and its new best score as well
    return(list(best_features = f_names[[f_best_index]],
                best_scores = f_scores[[f_best_index]],
                score = f_best_score))

  }else{

    # Don't change anything. Return the existing
    # best set of features and its corresponding score
    return(list(best_features = glob_f,
                best_scores = glob_f_scores,
                score = glob_s))

  }

}
