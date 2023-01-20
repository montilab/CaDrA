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
#' @param FS a SummarizedExperiment object containing binary features where
#' rows represent features of interest (e.g. genes, transcripts, exons, etc...)
#' and columns represent the samples.
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
#' @param metric a character string specifies a metric to search
#' for best features. \code{"pval"} or \code{"stat"} may be used which
#' corresponding to p-value or score statistic. Default is \code{pval}.
#' NOTE: \code{Revealer} method only utilized score statistics values
#' (no p-values).
#' @param weight if method is \code{ks}, specifies a vector of weight
#' will perform a weighted-KS testing. Default is \code{NULL}.
#' @param search_start a list of character strings (separated by commas)
#' which specifies feature names within the expression set object to start
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
#' @param best_score_only a logical value indicates whether or not the
#' function should return only the score corresponding to the search results.
#' Default is \code{FALSE}.
#' @param do_plot a logical value indicates whether or not to plot the
#' resulting meta-feature matrix. NOTE: plot can only be produced if resulting
#' meta-feature matrix contains more than 1 feature
#' (e.g. length(search_start) > 1  or top_N > 1). Default is \code{FALSE}.
#' @param warning a logical value indicates whether or not to print the
#' diagnostic messages. Default is \code{FALSE}.
#'
#' @return If \code{best_score_only} is set to \code{TRUE}, the function
#' returns a list object with the score corresponding to the union of the
#' search meta-feature. If \code{best_score_only} is set to \code{FALSE},
#' a list containing the FS object (pertaining to the returned meta-feature)
#' as well as its corresponding score and observed input scores are returned.
#' @examples
#'
#' # Load pre-computed expression set
#' data(sim_FS)
#'
#' # Load pre-computed input-score
#' data(sim_Scores)
#'
#' # Define additional parameters and run the function
#' candidate_search_result <- candidate_search(
#'   FS = sim_FS, input_score = sim_Scores, method = "ks",
#'   alternative = "less", weight = NULL, metric = "pval",
#'   search_start = NULL, top_N = 3, search_method = "both",
#'   max_size = 7, best_score_only = FALSE
#' )
#'
#' @export
#' @import SummarizedExperiment methods
candidate_search <- function(
    FS,
    input_score,
    method = c("ks", "wilcox", "revealer", "custom"),
    custom_function = NULL,
    custom_parameters = NULL,
    alternative = c("less", "greater", "two.sided"),
    metric = c("pval", "stat"),
    weight = NULL,
    search_start = NULL,
    top_N = 1,
    search_method = c("both", "forward"),
    max_size = 7,
    best_score_only = FALSE,
    do_plot = FALSE,
    warning = FALSE
){

  # Set up verbose option
  options(verbose = warning)

  method <- match.arg(method)
  metric <- match.arg(metric)
  alternative <- match.arg(alternative)
  search_method <- match.arg(search_method)

  # Check of FS and input_score are valid inputs
  check_data_input(FS = FS, input_score = input_score, warning=TRUE)
  
  # Select the appropriate method to compute scores based on
  # skewness of a given binary matrix
  s <- calc_rowscore(
    FS = FS,
    input_score = input_score,
    method = method,
    alternative = alternative,
    weight = weight,
    seed_names = NULL,
    custom_function = custom_function,
    custom_parameters = custom_parameters,
    warning = FALSE
  )

  # Check if the returning result has one or two columns:
  # score or p_value or both
  if("score" %in% colnames(s) & !"p_value" %in% colnames(s) & metric == "pval"){
    warning("metric = 'pval' is provided but the method function only ",
            "return score values. Thus, using 'stat' as metric to search ",
            "for best features.")
    metric <- "stat"
  }else if("p_value" %in% colnames(s) & !"score" %in% colnames(s) & metric == "stat"){
    warning("metric provided is 'stat' but the method function only ",
            "return p-values. Thus, using 'pval' as metric to search ",
            "for best features.")
    metric <- "pval"
  }

  # Score returned by the given method
  s.stat <- if("score" %in% colnames(s)){ s[,"score"] }
  s.pval <- if("p_value" %in% colnames(s)){ s[,"p_value"] }

  # compute the scores according to the provided metric
  score <- ifelse(rep(metric, nrow(FS)) %in% "pval", sign(s.stat)*s.pval, s.stat)

  verbose("Using ", metric, " as measure of improvement measure...\n")

  ###### FEATURE PRE-RANKING #####
  ##################################

  # Re-order FS in decreasing order of user-defined score (stat or pval)
  # This comes in handy when doing the top-N evaluation of
  # the top N 'best' features
  score.rank <- if(metric != "pval")
    order(score, decreasing=TRUE) else
      order(-sign(score), score)

  # Re-order the FS
  FS <- FS[score.rank,]

  # Re-order the computed scores
  score <- score[score.rank]

  ###### INITIALIZE VARIABLES ###########
  #######################################

  # Check if top_N is given and is numeric
  top_N <- as.integer(top_N)

  # Check if search_start is given
  if(is.null(search_start)){

    if(is.na(top_N) || length(top_N)==0 || top_N <= 0){
      stop("Please specify a NUMERIC top_N value to evaluate over top N ",
           "features (top_N must be >= 1).\n")
    }

    if(top_N > nrow(FS))
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

    if(length(search_start) == 0 || any(!search_start %in% rownames(FS)))
      stop("Provided starting feature does not exist among FS's rownames.\n\n")

    # Get the index of the search_start strings and start
    # the search with the defined indexes
    search_feature_index <- which(rownames(FS) %in% search_start)

  } # end else (!is.null)

  ## Check the search_method variable ####
  back_search <- ifelse(search_method == "both", TRUE, FALSE)

  ## Check the max_size variable ####
  max_size <- as.integer(max_size)

  if(is.na(max_size) || length(max_size)==0 || max_size <= 0 || max_size > nrow(FS)){
    stop("Please specify an integer value specifies a maximum size ",
         "that a meta-feature can extend to do for a given search ",
         "(max_size must be >= 1)",
         "and max_size must be lesser than the number of features in",
         "FS\n")
  }

  # Check the best_score_only variables
  if(!best_score_only %in% c(TRUE, FALSE)){
    stop("Please specify a logical value TRUE/FALSE ",
         "for best_score_only variable.\n")
  }

  # Performs the search over the feature indices
  topn_l <- lapply(seq_along(search_feature_index), function(x){
    # x=1;
    # Evaluate over top N features using their indexes
    best.s.index <- search_feature_index[x]

    # Print the featureData for this starting point feature so that we are aware
    # Here we assume the FS's fData is included as rownames
    start.feature <- rownames(FS)[best.s.index]
    best.feature <- start.feature
    best.s <- score[best.s.index]

    verbose("Top Feature ", x, ": ", start.feature, "\n")
    verbose("Score: ", best.s, "\n")

    #Fetch the vector corresponding to best score
    #Set this as the initial 'meta-feature'
    best.meta <- as.numeric(assay(FS)[best.s.index,])

    #counter variable for number of iterations
    i <- 0

    #Variable to store best score attained over all iterations
    #initialize this to the starting best score
    global.best.s <- best.s

    # Vector of features in the (growing) obtained meta-feature.
    # Begin with just the starting feature
    global.best.s.features <- c()

    ###### BEGIN ITERATIONS ###############
    #######################################

    verbose("\nBeginning candidate search...\n")

    while((ifelse(metric=="pval",
                  (sign(best.s) > 0 & (abs(best.s) < abs(global.best.s))),
                  best.s > global.best.s) | i == 0) &
          (length(global.best.s.features) < max_size)){

      verbose("\n\n")
      verbose("Iteration number ", (i+1), " ..\n")
      verbose("Global best score: ", global.best.s, "\n")
      verbose("Previous score: ", best.s, "\n")

      # Update scores and feature set since entry into the
      # loop means there is an improvement (iteration > 0)
      global.best.s <- best.s
      global.best.s.features <- c(global.best.s.features, best.feature)

      verbose("Current feature set: ", global.best.s.features, "\n")

      if(i != 0){

        verbose("Found feature that improves score!\n")

        # Update the new best meta feature (from meta mat)
        best.meta <- meta.mat[hit.best.s.index,]

        # Add that index to the group of indices to be excluded
        # for subsequent checks
        # Here we go off the rownames in the original matrix to
        # find which index to exclude from the FS in subsequent iterations
        best.s.index <- c(best.s.index, which(rownames(FS) == best.feature))

      }

      # Perform a backward check on the list of existing features and
      # update global scores/feature lists accordingly
      if(length(global.best.s.features) > 3 & back_search == TRUE){

        backward_search_results <- forward_backward_check(
          FS = FS,
          input_score = input_score,
          glob.f = global.best.s.features, # Global feature set so far
          glob.f.s = global.best.s,
          method = method,
          custom_function = custom_function,
          custom_parameters = custom_parameters,
          alternative = alternative,
          metric = metric,
          weight = weight
        )

        # Update globlal features, scores
        global.best.s.features <- backward_search_results[["best.features"]]
        global.best.s <- backward_search_results[["best.scores"]]

        # Update best.meta based on feature set
        best.meta <- as.numeric(ifelse(
          colSums(assay(FS)[global.best.s.features,]) == 0, 0, 1))

      }

      # Take the OR function between that feature and all
      # other features to see which gives the best score
      # Keep in mind, the number of rows in meta.mat keeps reducing by one
      # each time we find a hit that improves the score
      verbose("Forming meta-feature matrix with all other features in dataset.")

      # Here "*1" is used to convert the boolean back to integer 1's and 0's
      # Notice we remove anything in best.s.index from the original matrix
      # first to form the meta matrix.
      meta.mat <- base::sweep(assay(FS)[-best.s.index,], 2, best.meta, `|`)*1

      # Check if there are any features that are all 1's generated on
      # taking the union
      # We cannot compute statistics for such features and they thus need
      # to be filtered out
      if(any(rowSums(meta.mat) == ncol(meta.mat))){
        warning("Features with all 1's generated upon taking matrix union...",
                "Removing such features before progressing...\n")
        meta.mat <- meta.mat[rowSums(meta.mat) != ncol(meta.mat),]
      }

      # With the newly formed 'meta-feature' matrix, compute directional
      # scores and choose the feature that gives the best score
      s <- calc_rowscore(
        FS = meta.mat,
        input_score = input_score,
        method = method,
        alternative = alternative,
        weight = weight,
        seed_names = NULL,
        custom_function = custom_function,
        custom_parameters = custom_parameters,
        warning = FALSE
      )

      # Score returned by the given method
      s.stat <- if("score" %in% colnames(s)){ s[,"score"] }
      s.pval <- if("p_value" %in% colnames(s)){ s[,"p_value"] }

      # Take signed pval or stat depending on user-defined metric
      # This will be the same length as nrow(meta.mat)
      scores <- ifelse(rep(metric, nrow(meta.mat)) %in% "pval",
                       sign(s.stat)*s.pval, s.stat)

      # Find index of feature that gives lowest scores when
      # combined with chosen starting feature
      if(metric != "pval"){
        hit.best.s.index <- which.max(scores)
      } else { #If signed pvalues
        hit.best.s.index <- order(-sign(scores), scores)[1]
        #Top p-value ordered by sign and numerical value;
        #This is the index within the meta matrix
      }

      # The best score from the meta matrix
      best.s <- scores[hit.best.s.index]

      # Find which feature produced that score,
      # in combination with meta feature used
      # We go from index to rowname space here in the meta matrix
      # We can do this because rownames are preserved between the
      # original and meta features on using sweep()
      best.feature <- rownames(meta.mat)[hit.best.s.index]

      verbose("Feature that produced best score in combination ",
              "with previous meta-feature: ", best.feature )
      verbose("Score: ", best.s)

      # If no improvement (exiting loop)
      if(ifelse(metric == "pval", sign(best.s) < 0 |
                (abs(best.s) >= abs(global.best.s)), best.s <= global.best.s)){
        verbose("No further improvement in score has been found.")
      }

      #Increment counter
      i <- i+1

    } #########End of while loop

    verbose("\n\nFinished!\n\n")
    verbose("Number of iterations covered: ", i, "\n")
    verbose("Best score attained over iterations: ", global.best.s, "\n")

    if(length(global.best.s.features) == 1){
      verbose("No meta-feature that improves the enrichment was found...\n")
    }

    verbose("Features returned in FS: ", global.best.s.features, "\n")
    verbose("\n\n")

    # We don't just want the combination (meta-feature) at the end.
    # We want all the features that make up the meta-feature
    # This can be obtained using the list of indices that were progressively
    # excluded (if at all) in the step-wise procedure
    # If returning only those features that led to the best global score
    FS.best <- FS[global.best.s.features,]

    # Here, give the returned FS an annotation based on the starting
    # feature that gave these best results
    assayNames(FS.best) <- start.feature
    colnames(FS.best) <- colnames(FS)

    # Make a list contaning two elements.
    #The first will be the FS with the features that gave the best meta-feature
    #The second will be the score corresponding to the meta-feature
    # (named by the starting feature that led to that score)
    # Assign the name of the best meta-feature score to be the
    # starting feature that gave that score
    names(global.best.s) <- start.feature

    return(list("method" = method,
                "metric" = metric,
                "feature_set" = FS.best,
                "input_score" = input_score,
                "score" = global.best.s))

  })

  # length of search_feature must be greater 2 in order to generate topn_plot
  if(do_plot) topn_plot(topn_list = topn_l)

  # best_score_only
  if(best_score_only == TRUE){

    scores_l <- lapply(seq_along(topn_l), function(l){ topn_l[[l]][['score']] })

    # Working with scores for each top N run
    s <- unlist(scores_l)

    # Fetch the best score from the iterations
    # NEEDS UPDATING TO ACCOMODATE STATISTIC
    if(metric == "pval"){
      best_score <- s[order(s, decreasing = FALSE)][1]
      #Based on the p-values, the lowest value will be the most significant
    }else{
      best_score <- s[order(s, decreasing = TRUE)][1]
    }

    return(best_score)

  }

  # By Default, the function returns the top N candidate search results as a list of lists
  return(topn_l)

}


# Performance backward selection
# FS: # an Expression Set object with the same sample
# ordering and features as processed by the stepwise.search() function

# glob.f: a vector containing the feature (row) names whose
# union gives the best score (so far) in the search.
# Feature names should match those of the provided expression set object

# glob.f.s: score corresponding to the union of the specified vector of features

# metric: # a character string specifying which metric to use for
# stepwise search criteria. One of either 'pval' or 'stat' may be used,
# corresponding to the  p-value or statistic. Uses value passed in the
# candidate_search() function
forward_backward_check <- function
(
  FS,
  input_score,
  glob.f,
  glob.f.s,
  method,
  custom_function,
  custom_parameters,
  alternative,
  metric,
  weight
){

  verbose("Performing backward search...\n")
  verbose("Iterating over ", length(glob.f), " chosen features...\n")

  # Matrix of only global best features so far
  gmat <- assay(FS[glob.f,])
  rownames(gmat) <- glob.f

  # Here, we make a list that should store the features and their
  # corresponding meta-feature score for each leave-one-out run
  f.names <- list()
  f.scores <- c()

  # We want to see if leaving anyone feature out improves the overall
  # meta-feature score
  # Here we only consider previous features in the meta-feature to remove
  # (i.e. not the last one which was just added)
  for(n in seq_len(length(glob.f)-1)){
    #n=1;
    f.names[[n]] <- glob.f[-n]

    # Take leave-one-out union of features from matrix
    # This will result in a single vector to compute the scores on
    u <- ifelse(colSums(gmat[-n,]) == 0, 0, 1)

    # Compute scores for this meta feature
    # Here we suprress warnings just to avoid messages warning-related single
    # vector score computation (nrow(mat) < 2)
    u.mat <- matrix(t(matrix(u)),
                    nrow=1,
                    byrow=TRUE,
                    dimnames=list(c("sum"),
                                  colnames(FS)))
    
    s <- calc_rowscore(
      FS = u.mat,
      input_score = input_score,
      method = method,
      alternative = alternative,
      weight = weight,
      seed_names = NULL,
      custom_function = custom_function,
      custom_parameters = custom_parameters,
      warning = FALSE
    )

    # Check if the returning result has one or two columns:
    # score or p_value or both
    if(ncol(s) == 1){
      if(colnames(s) == "score" & metric == "pval"){
        warning("metric = 'pval' is provided but the method ",
                "function only return score values. ",
                "Thus, using 'stat' as metric to search for best features.")
        metric <- "stat"
      }else if(colnames(s) == "p_value" & metric == "stat"){
        warning("metric provided is 'stat' but the method function only ",
                "return p-values. Thus, using 'pval' as metric to ",
                "search for best features.")
        metric <- "pval"
      }
    }

    # Score returned by the given method
    s.stat <- if("score" %in% colnames(s)){ s[,"score"] }
    s.pval <- if("p_value" %in% colnames(s)){ s[,"p_value"] }

    score <- ifelse(metric %in% "pval", sign(s.stat)*s.pval, s.stat)

    f.scores <- c(f.scores, score)

  } # end for loop

  if(metric != "pval"){
    f.best.index <- which.max(f.scores)
  } else { #If signed pvalues
    f.best.index <- order(-sign(f.scores), f.scores)[1]
  }

  f.best.score <- f.scores[f.best.index]

  # Check if any one of the computed scores has a better score than the
  # entire meta-feature's score
  if(ifelse(metric == "pval",
            (sign(f.best.score) > 0 & abs(f.best.score) < abs(glob.f.s)),
            f.best.score > glob.f.s)){

    verbose("Found improvement on removing existing feature...\n")
    verbose("New feature set: ", f.names[[f.best.index]], "\n")
    verbose("New global best score: ", f.best.score, "\n")

    # Return the set of features that gave a better score than
    # the existing best score, and the score as well
    return(list(best.features=f.names[[f.best.index]], best.scores=f.best.score))

  } else{

    # Don't change anything. Return the existing
    # best set of features and the corresponding score
    return(list(best.features=glob.f, best.scores=glob.f.s))

  }

}


