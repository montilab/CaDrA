#'
#' Correlation Scoring Method
#'
#' Calculate correlation scores between the input score and a feature set
#'
#' @param FS a matrix of binary features where
#' rows represent features of interest (e.g. genes, transcripts, exons, etc...)
#' and columns represent the samples.
#' @param input_score a vector of continuous scores representing a phenotypic
#' readout of interest such as protein expression, pathway activity, etc.
#' The \code{input_score} object must have names or labels that match the column
#' names of FS object.
#' @param meta_feature a vector of one or more features representing known causes
#' of activation or features associated with a response of interest
#' (\code{e.g. input_score}). Default is NULL.
#' @param corr_method correlation method - Pearson or Spearman. Default is
#' \code{spearman}
#' @return  returns a vector of row-wise scores where its labels or names
#' must match the row names of \code{FS} object
#'
#' @noRd
#'
#' @examples
#' mat <- matrix(c(1,0,1,0,0,0,0,0,1,0,
#'                 0,0,1,0,1,0,1,0,0,0,
#'                 0,0,0,0,1,0,1,0,1,0), nrow=3)
#'
#' colnames(mat) <- 1:10
#' row.names(mat) <- c("TP_1", "TP_2", "TP_3")
#'
#' set.seed(42)
#' input_score = rnorm(n = ncol(mat))
#' names(input_score) <- colnames(mat)
#'
#' corr_rs <- corr_rowscore(
#'    FS = mat,
#'    input_score = input_score,
#'    meta_feature = NULL,
#'    corr_method = "spearman"
#' )
#'
corr_rowscore <- function
(
  FS,
  input_score,
  meta_feature = NULL,
  corr_method = c("spearman", "pearson")
)
{

  corr_method <- match.arg(corr_method)

  # Check if meta_feature is provided
  if(!is.null(meta_feature)){
    # Getting the position of the known meta features
    locs <- match(meta_feature, row.names(FS))

    # Taking the union across the known meta features
    if(length(locs) > 1) {
      meta_vector <- as.numeric(ifelse(colSums(FS[locs,]) == 0, 0, 1))
    }else{
      meta_vector <- as.numeric(FS[locs, , drop=FALSE])
    }

    # Remove the meta features from the binary feature matrix
    # and taking logical OR btw the remaining features with the meta vector
    FS <- base::sweep(FS[-locs, , drop=FALSE], 2, meta_vector, `|`)*1

    # Check if there are any features that are all 1s generated from
    # taking the union between the matrix
    # We cannot compute statistics for such features and thus they need
    # to be filtered out
    if(any(rowSums(FS) == ncol(FS))){
      verbose("Features with all 1s generated from taking the matrix union ",
              "will be removed before progressing...\n")
      FS <- FS[rowSums(FS) != ncol(FS), , drop=FALSE]
      # If no features remained after filtering, exist the function
      if(nrow(FS) == 0) return(NULL)
    }
  }

  # correlation
  scores <- apply(FS, 1, function (x) cor(x, input_score, method = corr_method))

  return(scores)

}
