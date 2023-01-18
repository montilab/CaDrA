#'
#' Wilcoxon Rank Sum Scoring Method
#'
#' Compute directional Wilcoxon rank sum score for each row of a
#' given binary feature matrix
#'
#' @param FS a feature set of binary features. It can be an expression matrix or
#' a \code{SummarizedExperiment} class object from SummarizedExperiment package
#' @param input_score a vector of continuous scores representing a phenotypic
#' readout of interest such as protein expression, pathway activity, etc.
#' The \code{input_score} object must have names or labels that match the column
#' names of FS object.
#' @param alternative a character string specifies an alternative
#' hypothesis testing (\code{"two.sided"} or \code{"greater"} or
#' \code{"less"}). Default is \code{less} for left-skewed significance testing.
#' @param warning a logical value indicates whether or not to print the 
#' diagnostic messages. Default is \code{TRUE}
#' 
#' @noRd
#'
#' @return A matrix with two columns: \code{score} and \code{p_value}
#' @import SummarizedExperiment
wilcox_rowscore <- function
(
  FS,
  input_score,
  alternative = c("less", "greater", "two.sided"),
  warning = TRUE
)
{

  alternative <- match.arg(alternative)
  
  # Check of FS and input_score are valid inputs
  if(warning == TRUE) check_data_input(FS = FS, input_score = input_score, warning=warning)
  
  # Get the feature names
  feature_names <- rownames(FS)
  
  # Wilcox is a ranked-based method
  # So we need to sort input_score from highest to lowest values
  input_score <- sort(input_score, decreasing=TRUE)
  
  # Re-order the matrix based on the order of input_score
  FS <- FS[, names(input_score)]
  
  # Extract the feature binary matrix
  if(class(FS)[1] == "SummarizedExperiment"){
    mat <- as.matrix(SummarizedExperiment::assay(FS))
  }else if(class(FS)[1] == "matrix"){
    mat <- as.matrix(FS)
  }else{
    mat <- matrix(t(FS), nrow=1, byrow=TRUE,
                  dimnames=list(feature_names, names(FS)))
  }

  # Since input_score is already ordered from largest to smallest
  # We can assign ranks as 1:N (N: number of samples)
  ranks <- seq(1, ncol(mat))

  #Compute the wilcox rank sum statitic and p-value per row in the matrix
  wilcox <- apply(X=mat, MARGIN=1, function(x){
    wilcox_score(
      x = ranks[which(x==1)],
      y = ranks[which(x==0)],
      alternative = alternative
    )
  })

  # Convert results as matrix
  # Wilcox method returns both score and p-values
  wilcox_mat <- matrix(NA, nrow=nrow(mat), 
                       ncol=2, byrow=TRUE, 
                       dimnames=list(rownames(mat), c("score", "p_value")))
  
  wilcox_mat[,1] <- wilcox[1,]
  wilcox_mat[,2] <- wilcox[2,]

  return(wilcox_mat)

}



#' Compute rank sum scores for a given binary feature
#'
#' @param x an integer ranked values for group 1
#' @param y an integer ranked values for group 2
#' @param mu a number uses as an optional parameter to form a null hypothesis.
#' Default is \code{0}.
#' @param alternative alternative hypothesis for p-value calculation
#' (\code{"two.sided"} or \code{"greater"} or \code{"less"}).
#' Default is \code{less} for left-skewed significance testing.
#' @param paired whether to perform paired test. Default is \code{FALSE}.
#' @param exact whether to compute exact p-value. Default is \code{FALSE}.
#' @param correct whether to consider continuity correction for p-value.
#' Default is \code{TRUE}.
#'
#' @noRd
#'
#' @return a vector with two values: \code{score} and \code{p_value}
#'
#' @importFrom stats pnorm pwilcox
wilcox_score <- function
(
  x,
  y,
  mu = 0,
  alternative = c("less", "greater", "two.sided"),
  paired = FALSE,
  exact = FALSE,
  correct = TRUE
)
{

  alternative <- match.arg(alternative)
  
  if (!missing(mu) && ((length(mu) > 1L) || !is.finite(mu)))
    stop("'mu' must be a single number")
  if (!is.numeric(x))
    stop("'x' must be numeric")
  if (!is.numeric(y))
    stop("'y' must be numeric")

  DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))

  x <- x[is.finite(x)]
  y <- y[is.finite(y)]

  if (length(x) < 1L) {
    #print(x)
    stop("not enough (finite) 'x' observations")}

  if (length(y) < 1L){
    #print(y)
    stop("not enough 'y' observations")}

  METHOD <- "Wilcoxon rank sum test"

  ##### Modification ######
  # Take input as ranks instead of continuous measures
  # (normally internally ranked: see below)
  r <- c(x,y)

  #r <- rank(c(x - mu, y))

  n.x <- as.double(length(x))
  n.y <- as.double(length(y))

  if (is.null(exact))
    exact <- (n.x < 50) && (n.y < 50)

  STATISTIC <- c(W = sum(r[seq_along(x)]) - n.x * (n.x + 1)/2)
  TIES <- (length(r) != length(unique(r)))

  if (exact && !TIES) {

    PVAL <- switch(alternative, two.sided = {
      p <- if (STATISTIC > (n.x * n.y/2))
        pwilcox(STATISTIC - 1, n.x, n.y,
                lower.tail = FALSE) else pwilcox(STATISTIC, n.x, n.y)
      min(2 * p, 1)
    }, greater = {
      pwilcox(STATISTIC - 1, n.x, n.y, lower.tail = FALSE)
    }, less = pwilcox(STATISTIC, n.x, n.y))

  } else {

    NTIES <- table(r)
    z <- STATISTIC - n.x * n.y/2
    SIGMA <- sqrt((n.x * n.y/12) * ((n.x + n.y + 1) -
                                      sum(NTIES^3 - NTIES)/
                                      ((n.x + n.y) * (n.x + n.y - 1))))

    if (correct) {
      CORRECTION <- switch(alternative, two.sided = sign(z) *
                             0.5, greater = 0.5, less = -0.5)
      METHOD <- paste(METHOD, "with continuity correction")
    }

    z <- (z - CORRECTION)/SIGMA

    PVAL <- switch(alternative, less = pnorm(z),
                   greater = pnorm(z, lower.tail = FALSE),
                   two.sided = 2 * min(pnorm(z),
                                       pnorm(z, lower.tail = FALSE)))

    if (exact && TIES)
      warning("cannot compute exact p-value with ties")

  }

  names(mu) <- ifelse (paired || !is.null(y), "location shift", "location")

  RVAL <- list(statistic = STATISTIC,
               parameter = NULL, p.value = as.numeric(PVAL),
               null.value = mu, alternative = alternative, method = METHOD,
               data.name = DNAME)

  class(RVAL) <- "htest"

  return(c(score=RVAL$statistic, p_value=RVAL$p.value))

}

