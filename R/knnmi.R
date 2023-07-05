#'
#' Mutual information estimation MI(X;Y) when the target (X) and features (Y) are continuous.
#'
#' This implements the algorithm described in: 
#' https://doi.org/10.1103/PhysRevE.69.066138
#' Alexander Kraskov, Harald Stogbauer, and Peter Grassberger
#' Phys. Rev. E 69, 066138 ?? Published 23 June 2004; Erratum Phys. Rev. E 83, 019903 (2011)
#'
#' Compute mutual information of \code{target} and \code{features}
#' where \code{target} and \code{features} are both continuous
#' @param target input vector of length N.
#' @param features input vector of length N or a matrix of size MxN.
#' @param k number of nearest neighbors.
#' @useDynLib CaDrA _mutual_inf_cc
#'
#' @return a double-precision value - mutual information estimation for
#' vectors \code{target} and \code{features}.
#' @examples
#'
#' data(mutual_info_df)
#' set.seed(654321)
#' mutual_inf_cc(mutual_info_df$Xc, t(mutual_info_df$Zc_XcYc))
#' ## 0
#'
#' mutual_inf_cc(mutual_info_df$Yc, t(mutual_info_df$Zc_XcYc))
#' ## 0.2738658
#'
#' @export
mutual_inf_cc <- function(target, features, k=3L){
  if (is.vector(features)) {
    stopifnot( "target and features vectors must have the same length"=length(target) == length(features) )
  } else {
    stopifnot( "number of columns in features matrix should be equal to length of the target vector"=
                 length(target) == ncol(features) )
  }
  stopifnot("k must be less than the length of target"=k < length(target))
  
  # Make sure the inputs are double precision floats
  if (!is.double(target)) {
    storage.mode(target) <- "double"
  }
  if (!is.double(features)) {
    storage.mode(features) <- "double"
  }
  
  res <- .Call('_mutual_inf_cc', target, features, as.integer(k))
  res
}



#' Mutual information estimation MI(X;Y) when the target (X) is continuous and the
#' features (Y) are discrete.
#'
#' This implements the algorithm described in:
#' https://doi.org/10.1371/journal.pone.0087357
#' Ross BC (2014) Mutual Information between Discrete and Continuous Data Sets. 
#' PLoS ONE 9(2): e87357. 
#' 
#' Compute mutual information of \code{target} and \code{y}
#' where the \code{target} is continuous and \code{features} are discrete.
#' @param target input vector of length N.
#' @param features input vector of length N or a matrix of size MxN.
#' @param k number of nearest neighbors.
#' @useDynLib CaDrA _mutual_inf_cd
#'
#' @return a double-precision vector - mutual information estimation for
#' vectors \code{target} and \code{features}.
#' @examples
#'
#' data(mutual_info_df)
#' set.seed(654321)
#' mutual_inf_cd(mutual_info_df$Zc_XdYd, t(mutual_info_df$Xd))
#' ## 0.128029
#'
#' M <- cbind(mutual_info_df$Xd, mutual_info_df$Yd)
#' mutual_inf_cd(mutual_info_df$Zc_XdYdWd, t(M))
#' ## 0.1070804 0.1041177
#' 
#' @export
mutual_inf_cd <- function(target, features, k=3L){
  if (is.vector(features)) {
    stopifnot( "target and features vectors must have the same length"=length(target) == length(features) )
  } else {
    stopifnot( "number of columns in features matrix should be equal to length of target"=
                 length(target) == ncol(features) )
  }
  stopifnot("k must be less than the length of target"=k < length(target))
  
  # Make sure the target is double precision floats
  if (!is.double(target)) {
    storage.mode(target) <- "double"
  }
  # Make sure this is an integer array - a numeric array containing
  # integer values isn't quite sufficient, they need to be actual
  # integers so the C/C++ code works properly.
  if (!is.integer(features)) {
    storage.mode(features) <- "integer"
  }
  res <- .Call('_mutual_inf_cd', target, features, as.integer(k))
  res
}


#'
#' Conditional mutual information estimation CMI(X;Y|Z) where X is a continuous vector.
#' The input Y and conditional input Z can be vectors or matrices. If Y and Z
#' are discrete then they must be numeric or integer valued.
#'
#' This implements the CMI algorithm described in: 
#' 
#' https:doi.org/10.1016/j.eswa.2012.05.014
#' 
#' Alkiviadis Tsimpiris, Ioannis Vlachos, Dimitris Kugiumtzis,
#' Nearest neighbor estimate of conditional mutual information in feature selection,
#' Expert Systems with Applications,
#' Volume 39, Issue 16, 2012, Pages 12697-12708
#' 
#' @param x vector of size N.
#' @param Y input vector of length N or a matrix of size NxM.
#' @param Z conditional input vector of length N or a matrix of size NxM.
#' @param k number of nearest neighbors.
#' @useDynLib CaDrA _cond_mutual_inf
#'
#' @return a double-precision vector
#'
#' @examples
#' data(mutual_info_df)
#' set.seed(654321)
#' cond_mutual_inf(mutual_info_df$Zc_XcYc,
#'                        mutual_info_df$Xc, t(mutual_info_df$Yc))
#' ## 0.2936858
#' 
#' M <- cbind(mutual_info_df$Xc, mutual_info_df$Yc)
#' ZM <- cbind(mutual_info_df$Yc, mutual_info_df$Wc)
#' cond_mutual_inf(mutual_info_df$Zc_XcYcWc, t(M), t(ZM))
#' ## 0.1171533 0.2192397
#'
#' @export
cond_mutual_inf <- function(x, Y, Z, k=3L){
  # conditional mutual information:  CMI(X, Y|Z)
  # x: vector of length N
  # Y: vector of length N or matrix of size MxN
  # Z: vector of length N or matrix of size MxN
  
  # Check that sizes match before continuing. 
  stopifnot("k must be less than the length of the input vector x"=k < length(x))
  
  # The "case" value is used in the C code to let it pick the right code
  # path without having to re-check for vector vs matrix in C.
  case <- 0
  if (is.vector(Y) && is.vector(Z)) {
    # When they're both vectors make sure their length is the same as x.  
    stopifnot( "x and Y must have the same length"=length(x) == length(Y) )
    stopifnot( "x and Z must have the same length"=length(x) == length(Z) )
  } else if (is.vector(Y) && is.matrix(Z)) {
    case <- 1
    # mixed vector & matrix
    stopifnot( "x and Y must have the same length"=length(x) == length(Y) )
    stopifnot( "Number of Z columns must be the same as the length of x"=length(x) == ncol(Z) )
  } else if (is.vector(Z) && is.matrix(Y)) {
    # mixed vector & matrix
    case <- 2
    stopifnot( "x and Z must have the same length"=length(x) == length(Z) )
    stopifnot( "Number of Y columns must be the same as the length of x"=length(x) == ncol(Y) )    
  }  else if (is.matrix(Z) && is.matrix(Y)) {
    # Both Y and Z are matrices.
    case <- 3
    stopifnot( "x and Y must have the same length"=length(x) == ncol(Y) )
    stopifnot( "x and Z must have the same length"=length(x) == ncol(Z) )
    stopifnot( "Y and Z must be the same size"=dim(Y) == dim(Z) )
  } else {
    # Incorrect arguments...
    stop("Y and Z must be vectors or matrices of the correct dimensionality.")
  }
  # If any of the input arguments are not double precision, convert them.
  if (!is.double(x)) {
    storage.mode(x) <- "double"
  }
  if (!is.double(Y)) {
    storage.mode(Y) <- "double"
  }
  if (!is.double(Z)) {
    storage.mode(Z) <- "double"
  }
  
  res <- .Call('_cond_mutual_inf', x, Y, Z, as.integer(k), as.integer(case))
  res
}


