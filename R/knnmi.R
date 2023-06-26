#'
#' Mutual information estimation for continuous case of 2 vectors.
#'
#' Compute mutual information of \code{x} and \code{y}
#' where \code{x} and \code{y} are both continuous
#' @param x input vector.
#' @param y input vector of the same length as x.
#' @param k number of nearest neighbors.
#' @param seed integer random seed. If seed is <=0 a random seed is used.
#' @useDynLib CaDrA _mutual_inf_cc_vec
#'
#' @return a double-precision value - mutual information estimation for
#' vectors \code{x} and \code{y}.
#' @examples
#'
#' data(mutual_info_df)
#' mutual_inf_cc_vec(mutual_info_df$Xc, mutual_info_df$Zc_XcYc)
#' ## 0
#'
#' mutual_inf_cc_vec(mutual_info_df$Yc, mutual_info_df$Zc_XcYc)
#' ## 0.2738658
#'
#' @export
mutual_inf_cc_vec <- function(x, y, k=3L, seed=0L){
  stopifnot( "x and y must have the same length"=length(x) == length(y) )
  stopifnot("k must be less than the length of x"=k < length(x))

  res <- .Call('_mutual_inf_cc_vec', x, y, as.integer(k), as.integer(seed))
  res
}



#'
#' Mutual information estimation for continuous case of an NxM matrix.
#'
#' Compute conditional mutual information of \code{x} and matrix \code{M}
#' where \code{x} and \code{M} are both continuous
#' @param x input vector.
#' @param M input matrix. The number of rows should be equal to the
#' length of vector x
#' @param k number of nearest neighbors.
#' @param seed integer random seed. If seed is <=0 a random seed is used.
#' @useDynLib CaDrA _mutual_inf_cc_mat
#'
#' @return a vector of length equal to the number of columns of matrix M.
#' @examples
#'
#' data(mutual_info_df)
#'
#' M <- cbind(mutual_info_df$Xc, mutual_info_df$Yc)
#' mutual_inf_cc_mat(mutual_info_df$Zc_XcYcWc, M)
#' ## 0.000000 0.199844
#'
#' @export
mutual_inf_cc_mat <- function(x, M, k=3L, seed=0L){

  stopifnot( "number of rows in matrix M should be equal to length of x"=
               length(x) == nrow(M) )
  stopifnot("k must be less than the length of x"=k < length(x))

  res <- .Call('_mutual_inf_cc_mat', x, M, as.integer(k), as.integer(seed))
  res

}



#'
#' Mutual information estimation for continuous/discrete case.
#'
#' Compute mutual information of \code{x} and \code{y}
#' where \code{x} is continuous vector and \code{y} are discrete
#' @param x input (continuous) vector.
#' @param y input (discrete) vector. It should have the same length as x.
#' @param k number of nearest neighbors.
#' @param seed integer random seed. If seed<=0 a random seed is used.
#' @param use_cc (logical) if TRUE the algorithm falls into continuous/continuous case
#' @useDynLib CaDrA _mutual_inf_cd_vec
#'
#' @return a double-precision value - mutual information estimation for
#' vectors \code{x} and \code{y}.
#' @examples
#'
#' data(mutual_info_df)
#'
#' mutual_inf_cd_vec(mutual_info_df$Zc_XdYd, mutual_info_df$Xd)
#' ## 0.128029
#'
#' @export
mutual_inf_cd_vec <- function(x, y, k=3L, seed=0L, use_cc=FALSE){
  stopifnot( "x and y must have the same length"=length(x) == length(y) )
  stopifnot("k must be less than the length of x"=k < length(x))

  if (!is.integer(y)) {
    y <- as.integer(y)
  }
  res <- .Call('_mutual_inf_cd_vec', x, y,
               as.integer(k), as.integer(seed), as.logical(use_cc))
  res
}


#'
#' Mutual information estimation for continuous/discrete case.
#'
#' Compute mutual information of \code{x} and \code{y}
#' where \code{x} is continuous vector and \code{y} are discrete
#' @param x input (continuous) vector.
#' @param M input (discrete) matrix. It should have the same number of rows as \code{length(x)}.
#' @param k number of nearest neighbors.
#' @param seed integer random seed. If seed<=0 a random seed is used.
#' @param use_cc (logical) if TRUE the algorithm falls into continuous/continuous case
#' @useDynLib CaDrA _mutual_inf_cd_mat
#'
#' @return vector of length m, where m is the number of columns in M
#' vectors \code{x} and \code{y}.
#' @examples
#'
#' data(mutual_info_df)
#'
#' M <- cbind(mutual_info_df$Xd, mutual_info_df$Yd)
#' mutual_inf_cd_mat(mutual_info_df$Zc_XdYdWd, M)
#' ## 0.1070804 0.1041177
#'
#' @export
mutual_inf_cd_mat <- function(x, M, k=3L, seed=0L, use_cc=FALSE){
  stopifnot( "x and M must have the same length"=length(x) == nrow(M) )
  stopifnot("k must be less than the length of x"=k < length(x))

  if (!is.integer(M)) {
    M <- matrix(as.integer(M), nrow=nrow(M))
  }

  res <- .Call('_mutual_inf_cd_mat', x, M,
               as.integer(k), as.integer(seed), as.logical(use_cc))
  res
}


#'
#' Conditional mutual information estimation for continuous case of 3 vectors.
#'
#' Compute conditional mutual information of \code{x},\code{y} given \code{z}
#' where \code{x}, \code{y} and \code{z} are all continuous
#' @param x input vector.
#' @param y input vector of the same length as x.
#' @param z conditional input vector of the same length as x.
#' @param k number of nearest neighbors.
#' @param seed integer random seed. If seed<=0 a random seed is used.
#' @useDynLib CaDrA _cond_mutual_inf_ccc_vec
#'
#' @return a double-precision value - mutual information estimation for
#' vectors \code{x} and \code{y}.
#'
#' @examples
#' data(mutual_info_df)
#'
#' cond_mutual_inf_ccc_vec(mutual_info_df$Zc_XcYc,
#'                        mutual_info_df$Xc, mutual_info_df$Yc)
#' ## 0.2936858
#'
#'
#' @export
cond_mutual_inf_ccc_vec <- function(x, y, z, k=3L, seed=0L){
  stopifnot( "x and y must have the same length"=length(x) == length(y) )
  stopifnot( "x and y must have the same length"=length(x) == length(z) )
  stopifnot("k must be less than the length of input vectors"=k < length(x))

  res <- .Call('_cond_mutual_inf_ccc_vec', x, y, z, as.integer(k), as.integer(seed))
  res

}


#'
#' Conditional mutual information estimation for continuous/discrete case of 3 vectors.
#'
#' Compute conditional mutual information of \code{x},\code{y} given \code{z}
#' where \code{x} is continuous, \code{y} and \code{z} are discrete
#' @param x input vector.
#' @param y input vector of the same length as x.
#' @param z conditional input vector of the same length as x.
#' @param k number of nearest neighbors.
#' @param seed integer random seed. If seed<=0 a random seed is used.
#' @useDynLib CaDrA _cond_mutual_inf_cdd_vec
#'
#' @return a double-precision value - mutual information estimation for
#' vectors \code{x} and \code{y}, given \code{z}.
#'
#' @examples
#' data(mutual_info_df)
#'
#' cond_mutual_inf_cdd_vec(mutual_info_df$Zc_XdYd, mutual_info_df$Xd,
#'                   mutual_info_df$Yd)
#' ## 0.1338664
#'
#' @export
cond_mutual_inf_cdd_vec <- function(x, y, z, k=3L, seed=0L){
  stopifnot( "x and y must have the same length"=length(x) == length(y) )
  stopifnot( "x and y must have the same length"=length(x) == length(z) )
  stopifnot("k must be less than the length of x"=k < length(x))

  res <- .Call('_cond_mutual_inf_cdd_vec', x, as.integer(y),
               as.integer(z), as.integer(k), as.integer(seed))
  res

}



#'
#' Conditional mutual information estimation for continuous case of a vector and matrix, given matrix.
#'
#' Compute conditional mutual information of vector\code{x}, matrix \code{M}
#' given matrix\code{Z}
#' where \code{x}, \code{M} and \code{Z} are all continuous
#' @param x input vector.
#' @param M input matrix of the same length as x.
#' @param Z conditional input matrix of the same length as x.
#' @param k number of nearest neighbors.
#' @param seed integer random seed. If seed<=0 a random seed is used.
#' @useDynLib CaDrA _cond_mutual_inf_ccc_mat
#'
#' @return a double-precision vector
#'
#' @examples
#' data(mutual_info_df)
#'
#' M <- cbind(mutual_info_df$Xc, mutual_info_df$Yc)
#' ZM <- cbind(mutual_info_df$Yc, mutual_info_df$Wc)
#' cond_mutual_inf_ccc_mat(mutual_info_df$Zc_XcYcWc, M, ZM)
#' ## 0.1171533 0.2192397
#'
#' @export
cond_mutual_inf_ccc_mat <- function(x, M, Z, k=3L, seed=0L){
  stopifnot( "M must be a matrix"= (class(M)[1]=="matrix"))
  stopifnot( "Z must be a matrix"= (class(Z)[1]=="matrix"))
  stopifnot( "x and M must have the same length"=length(x) == nrow(M) )
  stopifnot( "x and Z must have the same length"=length(x) == nrow(Z) )
  stopifnot( "M and Z must be the same size"=dim(M) == dim(Z) )
  stopifnot("k must be less than the length of x"=k < length(x))

  res <- .Call('_cond_mutual_inf_ccc_mat', x, M, Z, as.integer(k), as.integer(seed))
  res

}


#'
#' Conditional mutual information estimation for a continuous vector
#' and a discrete matrix, given another discrete matrix.
#'
#' Compute conditional mutual information of \code{x},\code{M} given \code{Z}
#' where \code{x} is continuous, \code{M} and \code{Z} are discrete
#' @param x input vector.
#' @param M input matrix of the same length as x.
#' @param Z conditional input matrix of the same length as x.
#' @param k number of nearest neighbors.
#' @param seed integer random seed. If seed<=0 a random seed is used.
#' @useDynLib CaDrA _cond_mutual_inf_cdd_mat
#'
#'
#' @return a double-precision vector - mutual information estimation for
#' vector \code{x} and matrix \code{M}, given matrix \code{Z}.
#'
#' @examples
#' data(mutual_info_df)
#'
#' M <- cbind(mutual_info_df$Xd, mutual_info_df$Yd)
#' ZM <- cbind(mutual_info_df$Yd, mutual_info_df$Wd)
#' cond_mutual_inf_cdd_mat(mutual_info_df$Zc_XdYdWd, M, ZM)
#' ## 0.1757598 0.1086227
#'
#' @export
cond_mutual_inf_cdd_mat <- function(x, M, Z, k=3L, seed=0L){

  stopifnot( "M must be a matrix"= (class(M)[1]=="matrix"))
  stopifnot( "Z must be a matrix"= (class(Z)[1]=="matrix"))
  stopifnot( "x and M must have the same length"=length(x) == nrow(M) )
  stopifnot("k must be less than the length of x"=k < length(x))
  stopifnot( "M and Z must be the same size"=dim(M) == dim(Z) )

  if (!is.integer(M)) {
    M <- matrix(as.integer(M), nrow=nrow(M))
  }
  if (!is.integer(Z)) {
    Z <- matrix(as.integer(Z), nrow=nrow(Z))
  }

  if (!is.integer(Z)) {
    Z <- matrix(as.integer(Z), nrow=nrow(Z))
  }
  res <- .Call('_cond_mutual_inf_cdd_mat', x, M, Z, as.integer(k), as.integer(seed))
  res
}


