#' kde2d_wrap_ wrapper
#'
#' Simplified version of Two-Dimensional Kernel Density Estimation 
#' @param n length 
#' @param x x coordinate of data 
#' @param y y coordinate of data
#' @param h vector of bandwidths for x and y directions.
#' @param lims The limits of the rectangle covered by the grid as c(xl, xu, yl, yu).
#' @useDynLib CaDrA kde2d_wrap_
kde2d_wrap <- function(x, y, h, n=25, lims = c(range(x), range(y))) {

  res <- .Call( kde2d_wrap_, as.integer(n), x, y, h, lims)
  res
}