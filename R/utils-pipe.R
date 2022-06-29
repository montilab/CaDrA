#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @examples 
#' 
#' # Load R library
#' library(Biobase)
#' 
#' # Load pre-computed expression set
#' data(sim.ES)
#' 
#' # filter mtcars with mpg > 20
#' ESet <- sim.ES %>% exprs(.)
#' 
#' @return The result of calling `rhs(lhs)`.
NULL
