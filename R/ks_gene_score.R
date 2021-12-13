
#' Compute directional KS scores for each row of a given binary matrix
#'
#' @param x ranked list
#' @param y positions of geneset items in ranked list (ranks)
#' @param weight weights for weighted score (see Subramanian et al.) (usually, sort(score))
#' @param weight_p weights exponent
#' @param alternative alternative hypothesis for p-value calculation
#' @param do_pval compute asymptotic p-value
#' @param absolute takes max - min score rather than the maximum deviation from null
#' @param exact compute exact p-value
#'
#' @return A data frame
#' @export
#' @importFrom stats ks.test
ks_gene_score <- function
(
  x,                                             
  y,                                             
  weight = NULL,
  weight_p = 1,
  alternative = "less",  
  do_pval = TRUE, 
  absolute = FALSE, 
  exact = NULL
)
{
  
  # Get the length of ranked list
  n.x = length(x);                            
  
  # efficient version of ks.score (should give same results as ks.test, when weight=NULL)
  #
  alternative <- match.arg(alternative, c("two.sided", "greater", "less"))
  DNAME <- paste( "1:", n.x, " and ", deparse(substitute(y)), sep="" )
  METHOD <- "Two-sample Kolmogorov-Smirnov test"
  n.y <- length(y)
  if ( n.y < 1 )  stop("Not enough y data")
  if ( any(y>n.x) ) stop( "y must be <= n.x: ", max(y) )
  if ( any(y<1) ) stop( "y must be positive: ", min(y) )
  if ( do_pval && !is.null(weight) ) warning("p-value meaningless w/ weighted score")
  if ( !is.null(weight) && length(weight)!=n.x ) stop("weights must be same length as ranked list: ", length(weight), " vs ", n.x)
  x.axis <- y.axis <- NULL
  
  # weighted GSEA score
  #
  if ( !is.null(weight) )
  {
    weight <- abs(weight[y])^weight_p
    
    Pmis <- rep(1, n.x); Pmis[y] <- 0; Pmis <- cumsum(Pmis); Pmis <- Pmis/(n.x-n.y)
    Phit <- rep(0, n.x); Phit[y] <- weight; Phit <- cumsum(Phit); Phit <- Phit/Phit[n.x]
    z <- Phit-Pmis
    
    score <- if (absolute) max(z)-min(z) else z[which.max(abs(z))]
    
    x.axis <- 1:n.x
    y.axis <- z
  }
  # KS score
  #
  else
  {
    y <- sort(y)
    n <- n.x * n.y/(n.x + n.y)
    hit <- 1/n.y
    mis <- 1/n.x
    
    # to compute score, only the y positions and their immediate preceding
    # ..positions are needed
    #
    Y <- sort(c(y-1,y)); Y <- Y[diff(Y)!=0]; y.match <- match(y,Y); if ( any(is.na(y.match)) ) browser()
    D <- rep( 0, length(Y) ); D[y.match] <- (1:n.y)
    zero <- which(D==0)[-1]; D[zero] <- D[zero-1]
    
    z <- D*hit - Y*mis
    
    score <- if (absolute) max(z)-min(z) else z[which.max(abs(z))]
    
  }
  
  # Here, choose suppressWarnings simply because you will generally have ties for binary data matrix
  PVAL <- suppressWarnings(ks.test(x=1:n.x, y=y, alternative=alternative, exact=exact)$p.value)
  
  return(data.frame(score=score, p_value=PVAL))
  
}
