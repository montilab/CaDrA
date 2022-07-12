#' Compute a Kolmogorov-Smirnov score for a given ranked list
#'
#' @param n.x length of a ranked list
#' @param y positions of interested geneset items in the ranked list 
#' @param weights weights for weighted score (see Subramanian et al.) (usually, sort(score))
#' @param weight_p weights exponent
#' @param alternative alternative hypothesis for p-value calculation
#' @param do_pval compute asymptotic p-value
#' @param absolute takes max - min score rather than the maximum deviation from null
#' @param exact compute exact p-value
#' @param plot_dat return dataframe of x and y axis data for ES plot
#' 
#' @examples
#' 
#' # Load R library
#' library(Biobase)
#' 
#' # Load pre-computed Top-N list generated for sim.ES dataset
#' data(topn.list)
#' 
#' # With the results obtained from top-N evaluation,
#' # We can find the combination of features that gives the best score in top N searches
#' topn_best_meta <- topn_best(topn_list=topn.list) 
#' 
#' # Extract the meta-feature set
#' ESet =  topn_best_meta[["ESet"]]                    
#' 
#' # Make sure mat variable is a matrix
#' mat <- as.matrix(exprs(ESet))
#' 
#' # For cases when matrix only has one row
#' if(ncol(mat) == 1){
#'   mat <- matrix(t(mat), nrow=1, byrow=TRUE, dimnames = list(rownames(ESet), rownames(mat))) 
#' }
#' 
#' # Add on the OR function of all the returned entries
#' or <- ifelse(colSums(mat)==0, 0, 1)
#' 
#' mat <- rbind(mat, or)
#' 
#' # Get x and y axis data for ES plot of 
#' # cumulative function of individual features (i.e. the OR function)
#' ES_dat <- ks_gene_score(
#'    n.x=length(or), y=which(or==1), plot_dat = TRUE, alternative = "less"
#' )
#'
#' @return a data frame with two columns: \code{score} and \code{p_value}
#' @export
#' @importFrom stats ks.test
ks_gene_score <- function
(
  n.x,                                             
  y,                                             
  weights = NULL,
  weight_p = 1,
  alternative = "less",  
  do_pval = TRUE, 
  absolute = FALSE, 
  exact = NULL,
  plot_dat = FALSE
)
{
  # efficient version of ks.score (should give same results as ks.test, when weights=NULL)
  #
  alternative <- match.arg(alternative, c("two.sided", "greater", "less"))
  DNAME <- paste( "1:", n.x, " and ", deparse(substitute(y)), sep="" )
  METHOD <- "Two-sample Kolmogorov-Smirnov test"
  n.y <- length(y)
  if ( n.y < 1 )  stop("Not enough y data")
  if ( any(y > n.x) ) stop( "y must be <= n.x: ", max(y) )
  if ( any(y < 1) ) stop( "y must be positive: ", min(y) )
  if ( do_pval && !is.null(weights) ) warning("p-value is meaningless w/ weighted score")
  if ( !is.null(weights) && length(weights) != n.x ) stop("weights must be same length as ranked list: ", length(weights), " vs ", n.x)
  x.axis <- y.axis <- NULL
  
  # weighted GSEA score
  #
  if ( !is.null(weights) )
  {
    weights <- abs(weights[y])^weight_p
    
    Pmis <- rep(1, n.x); Pmis[y] <- 0; Pmis <- cumsum(Pmis); Pmis <- Pmis/(n.x-n.y)
    Phit <- rep(0, n.x); Phit[y] <- weights; Phit <- cumsum(Phit); Phit <- Phit/Phit[n.x]
    z <- Phit-Pmis
    
    score <- if (absolute) max(z)-min(z) else z[which.max(abs(z))]
    
    x.axis <- seq_len(n.x)
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
    Y <- sort(c(y-1,y)); Y <- Y[diff(Y)!=0]; y.match <- match(y,Y); #if ( any(is.na(y.match)) ) browser()
    D <- rep( 0, length(Y) ); D[y.match] <- seq_len(n.y)
    zero <- which(D==0)[-1]; D[zero] <- D[zero-1]
    
    z <- D*hit - Y*mis
    
    score <- if (absolute) max(z)-min(z) else z[which.max(abs(z))]
    
    if (plot_dat) {
      x.axis <- Y;
      y.axis <- z;
      if(Y[1] > 0) {
        x.axis <- c(0, x.axis);
        y.axis <- c(0, y.axis);
      }
      if ( max(Y) < n.x ) {
        x.axis <- c(x.axis,n.x)
        y.axis <- c(y.axis,0)
      }
    }
  }
  if(plot_dat){
    d <- data.frame("x"=x.axis, "y"=y.axis)
    return(d)
  }
  # Here, choose suppressWarnings simply because you will generally have ties for binary data matrix
  PVAL <- ks.test(x=seq_len(n.x), y=y, alternative=alternative, exact=exact)$p.value
  
  return(c(score=score, p_value=PVAL))
}