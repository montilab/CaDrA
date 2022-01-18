#' @useDynLib CaDrA ks_genescore_mat_
ks_genescore_mat <- function(mat,alt, weight) .Call(ks_genescore_mat_, mat, alt, weight)


#' Compute skewdness scores per feature
#'
#' Compute scores based on skewdness of a given binary matrix (based on sample ordering) used in the stepwise search
#' @param mat matrix of binary features to compute row-wise ks scores for
#' @param method a character string specifying the method used to score features, must be one of "ks" or "wilcox"
#' @param ... additional arguments provided to either ks.genescore() or wilcox.genescore() functions depending on method of choice
#' @export
compute_score <- function
(
mat,                     # Matrix to apply scoring function over
method=c("ks","wilcox"), # Scoring method to apply over each row in matrix
...                      # Additional parameters passed to either scoring function (Ex: alternative, ranks, weights etc.)
)
{
  # Additional arguments for ks or wilcox.genescore.mat functions
  compute_score_args <- list(...)
  
  if(nrow(mat) < 2)
    warning("You are computing a row-wise statistic over a matrix with nrow < 2\n")
  
  # If no alternative is specified, we use "less" as default.
  if(is.null(compute_score_args$alt)){
    warning("No alternative hypothesis specified. Using 'less' by default ..\n")
    compute_score_args$alt <- "less"
  }
  
  # If invalid method specification
  if(!(method %in% c("ks","wilcox")))
     stop("Invalid method specification to compute scores.. please specify either 'ks' or 'wilcox'..\n")
   
  # Check if there are rows with ALL 0's or 1's (cannot use wilcox or KS test if this is the case)
  # This may occur when taking unions of meta-feature with existing features (feature 'saturation') where they all become 1 upon taking the union
  if(any(rowSums(mat)==ncol(mat))){
    num_1_f <- sum(rowSums(mat)==ncol(mat))
    stop(num_1_f, " entries with all 1's present.. Cannot compute statistics for such features!\n")
  } 
  
  if (any(rowSums(mat)==0)){
    num_0_f <- sum(rowSums(mat)==0)
    stop(num_0_f, " entries with all 0's present.. Cannot compute statistics for such features!\n") 
  }
  
  s <- switch(method,
                 ks=ks_genescore_mat(mat,
                                      alt=compute_score_args$alt, 
                                      weight=compute_score_args$wts),
                 #ks = ks_genescore_mat(mat),
                 wilcox=wilcox.genescore.mat(mat,
                                             alt=compute_score_args$alt,
                                             ranks=compute_score_args$rnks)) 
     
  # Score returned by either ks or wilcox-based functions
                 
  return(s) 
}

#' Row-wise matrix Kolmogorov-Smirnov scoring
#' 
#' Compute directional KS scores for each row of a given binary matrix
#' @param mat matrix of binary features to compute row-wise scores for based on the  Kolmogorov-Smirnov test
#' @param alt a character string specifying the alternative hypothesis, must be one of "two.sided","less" or "greater". Value passed to ks.genescore() function
#' @param weight a vector of weights to use if performing a weighted-KS test. Default is NULL. Value passed to ks.genescore() function 
#' @return A 2 x N matrix (N = number of rows in matrix) with the first and second rows corresponding to the KS scores and p-values, respectively 
#' @export
ks.genescore.mat<-function
(
  mat, 
  alt=c("two.sided","less","greater"), 
  weight=NULL 
)
{
  
  #Compute the ks statitic and p-value per row in the matrix
  ks<-apply(X = mat, 1, FUN = function(x,w=weight){
    ks.genescore(n.x=length(x),y=which(x==1),alternative=alt,weight=w,bare=TRUE)
    })
  #bare=T will simply return a 'two-tuple' containing the ks statistic and the p-value
  #Applying this across the matrix will create a two row matrix, with the first the statistic and the second the p-value
  
  return(ks)
}

#' Row-wise matrix Wilcoxon rank sum scoring
#' 
#' Compute rank sum scores for each row of a given binary matrix
#' @param mat matrix of binary features to compute row-wise scores for based on the Wilcoxon rank sum test
#' @param alt a character string specifying the alternative hypothesis, must be one of "two.sided","less" or "greater". Value passed to wilcox.genescore() function
#' @param ranks a vector of ranks to use when performing the Wilcoxon test. Default is NULL. If NULL, then samples are assumed to be ordered by increasing ranking. Value passed to wilcox.genescore() function 
#' @return A 2 x N matrix (N = number of rows in matrix) with the first and second rows corresponding to the rank sum statistic scores and p-values, respectively 
#' @export
wilcox.genescore.mat <- function
(
  mat,
  alt=c("two.sided","less","greater"), 
  ranks=NULL
)
{
  
  # If ranks for samples are not provided, assume it's ordered by decreasing rank and assign ranks 1: N (N: number of samples)
  if(is.null(ranks))
    ranks <- seq(1,ncol(mat))
  
  #Compute the wilcox rank sum statitic and p-value per row in the matrix
  w <- apply(X = mat, 1, FUN = function(x,r=ranks){
    # Note that the first and second group being compared are the ones with and without the alteration, respectively. 
    # Alternative should be set to "less" if doing a one-tailed test, in this case to find features that are left-skewed (lower ranked) 
    wilcox.genescore(x=r[x==1],y=r[x==0],alternative=alt,bare=TRUE,exact=FALSE)
  })
  
  #bare=TRUE will simply return a 'two-tuple' containing the statistic and the p-value
  #Applying this across the matrix will create a two row matrix, with the first the statistic and the second the p-value (if bare=TRUE)
  
  return(w)
}

#' @importFrom stats ks.test
ks.genescore <- function
(
  n.x,                                         # length of ranked list
  y,                                           # positions of geneset items in ranked list (ranks)
  do.pval=TRUE,                                # compute asymptotic p-value
  alternative=c("two.sided","less","greater"), # alternative hypothesis for p-value calculation
  do.plot=FALSE,                               # draw the ES plot
  plot.dat=FALSE,                              # return dataframe of x and y axis data for ES plot
  bare=FALSE,                                  # return score & p-value only (a 2-tuple)
  weight=NULL,                                 # weights for weighted score (see Subramanian et al.) (usually, sort(score))
  weight.p=1,                                  # weights' exponent
  cls.lev=c(0,1),                              # class labels to display
  absolute=FALSE,                              # takes max - min score rather than the maximum deviation from null
  plot.labels=FALSE,                           # hits' labels
  exact=NULL,                                  # compute exact p-value
  ...                # additional plot arguments
)
{
  # efficient version of ks.score (should give same results as ks.test, when weight=NULL)
  #
  alternative <- match.arg(alternative)
  DNAME <- paste( "1:", n.x, " and ", deparse(substitute(y)), sep="" )
  METHOD <- "Two-sample Kolmogorov-Smirnov test"
  n.y <- length(y)
  if ( n.y < 1 )  stop("Not enough y data")
  if ( any(y>n.x) ) stop( "y must be <= n.x: ", max(y) )
  if ( any(y<1) ) stop( "y must be positive: ", min(y) )
  if ( do.pval && !is.null(weight) ) warning("p-value meaningless w/ weighted score")
  if ( !is.null(weight) && length(weight)!=n.x ) stop("weights must be same length as ranked list: ", length(weight), " vs ", n.x)
  x.axis <- y.axis <- NULL
  
  # weighted GSEA score
  #
  if ( !is.null(weight) )
  {
    weight <- abs(weight[y])^weight.p
    
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
    
    if (do.plot) {
      x.axis <- Y;
      y.axis <- z;
      if(Y[1]>0) {
        x.axis <- c(0,x.axis);
        y.axis <- c(0,y.axis);
      }
      if ( max(Y)<n.x ) {
        x.axis <- c(x.axis,n.x)
        y.axis <- c(y.axis,0)
      }
    }
  }
  
  if ( plot.dat )
  {
    d<-data.frame("x"=x.axis,"y"=y.axis)
    return(d)
  }
# if ( do.plot )
#  {
#    plot( x.axis, y.axis, type="l",
#          xlab=paste("up-regulated for class ", cls.lev[2], " (KS>0) vs ",
#                     "up-regulated for class ", cls.lev[1], " (KS<0)", sep="" ),
#          ylab="gene hits",...)
#    abline(h=0)
#    abline(v=n.x/2,lty=3)
#    axis(1,at=y,labels=plot.labels,tcl=0.25,las=2)
#    i.max <- which.max(abs(y.axis))
#    points( x.axis[i.max], y.axis[i.max], pch=20, col="red")
#    text(x.axis[i.max]+n.x/20,y.axis[i.max],round(y.axis[i.max],2))
#  }
  if (!do.pval)
    return(score)
  
  # ELSE, compute asymptotic p-value
  #
  names(score) <- switch(alternative, two.sided="D", greater="D^+", less="D^-")
  # Here, choose suppressWarnings simply because you will generally have ties for binary data matrix
  PVAL <- suppressWarnings(ks.test(1:n.x,y=y,alternative=alternative,exact=exact)$p.value)
  
  if ( bare ) {
    return( c(score=score, p.value=PVAL) )
  }
  RVAL <- list(statistic = score,
               p.value = PVAL, alternative = alternative, 
               method = METHOD, data.name = DNAME)
  class(RVAL) <- "htest"
  
  if(!plot.dat)
    return(RVAL)
}

#' @importFrom stats pnorm pwilcox
wilcox.genescore <- function 
(x,                                                  # ranks for group 1
 y,                                                  # ranks for group 2
 alternative = c("two.sided", "less", "greater"),    # alternative hypothesis for p-value calculation
 paired = FALSE,                                     # paired test
 exact = NULL,                                       # comput exact p-value
 correct = TRUE,                                     # for contuity correction (p-value)
 bare = FALSE,                                       # only return statistic and p-value tuple
 mu=0,                                               # a number specifying an optional parameter used to form null hypothesis
 ...                                                 # additional arguments passed to base wilcox.test() function
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
    print(x)
    stop("not enough (finite) 'x' observations")}
  
  if (length(y) < 1L){
    print(y)
    stop("not enough 'y' observations")}
  
  METHOD <- "Wilcoxon rank sum test"
  
  ##### Modification ######
  # Take input as ranks instead of continuous measures (normally internally ranked: see below)
  r <- c(x,y)
  
  #r <- rank(c(x - mu, y))


n.x <- as.double(length(x))
n.y <- as.double(length(y))
if (is.null(exact)) 
  exact <- (n.x < 50) && (n.y < 50)
STATISTIC <- c(W = sum(r[seq_along(x)]) - n.x * (n.x + 
                                                   1)/2)
TIES <- (length(r) != length(unique(r)))
if (exact && !TIES) {
  PVAL <- switch(alternative, two.sided = {
    p <- if (STATISTIC > (n.x * n.y/2)) pwilcox(STATISTIC - 
                                                  1, n.x, n.y, lower.tail = FALSE) else pwilcox(STATISTIC, 
                                                                                                n.x, n.y)
    min(2 * p, 1)
  }, greater = {
    pwilcox(STATISTIC - 1, n.x, n.y, lower.tail = FALSE)
  }, less = pwilcox(STATISTIC, n.x, n.y))
} else {
  NTIES <- table(r)
  z <- STATISTIC - n.x * n.y/2
  SIGMA <- sqrt((n.x * n.y/12) * ((n.x + n.y + 1) - 
                                    sum(NTIES^3 - NTIES)/((n.x + n.y) * (n.x + n.y - 
                                                                           1))))
  if (correct) {
    CORRECTION <- switch(alternative, two.sided = sign(z) * 
                           0.5, greater = 0.5, less = -0.5)
    METHOD <- paste(METHOD, "with continuity correction")
  }
  z <- (z - CORRECTION)/SIGMA
  PVAL <- switch(alternative, less = pnorm(z), greater = pnorm(z, 
                                                               lower.tail = FALSE), two.sided = 2 * min(pnorm(z), 
                                                                                                        pnorm(z, lower.tail = FALSE)))
  
  if (exact && TIES) 
    warning("cannot compute exact p-value with ties")
}
names(mu) <- ifelse (paired || !is.null(y),"location shift","location") 

RVAL <- list(statistic = STATISTIC, parameter = NULL, p.value = as.numeric(PVAL), 
             null.value = mu, alternative = alternative, method = METHOD, 
             data.name = DNAME)

class(RVAL) <- "htest"

if(bare){	  
   return(c(score=RVAL$statistic, p.value=RVAL$p.value))         
} else{
  RVAL
}
} #end function