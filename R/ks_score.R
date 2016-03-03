#' Row-wise matrix Kolmogorov-Smirnov scoring
#' 
#' Compute directional KS scores for each row of a given matrix
#' @param mat matrix of binary features to compute row-wise ks scores for
#' @param alt a character string specifying the alternative hypothesis, must be one of "two.sided","greater" or "less". Value passed to ks.genescore() function
#' @param weight a vector of weights to use if performing a weighted-KS test. Default is NULL. Value passed to ks.genescore() function 
#' @return A 2 x N matrix (N = number of rows in matrix) with the first and second rows corresponding to the KS scores and p-values, respectively 
#' @export
ks.genescore.mat<-function(mat, 
                           alt="two.sided", #One of "two.side","greater" or "less" for alternative hypotheses p-value computation. Value passed to ks.genescore() function
                           weight=NULL #If performing a weighted-KS test, requires ordered weights (vector of values) per sample. Value passed to ks.genescore() function
){
  
  #Compute the ks statitic and p-value per row in the matrix
  ks<-apply(X = mat, 1, FUN = function(x,w=weight){
    ks.genescore(n.x=length(x),y=which(x==1),alternative=alt,weight=w,bare=T)})
  #Bare=T will simply return a 'two-tuple' containing the ks statistic and the p-value
  #Applying this across the matrix will create a two row matrix, with the first the statistic and the second the p-value
  
  return(ks)
}

#' Directional Kolmogorov-Smirnov scoring
#'
#' Compute two-sample directional KS scores 
#' @param n.x length of ranked list
#' @param y positions (ranks) of geneset items in ranked list
#' @param do.pval a logical indicating whether or not to compute asymptotic p-value or not
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided","greater" or "less". 
#' @param do.plot a logical indicating whether or not to draw the Enrichment Score (ES) plot
#' @param plot.dat a logical indicating whether or not to return the x and y axis data for the ES plot
#' @param bare a logical indicating whether or not to return the KS score & p-value only (a 2-tuple)
#' @param weight a vector of weights to use if performing a weighted-KS test. Default is NULL. 
#' @param weight.p Exponent for weights
#' @param cls.lev class labels to display
#' @param absolute logical indicating whether or not to take the (max - min) score rather than the maximum deviation from null
#' @param plot.labels logical indicating whether or not to plot hits' labels
#' @param exact 
#' @param ... additional arguments for plot() function
#' @return If plot.dat is set to TRUE, this function returns a data frame object containing x and y-axis data corresponding to the ES plot. If bare is set to TRUE, this function returns a tuple vector containing the directional ks score and p-value, respectively. If both are set to FALSE, this function returns a list object containing values corresponding to the ks statistic, p-value, alternative, method and data name 
#' @export
ks.genescore <- function
(
  n.x,               # length of ranked list
  y,                 # positions of geneset items in ranked list (basically, ranks)
  do.pval=T,         # compute asymptotic p-value
  alternative=c("two.sided","greater","less"),
  do.plot=F,         # draw the ES plot
  plot.dat=F,        # return dataframe of x and y axis data for ES plot
  bare=F,            # return score & p-value only (a 2-tuple)
  weight=NULL,       # weights for weighted score (see Subramanian et al.) (usually, sort(score))
  weight.p=1,        # weights' exponent
  cls.lev=c(0,1),    # class labels to display
  absolute=F,        # takes max - min score rather than the maximum deviation from null
  plot.labels=FALSE, # hits' labels
  exact=NULL,
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
  if ( do.plot )
  {
    plot( x.axis, y.axis, type="l",
          xlab=paste("up-regulated for class ", cls.lev[2], " (KS>0) vs ",
                     "up-regulated for class ", cls.lev[1], " (KS<0)", sep="" ),
          ylab="gene hits",...)
    abline(h=0)
    abline(v=n.x/2,lty=3)
    axis(1,at=y,labels=plot.labels,tcl=0.25,las=2)
    i.max <- which.max(abs(y.axis))
    points( x.axis[i.max], y.axis[i.max], pch=20, col="red")
    text(x.axis[i.max]+n.x/20,y.axis[i.max],round(y.axis[i.max],2))
  }
  if ( !do.pval )
    return(score)
  
  # ELSE, compute asymptotic p-value
  #
  names(score) <- switch(alternative, two.sided="D", greater="D^+", less="D^-")
  PVAL <- ks.test(1:n.x,y=y,alternative=alternative,exact=exact)$p.value
  
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