
#' Top 'N' Plot
#' 
#' Plots a heatmap representation of overlapping features given a list of top N candidate search results
#' @param topN_list a list of lists, where each list entry is one that is returned by the candidate search run for a given starting index. This is computed within and can be returned by the topn.eval() function.
#' 
#' @return a heatmap of the top N evaluation for a given top N search evaluation
#' @examples
#' # Load pre-computed Top-N list generated for sim.ES dataset
#' data(topn.list)
#' topn_plot(topn.list)
#' 
#' @export
#' @import gplots
#' @importFrom graphics legend
topn_plot <- function(topN_list){
  
  eset_l <- lapply(topN_list, "[[", 1)
  scores_l <- lapply(topN_list, "[[", 2)
  
  f_list <- lapply(eset_l,featureNames)  #Get the list of feature names from each ESet
  
  f_union <- Reduce(f = union,f_list) #Get the union of all features that were returned across all top N runs
  
  f_checklist<-lapply(f_list,function(x,ref=f_union){
    return(f_union %in% x)
  })
  
  # Make a matrix indicating which features are found across each top n run
  m <- do.call(cbind,f_checklist)*1   #Multiplying by 1 is just to convert boolean values into 1's and 0's
  rownames(m) <- f_union
  
  # Working with scores for each top N run
  s <- unlist(scores_l)
  colnames(m) <- names(s)
  
  # Order matrix in increasing order of KS score p-values
  # Add labels of which rank it was originally, and what the meta-feature p-value is
  # Here we take the negative log transform of the p-value just to avoid 0s (if p-values are too small)
  # Note that this means the HIGHER the transformed score, the more significant
  s.log <- -log(s) 
  
  colnames(m) <- paste(colnames(m)," [",seq(1,ncol(m)),"] ",round(s.log,3),sep="")
  m <- m[,order(s.log,decreasing = T)] #We order matrix columns in increasing order of search p-value (i.e. decreasing negative-log p-value)
  
  colcode <- if (all(m==1)) c("firebrick2","white") else c("white","firebrick2")
  
  if(ncol(m)>=2){
    
    cat("Generating top N overlap heatmap..\n\n")
    heatmap.2(x = m,
              col=colcode,
              Colv=FALSE,
              dendrogram="none",
              margins=c(10,10),
              cexRow=0.7,
              cexCol=0.7,
              cex.main=0.8,
              key=F,
              trace="none",
              sepwidth=c(0.1,0.1),
              sepcolor="grey90",
              colsep=1:ncol(m),
              rowsep=1:nrow(m))
    
    legend("topleft",
           legend=c("Present","Absent"),
           fill=c("firebrick2","white"),
           bty="n")
    
  } else{
    warning("Cannot plot overlap matrix for N=1. Please use a larger N value for top N evaluation visualization..\n\n")
  }
  
}


