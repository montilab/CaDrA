# Function to convert units depending on the version of R (for adjusting the dimensions of the ggplots)
.setup.tab.heights <- function(plot.tab, dims) {
  panels <- plot.tab$layout$t[grep("panel", plot.tab$layout$name)]
  has.unit.vec <- getRversion() > '3.3'
  if(has.unit.vec) {
    plot.tab$heights[panels] <- unit(dims, units="null")
  }
  else {
    plot.tab$heights[panels] <- lapply(dims, unit, "null") 
  }
  return(plot.tab)
}


# Solution to drawing arbitrary number of plots stacked above one another
# Each plot will be left-aligned
# Each plot should first be made a ggplotGrob object
# Ex: grid.draw(rbind_gtable_max(ggplotGrob(ggplot1),ggpotGrob(ggplot2)))
#' @import gtable
#' @importFrom grid unit.pmax
rbind_gtable_max <- function(...){
  
  gtl <- list(...)
  stopifnot(all(sapply(gtl, is.gtable)))
  bind2 <- function (x, y) 
  {
    stopifnot(ncol(x) == ncol(y))
    if (nrow(x) == 0) 
      return(y)
    if (nrow(y) == 0) 
      return(x)
    y$layout$t <- y$layout$t + nrow(x)
    y$layout$b <- y$layout$b + nrow(x)
    x$layout <- rbind(x$layout, y$layout)

    # Katia This line causes a warning during the package check() and
    # and considered to be not a good practice
    # See the note in ?`:::` about the use of this operator
    #x$heights <- gtable:::insert.unit(x$heights, y$heights)
    x$heights <- grid::unit.c(x$heights, y$heights)
    x$rownames <- c(x$rownames, y$rownames)
    x$widths <- grid::unit.pmax(x$widths, y$widths)
    x$grobs <- append(x$grobs, y$grobs)
    x
  }
  
  Reduce(bind2, gtl)
}


#' Enrichment Score (ES) plot
#' 
#' Plot the ES running sum statistic, including the max score (KS) using x and y-axis data returned by ks.genescore() when plot.dat is set to TRUE for a given dataset
#' @param d a data frame object with columns 'x' and 'y' containing the necessary data for an ES plot (returned by ks.genescore() when plot.dat is set to TRUE for a given dataset)
#' @return A plot graphic of the Enrichment Score (ES) for a given distribution
#' @export
#' @import ggplot2
plot_ES <- function(d 
){
  #Katia: adding ".data" to avoid a warning during check:
  # no visible binding for global variable 
  g <- ggplot(data = d,aes(x=.data$x, y=.data$y))
  g <- g+
    #geom_line(size=1.25,colour="blueviolet")+
    geom_line(size=1.25,colour="darkgoldenrod1")+
    geom_hline(yintercept=0,linetype=2)+
    #Katia: adding ".data" to avoid a warning during check:
    # no visible binding for global variable
    geom_point(data=d[which.max(d$y),],aes(x=.data$x,y=.data$y),colour="black",fill="red",size=3,shape=21)+
    annotate("text",x=d[which.max(d$y),1]+8,y=max(d$y),label=as.character(round(max(d$y),3)),size=3)+
    scale_x_continuous(expand = c(0,0))+
    scale_y_reverse()+ # This inverts the ES score statistic line
    theme_classic()+
    theme(panel.border=element_rect(colour="black",fill=NA))+
    labs(x="",y="Enrichment Score (ES)")
  
  g
}


#' Step-wise heuristic search plot
#' 
#' Plot the stepwise search results for a given CaDrA run. Plot will include an optional bar plot of the continuous ranking variable (top), 
#' @param ESet an ESet object containing only the features returned from the stepwise.search() function (if any) 
#' @param var.score an optional integer vector of continuous measures used to rank samples for the stepwise search (assumed in matching order)
#' @param var.name a string object describing the name of the continuous measure used, which will be used as the y-axis label for the metric plot
#' @return A plot graphic with the ranked metric plot (optional), a tile plot of the features within the provided ESet, and the corresponding Enrichment Score (ES) for a given distribution (here, this will correspond to the logical OR of the features)
#' @examples
#' data(sim.ES)
#' data(topn.list)
#' 
#' # Plot the results from a top-N evaluation by passing the resulting ESet from a specific run
#' # To find the combination of features that had the best score
#' best.meta <- topn.best(topn.list) 
#' 
#' # Now we can plot this set of features
#' meta.plot(best.meta$ESet)
#' 
#' # If a continuous ranking variable was used for the sample-ranking, we can visualize it together
#' # Just for illustation purposes, we simulate random sample scores
#' sample.scores <- sort(runif(ncol(sim.ES),0,1),decreasing=TRUE)
#' meta.plot(best.meta$ESet,var.score=sample.scores,var.name="My ranking variable")
#' @export
#' @import ggplot2 reshape2
#' @importFrom grid unit.pmax grid.draw
meta.plot<-function(ESet, #ExpressoinSet containing somatic mutation/CNA data with samples ordered by a given measure
                    var.score=NULL, #Ordered vector of measure used for search
                    var.name=""){
  
  #Plot y axis label
  y_lab <- var.name
  #Plot x axis label
  x_lab <- paste("Samples (n = ",ncol(ESet),")",sep="")
  
  
  # Plot for continuous metric used to rank samples (Ex: ASSIGN scores)
  if(!is.null(var.score)){
    
    # Make dataframe of metric and sample information, with rows arranged in specific order of measure used for search
    # Here we assume that var.score provided is ordered in the order that the stepwise search was run for the dataset
    var.d<-data.frame("measure"=var.score,"sample"=factor(colnames(ESet),levels=colnames(ESet)))
    
    
    # This plot assumes that there are two columns in the data frame called 'sample' and 'measure'
    # The 'sample' variable has to be a factor with ordered levels for displaying in the correct sample order
    
    #Katia: adding ".data" to avoid a warning during check:
    # no visible binding for global variable 
    m.plot<-ggplot(var.d,aes(x=.data$sample,y=.data$measure,group=1))+
      geom_area(alpha=0.6,fill="deepskyblue4",linetype=1,size=0.5,color="black")+
      scale_y_continuous(expand = c(0,0))+
      scale_x_discrete(expand = c(0,0))+
      theme_classic()+labs(x=x_lab,y=y_lab)+
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.line.x=element_line(color="black"),
            axis.line.y=element_line(color="black"))
  }
  
  
  # Start working with the feature set for plots
  # We are going to fetch the logical OR summary for the set of feature
  # Along with the corresponding ES running score statistic for the OR summary
  mat<-exprs(ESet)
  
  # Add on the OR function of all the returned entries
  or<-ifelse(colSums(mat)==0,0,1)
  mat<-rbind(mat,or)
  
  # Get x and y axis data for ES plot of cumulative function of individual features (i.e. the OR function)
  ES.dat <- ks.genescore(n.x=length(or),y=which(or==1),plot.dat = TRUE, do.plot = TRUE,alternative = "less")
  
  # Plot for ES scores
  ES.plot <- plot_ES(d = ES.dat)
  
  # Give the last row no row name (this is just for the purpose of the plot)
  rownames(mat)[nrow(mat)] <- ""
  
  
  mat[nrow(mat),]<-2*(mat[nrow(mat),]) #Make the OR function have higher values for a different color (red)
  
  
  m <- ExpressionSet(assayData = mat,featureData = AnnotatedDataFrame(data.frame("Name"=rownames(mat),row.names = rownames(mat))))
  
  
  x<-exprs(m)
  
  x.m<-melt(x)
  
  # Make factor for genes so it stays in order of dataset (feature ESet)
  # Here we need to reverse the order of the levels because geom_tile plots ascending from bottom to top
  x.m$Var1<-factor(as.character(x.m$Var1),levels=rev(unique(as.character(x.m$Var1))))
  
  # Plot for mutation/CNA feature profile (with summary OR)
  #Katia: adding ".data" to avoid a warning during check:
  # no visible binding for global variable 
  feature.plot <- ggplot(x.m,aes(x=factor(.data$Var2),y=.data$Var1,fill=.data$value))+
    geom_tile(colour=NA)+ #So that there are no line borders around each cell
    scale_fill_gradient2( #This is only if we want different colors for different values
      high="red", #This will be for our OR function values (which have 2's for 1's, based on our pre-processing) 
      mid="black", #This is for our individual feature 1 values
      low = "white", #This is for the background (0) values
      midpoint = 1)+ 
    theme(line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_text(size=8,colour="black",face="bold"),
          panel.border=element_rect(colour="black",fill=NA))+
    labs(x="",y="Feature") 
  
  
  if(!is.null(var.score)){
    # Align the three plots, adjusting the widths to match
    plot.tab <- rbind_gtable_max(ggplotGrob(m.plot),ggplotGrob(feature.plot+theme(legend.position="none")),ggplotGrob(ES.plot))
    
    # Set heights for each panel
    plot.tab <- .setup.tab.heights(plot.tab, c(2,1.25,3))
    
  } else{
    
    plot.tab <- rbind_gtable_max(ggplotGrob(feature.plot+theme(legend.position="none")),ggplotGrob(ES.plot))
    
    # Set heights for each panel
    plot.tab <- .setup.tab.heights(plot.tab, c(1.25, 3))
  }
  
  
  # Draw combined plot to canvas
  grid.draw(plot.tab)
  
}

#' Top 'N' plot
#' 
#' Plots a heatmap representation of overlapping features given a list of top N stepwise search results
#' @param topN.list a list of lists, where each list entry is one that is returned by the stepwise search run for a given starting index (See ks.stepwise()). This is computed within, and can be returned by the topn.eval() function.
#' @return a heatmap of the top N evaluation for a given top N search evaluation
#' @examples
#' # Load pre-computed Top-N list generated for sim.ES dataset
#' data(topn.list)
#' topn.plot(topn.list)
#' @export
#' @import gplots
#' @importFrom graphics legend
topn.plot <- function(topN.list){
  
  eset.l <- lapply(topN.list, "[[", 1)
  scores.l <- lapply(topN.list, "[[", 2)
  
  
  f_list <- lapply(eset.l,featureNames)  #Get the list of feature names from each ESet
  
  f_union <- Reduce(f = union,f_list) #Get the union of all features that were returned across all top N runs
  
  f_checklist<-lapply(f_list,function(x,ref=f_union){
    return(f_union %in% x)
  })
  
  # Make a matrix indicating which features are found across each top n run
  m <- do.call(cbind,f_checklist)*1   #Multiplying by 1 is just to convert boolean values into 1's and 0's
  rownames(m) <- f_union
  
  # Working with scores for each top N run
  s <- unlist(scores.l)
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
