
#' Candidate heuristic search plot
#' 
#' Plot the candidate search results for a given CaDrA run. Plot will include an optional bar plot of the continuous ranking variable (top), 
#' @param ESet an ESet object containing only the features returned from the candidate_search() function (if any) 
#' @param var_score an optional integer vector of continuous measures used to rank samples for the stepwise search (assumed in matching order)
#' @param var_name a string object describing the name of the continuous measure used, which will be used as the y-axis label for the metric plot
#' @return A plot graphic with the ranked metric plot (optional), a tile plot of the features within the provided ESet, and the corresponding Enrichment Score (ES) for a given distribution (here, this will correspond to the logical OR of the features)
#' @examples
#' data(sim.ES)
#' data(topn.list)
#' 
#' # Plot the results from a top-N evaluation by passing the resulting ESet from a specific run
#' # To find the combination of features that had the best score
#' best_meta <- topn_best(topn_list=topn.list) 
#' 
#' # Now we can plot this set of features
#' meta_plot(best_meta$ESet)
#' 
#' # If a continuous ranking variable was used for the sample-ranking, we can visualize it together
#' # Just for illustration purposes, we simulate random sample scores
#' sample_scores <- sort(runif(ncol(sim.ES),0,1), decreasing=TRUE)
#' meta_plot(ESet=best_meta$ESet, var_score=sample_scores, var_name="My ranking variable")
#' @export
#' @import ggplot2 reshape2
#' @importFrom grid unit.pmax grid.draw
meta_plot <- function(
  ESet,                       #ExpressoinSet containing somatic mutation/CNA data with samples ordered by a given measure
  var_score = NULL,           #Ordered vector of measure used for search
  var_name = ""
){
  
  # Plot y axis label
  y_lab <- var_name
  
  # Plot x axis label
  x_lab <- paste("Samples (n = ", ncol(ESet),")", sep="")
  
  # Plot for continuous metric used to rank samples (Ex: ASSIGN scores)
  if(!is.null(var_score)){
    
    # Make dataframe of metric and sample information, with rows arranged in specific order of measure used for search
    # Here we assume that var_score provided is ordered in the order that the stepwise search was run for the dataset
    var_d <- data.frame("measure"=var_score, "sample"=factor(colnames(ESet), levels=colnames(ESet)))
    
    # This plot assumes that there are two columns in the data frame called 'sample' and 'measure'
    # The 'sample' variable has to be a factor with ordered levels for displaying in the correct sample order
    
    # Katia: adding ".data" to avoid a warning during check:
    # no visible binding for global variable 
    m_plot <- ggplot(data=var_d, aes(x=.data$sample, y=.data$measure, group=1)) +
      geom_area(alpha=0.6, fill="deepskyblue4", linetype=1, size=0.5, color="black") +
      scale_y_continuous(expand = c(0,0)) +
      scale_x_discrete(expand = c(0,0)) +
      theme_classic() +
      labs(x=x_lab, y=y_lab) +
      theme(
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x=element_line(color="black"),
        axis.line.y=element_line(color="black")
      )
    
  }
  
  # Start working with the feature set for plots
  # We are going to fetch the logical OR summary for the set of feature
  # Along with the corresponding ES running score statistic for the OR summary
  mat <- exprs(ESet)
  
  # Add on the OR function of all the returned entries
  or <- ifelse(colSums(mat)==0, 0, 1)
  mat <- rbind(mat, or)
  
  # Get x and y axis data for ES plot of cumulative function of individual features (i.e. the OR function)
  ES_dat <- ks_gene_score(n.x=length(or), y=which(or==1), plot_dat = TRUE, alternative = "less")
  
  # Plot for ES scores
  ES_plot <- plot_ESet(d = ES_dat)
  
  # Give the last row no row name (this is just for the purpose of the plot)
  rownames(mat)[nrow(mat)] <- ""
  
  mat[nrow(mat),]<-2*(mat[nrow(mat),]) #Make the OR function have higher values for a different color (red)
  
  m <- ExpressionSet(assayData = mat, featureData = AnnotatedDataFrame(data.frame("Name"=rownames(mat), row.names = rownames(mat))))
  
  x <- exprs(m)
  
  x_m <- melt(x)
  
  # Make factor for genes so it stays in order of dataset (feature ESet)
  # Here we need to reverse the order of the levels because geom_tile plots ascending from bottom to top
  x_m$Var1 <- factor(as.character(x_m$Var1), levels=rev(unique(as.character(x_m$Var1))))
  
  # Plot for mutation/CNA feature profile (with summary OR)
  # Katia: adding ".data" to avoid a warning during check:
  # no visible binding for global variable 
  feature_plot <- ggplot(data=x_m, aes(x=factor(.data$Var2), y=.data$Var1, fill=.data$value)) +
    geom_tile(colour=NA)+ #So that there are no line borders around each cell
    scale_fill_gradient2( #This is only if we want different colors for different values
      high="red", #This will be for our OR function values (which have 2's for 1's, based on our pre-processing) 
      mid="black", #This is for our individual feature 1 values
      low = "white", #This is for the background (0) values
      midpoint = 1
    ) + 
    theme(
      line=element_blank(),
      axis.text.x=element_blank(),
      axis.text.y=element_text(size=8,colour="black",face="bold"),
      panel.border=element_rect(colour="black",fill=NA)
    ) +
    labs(x="", y="Feature") 
  
  if(!is.null(var_score)){
    # Align the three plots, adjusting the widths to match
    plot_tab <- stacked_gtable_max(ggplotGrob(m_plot), ggplotGrob(feature_plot + theme(legend.position="none")), ggplotGrob(ES_plot))
    
    # Set heights for each panel
    plot_tab <- setup_tab_heights(plot_tab, c(2,1.25,3))
    
  } else{
    # Align the three plots, adjusting the widths to match
    plot_tab <- stacked_gtable_max(ggplotGrob(feature_plot + theme(legend.position="none")), ggplotGrob(ES_plot))
    
    # Set heights for each panel
    plot_tab <- setup_tab_heights(plot_tab, c(1.25, 3))
  }
  
  
  # Draw combined plot to canvas
  grid.draw(plot_tab)
  
}

# Function to convert units depending on the version of R (for adjusting the dimensions of the ggplots)
setup_tab_heights <- function(plot_tab, dims) {
  
  panels <- plot_tab$layout$t[grep("panel", plot_tab$layout$name)]
  has_unit_vec <- getRversion() > '3.3'
  
  if(has_unit_vec) {
    plot_tab$heights[panels] <- unit(dims, units="null")
  }else {
    plot_tab$heights[panels] <- lapply(dims, unit, "null") 
  }
  
  return(plot_tab)
  
}


# Solution to drawing arbitrary number of plots stacked above one another
# Each plot will be left-aligned
# Each plot should first be made a ggplotGrob object
# Ex: grid.draw(stacked_gtable_max(ggplotGrob(ggplot1), ggpotGrob(ggplot2)))
#' @import gtable
#' @importFrom grid unit.pmax
stacked_gtable_max <- function(...){
  
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
plot_ESet <- function(
  d 
){
  
  # Katia: adding ".data" to avoid a warning during check:
  # no visible binding for global variable 
  g <- ggplot(data = d, aes(x=.data$x, y=.data$y))
  g <- g +
    #geom_line(size=1.25,colour="blueviolet")+
    geom_line(size=1.25, colour="darkgoldenrod1") +
    geom_hline(yintercept=0, linetype=2) +
    # Katia: adding ".data" to avoid a warning during check:
    # no visible binding for global variable
    geom_point(data=d[which.max(d$y),], aes(x=.data$x, y=.data$y), colour="black", fill="red", size=3, shape=21) +
    annotate("text", x=d[which.max(d$y),1]+8, y=max(d$y),label=as.character(round(max(d$y),3)), size=3) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_reverse() + # This inverts the ES score statistic line
    theme_classic() +
    theme(panel.border=element_rect(colour="black", fill=NA)) +
    labs(x="", y="Enrichment Score (ES)")
  
  g
  
}




