
#Solution to drawing arbitrary number of plots stacked above one another
#Each plot will be left-aligned
#Each plot should first be made a ggplotGrob object
#Ex: grid.draw(rbind_gtable_max(ggplotGrob(ggplot1),ggpotGrob(ggplot2)))

#' Adjusted row-bind of ggplotGrob objects
#' 
#' Row-bind ggplotGrob objects so that multiple, stacked plots are aligned based on the maximum width of any plots. in the table. Each plot will be left-aligned. Each plot should be a ggplotGrob object. Ex: grid.draw(rbind_gtable_max(ggplotGrob(ggplot1),ggpotGrob(ggplot2)))
#' @param ... multiple ggplotGrob objects using which a list is made
#' @return a gtable object with the width-adjusted stack of ggplots
#' @export
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
    x$heights <- gtable:::insert.unit(x$heights, y$heights)
    x$rownames <- c(x$rownames, y$rownames)
    x$widths <- grid::unit.pmax(x$widths, y$widths)
    x$grobs <- append(x$grobs, y$grobs)
    x
  }
  
  Reduce(bind2, gtl)
}

#' Standardize matrix by row
#' 
#' Standardize a matrix by row, converting the values to a Z score
#' @param mat a matrix object of continuous values to be row-standardized
#' @return A Z-score transformed matrix with the same dimensions as the input matrix
#' @export
row_standard<-function(mat){
  sds<-apply(mat,1,sd)
  means<-apply(mat,1,mean)
  
  mat.std<-(mat-means)/sds
  mat.std
}


#' Enrichment Score (ES) plot
#' 
#' Plot the ES running sum statistic, including the max score (KS) using x and y-axis data returned by ks.genescore() when plot.dat is set to TRUE for a given dataset
#' @param d a data frame object with columns 'x' and 'y' containing the necessary data for an ES plot
#' @return A plot graphic of the Enrichment Score (ES) for a given distribution
#' @export
plot_ES<-function(d 
){
  g<-ggplot(data = d,aes(x=x,y=y))
  g<-g+
    #geom_line(size=1.25,colour="blueviolet")+
    geom_line(size=1.25,colour="darkgoldenrod1")+
    geom_hline(yintercept=0,linetype=2)+
    geom_point(data=d[which.max(d$y),],aes(x=x,y=y),colour="black",fill="red",size=3,shape=21)+
    annotate("text",x=d[which.max(d$y),1]+8,y=max(d$y),label=as.character(round(max(d$y),3)),size=2.5)+
    scale_x_continuous(expand = c(0,0))+
    theme_classic()+
    theme(panel.border=element_rect(colour="black",fill=NA))+
    labs(x="",y="Enrichment Score (ES)")
  
  g
}

#' Step-wise heuristic search plot
#' 
#' Plot together the continuous variable used to rank samples in ks.stepwise()
#' @param ESet an ESet object containing only the features returned from the ks.stepwise() function (if any) 
#' @param d a data frame object with columns 'sample' and 'measure' containing the continuous measures used to rank samples for the stepwise search (in matching order)
#' @param ES.dat a data frame object with columns 'x' and 'y' containing the necessary data for an ES plot returned by the ks.genescore() function when plot.dat is set to TRUE for the given dataset 
#' @param variable_name a string object describing the name of the continuous measure used, which will be used as the y-axis label for the metric plot
#' @param row_std logical indicating whether or not to perform row Z-score transformation on the data before visualizing (only for continuous data). Default is FALSE
#' @return A plot graphic with the ranked metric plot, a tile plot of the features within the provided ESet, and the corresponding Enrichment Score (ES) for a given distribution (here, this will correspond to the union of the features, which we assume is the last row in the ESet feature space)
#' @export
plot_comb<-function(ESet, #ExpressoinSet containing somatic mutation/CNA data with samples ordered by a given measure
                    d, #Dataframe of metric and sample information, in decreasing order of measure
                    ES.dat, #Dataframe with x and y axis information for KS plot
                    variable_name="",
                    row_std=F #Whether or not to perform row Z-score transformation on the data before visualizing
                    
){
  
  #Plot y axis label
  y_lab<-paste(variable_name,"\nscore",sep="")
  #Plot x axis label
  x_lab<-paste("Samples (n = ",ncol(ESet),")",sep="")
  
  
  
  #Plot for metric (Ex: Assign scores)
  m.plot<-ggplot(d,aes(x=sample,y=measure,group=1))+
    geom_area(alpha=0.6,fill="deepskyblue4",linetype=1,size=0.5,color="black")+
    #geom_line(aes(x=sample,y=measure),colour="black")+
    scale_y_continuous(limits=c(0,1),expand = c(0,0))+
    theme_classic()+labs(x=x_lab,y=y_lab)+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
  
  
  #Plot for ES scores
  ks.plot<-plot_ES(d = ES.dat)
  
  x<-exprs(ESet)
  
  
  #Row standardize data (NOTE: ONLY DO THIS IF THE DATA IS CONTINUOUS)
  if(row_std==T)
    x<-row_standard(x)
  x.m<-melt(x)
  
  #Make factor for genes so it stays in order of dataset
  #Here we need to reverse the order of the levels because geom_tile plots ascending from bottome to top
  x.m$Var1<-factor(as.character(x.m$Var1),levels=rev(unique(as.character(x.m$Var1))))
  
  #Plot for mutation/CNA profile
  mut.plot<-ggplot(x.m,aes(x=factor(Var2),y=Var1,fill=value))+
    #geom_tile(colour="black")+
    geom_tile(colour=NA)+
    scale_fill_gradient2( #This is only if we want different colors for different values
      high="red", #This will be for our OR function values (which have 2's for 1's, based on our pre-processing) 
      mid="black", #This is for our individual feature 1 values
      low = "white", #This is for the background (0) values
      midpoint = 1)+ 
    theme(line=element_blank(),
          #      axis.text.x=element_text(size=10,angle=90),
          #      axis.text.y=element_text(size=5),
          axis.text.x=element_blank(),
          axis.text.y=element_text(size=5,colour="black"),
          #  axis.title.y=element_text(face="italic"),
          panel.border=element_rect(colour="black",fill=NA))+
    #labs(x="",y=paste(nrow(x)," mutations",sep="")) 
    labs(x="",y="Feature") 
  #theme(line=element_blank(),axis.text.y=element_blank(),axis.text.x=element_text(angle=90),panel.border=element_rect(colour="black",fill=NA),legend.position="none")
  
  # Extract the legend from cna plot
  legend = gtable_filter(ggplot_gtable(ggplot_build(mut.plot)), "guide-box")
  #grid.draw(legend)
  
  plot.tab<-rbind_gtable_max(ggplotGrob(m.plot),ggplotGrob(mut.plot+theme(legend.position="none")),ggplotGrob(ks.plot))
  
  #Set heights for each panel
  panels <- plot.tab$layout$t[grep("panel", plot.tab$layout$name)]
  plot.tab$heights[panels] <- lapply(c(2,1,3), unit, "null")
  
  grid.draw(plot.tab)
  #grid.arrange(m.plot,cna.plot+theme(legend.position="none"),legend = legend,nrow=2,ncol=1)
  
}
