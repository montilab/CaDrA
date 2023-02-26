#'
#' Candidate Drivers Search Plot
#'
#' By utilizing the top N results obtained from \code{candidate_search()},
#' we can find the best meta-feature among the top N searches using
#' \code{topn_best()}. \code{meta_plot()} is then used to produce graphics
#' including a tile plot for the top meta-features that associated with a molecular 
#' phenotype of interest (e.g. input_score), the KS enrichment plot of the meta-features, 
#' and lastly, a density diagram of the distribution of the observed
#' input scores sorted from largest to smallest at the top.
#' @param topn_best_list a list of objects returned from \code{candidate_search()} 
#' corresponding to the search of a set of top N features given by top_N value.
#' The topn_best_list contains a list of returned meta-feature set, its corresponding best score, 
#' and observed input scores.
#' @param input_score_label a label that references to the \code{input_score}
#' variable that was used to compute the top N best features. Default is \code{NULL}.
#'
#' @return a density diagram of observed input scores (sorted from largest to lowest), 
#' a tile plot of the top features within the meta-feature set, and a KS enrichment plot 
#' of the meta-feature set for a given distribution (here, this will correspond 
#' to the logical OR of the features)
#' @examples
#'
#' # Load pre-computed Top-N list generated for sim_FS dataset
#' data(topn_list)
#'
#' # With the results obtained from top-N evaluation,
#' # We can find the combination of features that gives the best score in
#' # top N searches
#' topn_best_meta <- topn_best(topn_list = topn_list)
#'
#' # Now we can plot this set of best meta-feature
#' meta_plot(topn_best_list = topn_best_meta)
#'
#' @export
#' @import SummarizedExperiment ggplot2 reshape2
#' @importFrom grid unit.pmax grid.draw
meta_plot <- function(topn_best_list, input_score_label=NULL){

  # Get best meta-feature set and observed input score from top N best results
  best_feature_set <- topn_best_list[["feature_set"]]
  obs_input_score <- topn_best_list[["input_score"]]
  
  ## Order input score from largest to lowest values
  input_score <- obs_input_score[order(obs_input_score, decreasing = TRUE)]
  
  # Re-order the matrix based on the order of input_score
  feature_set <- best_feature_set[,names(input_score)]
    
  # Plot for continuous metric used to rank samples (Ex: ASSIGN scores)
  # Get the input score label
  var_name <- ifelse(is.null(input_score_label), "input_score",
                     input_score_label)
  
  # Plot y axis label
  y_lab <- var_name
  
  # Plot x axis label
  x_lab <- paste("Samples (n = ", ncol(feature_set),")", sep="")
  
  # Make dataframe of measure and sample information,
  # with rows arranged in specific order of measure used for search
  # Here we assume that input_score provided is ordered in the order
  # that the stepwise search was run for the dataset
  var_d <- data.frame(
    "measure" = input_score,
    "sample" = factor(colnames(feature_set), levels=colnames(feature_set))
  )
  
  # This plot assumes that there are two columns in the data frame
  # called 'sample' and 'measure'
  # The 'sample' variable has to be a factor with ordered levels
  # for displaying in the correct sample order
  
  # Katia: adding ".data" to avoid a warning during check:
  # no visible binding for global variable
  m_plot <- ggplot(data=var_d,
                   aes(x=.data$sample, y=.data$measure, group=1)) +
    geom_area(alpha=0.6, fill="deepskyblue4",
              linetype=1,  color="black") +
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
  
  # Extract the feature binary matrix
  mat <- as.matrix(SummarizedExperiment::assay(feature_set))
  
  # Add on the OR function of all the returned entries
  or <- ifelse(colSums(mat)==0, 0, 1)
  
  # combine mat with the OR function
  mat <- rbind(mat, or)
  
  # Get x and y axis data for ks enrichment plot of 
  # cumulative function of individual features (i.e. the OR function)
  enrichment_dat <- ks_plot_coordinates(
    n_x = length(or),
    y = which(or==1),
    weight = NULL,
    alt = "less"
  )
  
  # Plot for feature_set scores
  enrichment_plot <- ks_plot(df = enrichment_dat)
  
  # Give the last row no row name (this is just for the purpose of the plot)
  rownames(mat)[nrow(mat)] <- ""
  
  # Make the OR function have higher values for a different color (red)
  mat[nrow(mat),] <- 2*(mat[nrow(mat),])
  
  # create the meta plot
  x <- mat
  x_m <- reshape2::melt(x)
  
  # Make factor for genes so it stays in order of dataset (feature feature_set)
  # Here we need to reverse the order of the levels because geom_tile
  # plots ascending from bottom to top
  x_m$Var1 <- factor(as.character(x_m$Var1),
                     levels=rev(unique(as.character(x_m$Var1))))
  
  # Plot for mutation/CNA feature profile (with summary OR)
  feature_plot <- ggplot(data=x_m, aes(x=factor(.data$Var2),
                                       y=.data$Var1, fill=.data$value)) +
    geom_tile(colour=NA)+
    scale_fill_gradient2(
      high="red",
      mid="black",
      low = "white",
      midpoint = 1
    ) +
    theme(
      line=element_blank(),
      axis.text.x=element_blank(),
      axis.text.y=element_text(size=8,colour="black",face="bold"),
      panel.border=element_rect(colour="black",fill=NA)
    ) +
    labs(x="", y="Feature")
  
  # Align the three plots, adjusting the widths to match
  plot_tab <- stacked_gtable_max(ggplotGrob(m_plot),
                                 ggplotGrob(feature_plot +
                                              theme(legend.position="none")),
                                 ggplotGrob(enrichment_plot))
  
  # Set heights for each panel
  plot_tab <- setup_tab_heights(plot_tab, c(2,1.25,3))
  
  # Draw combined plot to canvas
  grid.draw(plot_tab)
  
}

# Function to convert units depending on the version of R
# (for adjusting the dimensions of the ggplots)
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

  stopifnot(
    all(lapply(gtl, is.gtable) |> unlist())
  )

  bind2 <- function(x, y){

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

#' KS Enrichment Plot
#'
#' An enrichment plot of a sum statistic across a given meta-feature using KS method. 
#' The x and y-axis data for the KS enrichment plot can be returned from ks_plot_coordinates()
#' @param df a data frame object with columns 'X' and 'Y' containing the
#' x, y coordinates for the KS enrichment plot
#' @return A graphic of enrichment scores for a given distribution using KS method
#' 
#' @noRd
#'
#' @examples
#'
#' # Load R library
#' library(SummarizedExperiment)
#'
#' # Load pre-computed Top-N list generated for sim_FS dataset
#' data(topn_list)
#'
#' # With the results obtained from top-N evaluation,
#' # We can find the combination of features that gives
#' # the best score in top N searches
#' topn_best_meta <- topn_best(topn_list=topn_list)
#'
#' # Extract the best meta-feature set
#' best_meta_feature <-  topn_best_meta[["feature_set"]]
#'
#' # Get the meta-feature matrix
#' mat <- as.matrix(SummarizedExperiment::assay(best_meta_feature))
#'
#' # Compute logical OR of all returned features in meta-feature matrix
#' or <- ifelse(colSums(mat)==0, 0, 1)
#'
#' # Extract x and y axis data for KS enrichment plot of
#' # cumulative function of individual features (e.g. logical OR function)
#' ks_enrichment_coordinates <- ks_plot_coordinates(
#'    n_x = length(or),
#'    y = which(or==1),
#'    weight = NULL,
#'    alt = "less"
#' )
#'
#' # Plot the enrichment scores
#' ks_plot(df = ks_enrichment_coordinates)
#'
#' @import ggplot2
ks_plot <- function(df){

  # Katia: adding ".data" to avoid a warning during check:
  # no visible binding for global variable
  g <- ggplot(data = df, aes(x=.data$X, y=.data$Y))

  g <- g +
    #geom_line(size=1.25,colour="blueviolet")+
    geom_line(size=1.25, colour="darkgoldenrod1") +
    geom_hline(yintercept=0, linetype=2) +
    # Katia: adding ".data" to avoid a warning during check:
    # no visible binding for global variable
    geom_point(data=df[which.max(df$Y),], aes(x=.data$X, y=.data$Y),
               colour="black", fill="red", size=3, shape=21) +
    annotate("text", x=df[which.max(df$Y),1]+8, y=max(df$Y),
             label=as.character(round(max(df$Y),3)), size=3) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_reverse() + # This inverts the feature_set score statistic line
    theme_classic() +
    theme(panel.border=element_rect(colour="black", fill=NA)) +
    labs(x="", y="Enrichment Score (feature_set)")

  g

}




#' x, y coordinates for creating ks plot
#'
#' Return a data frame with x, y coordinates using ks method
#' which will later used to generate the enrichment plot
#' @param n_x length of ranked list
#' @param y positions of geneset items in ranked list (ranks)
#' @param weight a vector of weights
#' @param alt alternative hypothesis for p-value calculation
#' (\code{"two.sided"} or \code{"greater"} or \code{"less"}).
#' Default is \code{less} for left-skewed significance testing.
#'
#' @useDynLib CaDrA ks_plot_wrap_
#' @noRd
#'
#' @return return a data frame with x, y coordinates to pass to ks_plot()
ks_plot_coordinates <- function(n_x, y, weight, alt=c("less", "greater", "two.sided")){
  
  if(length(alt) > 0){
    alt_int<- switch(alt, two.sided=0L, less=1L, greater=-1L, 1L)
  } else {
    alt_int <- 1L
  }
  y <- as.integer(y)
  n_x <- as.integer(n_x)
  res <- .Call(ks_plot_wrap_, n_x, y, weight, alt_int)
  res <- res[!is.na(res$X), ]
  res
  
}


