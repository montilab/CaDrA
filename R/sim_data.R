#' Simulated Genomic Data
#'
#' A simulated ExpressionSet that comprises of 1000 genomic features (rows) and 200 sample profiles (columns). 
#' Each row feature is represented by a vector of binary number (1/0) indicating the presence/absence of the feature in the samples. 
#' This simulated data includes 10 left-skewed (i.e. True Positive or TP) and 990 uniformly-distributed (i.e. True Null or TN) features.
#'
#' @docType data
#'
#' @usage data(sim.ES)
#'
#' @format An object of class \code{ExpressionSet} from \code{Biobase} package
#' containing a matrix of 1000 rows (features) and 200 columns (samples). 
#' See \href{http://bit.ly/2rjy29l}{ExpressionSet documentation} for more details.
#'
#' @references Kartha VK, Kern JG, Sebastiani P, Zhang L,
#' Varelas X, Monti S (2017) CaDrA: A computational framework for performing 
#' candidate driver analyses using binary genomic features. 
#' (\href{https://www.frontiersin.org/articles/10.3389/fgene.2019.00121/full}{Frontiers in Genetics})
#'
"sim.ES"
