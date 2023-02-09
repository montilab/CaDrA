#' Simulated Genomic Data
#'
#' A simulated SummarizedExperiment object that comprises of 1000 genomic
#' features (rows) and 100 sample profiles (columns).
#' Each row is represented by a vector of binary values (1/0)
#' indicating the presence/absence of the feature in the samples.
#' This simulated data includes 10 left-skewed (i.e. True Positive or TP)
#' and 990 uniformly-distributed (i.e. True Null or TN) features.
#'
#' @docType data
#'
#' @usage data(sim_FS)
#'
#' @format An object of class \code{SummarizedExperiment} from
#' \code{SummarizedExperiment} package containing an assay of 1000 rows (features)
#' and 100 columns (samples). Each row is represented by a vector of binary values (1/0)
#' indicating the presence/absence of the feature in the samples.
#' 
#' See \code{?SummarizedExperiment} for more details.
#'
#' @references Kartha VK, Kern JG, Sebastiani P, Zhang L,
#' Varelas X, Monti S (2019) CaDrA: A computational framework for performing
#' candidate driver analyses using binary genomic features.
#' (\href{https://www.frontiersin.org/articles/10.3389/fgene.2019.00121/full}{Frontiers in Genetics})
#'
"sim_FS"
