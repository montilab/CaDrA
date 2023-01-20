#' Simulated Input Scores
#'
#' A simulated vector of continuous scores generated from rnorm(n=ncol(sim_FS),
#' mean=0, sd=1) with set.seed(123) based on the number of samples in the
#' simulated dataset (sim_FS)
#'
#' @docType data
#'
#' @usage data(sim_Scores)
#'
#' @format A vector of continuous scores randomly generated from rnorm(n=ncol(sim_FS),
#' mean=0, sd=1) with set.seed(123) based on the number of samples in the
#' simulated dataset (sim_FS)
#'
#' @references Kartha VK, Kern JG, Sebastiani P, Zhang L,
#' Varelas X, Monti S (2019) CaDrA: A computational framework for performing
#' candidate driver analyses using binary genomic features.
#' (\href{https://www.frontiersin.org/articles/10.3389/fgene.2019.00121/full}{Frontiers in Genetics})
#'
"sim_Scores"
