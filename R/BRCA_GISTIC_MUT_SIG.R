#' Genomic Data from the TCGA BRCA
#'
#' Binary matrix of mutation and SCNA features comprises of 16873 features where 1/0 vectors 
#' indicating the presence/absence of the feature in the samples. There
#' are 963 samples in this cohort study.
#'
#' @docType data
#'
#' @usage data(BRCA_GISTIC_MUT_SIG)
#'
#' @format An object of class \code{ExpressionSet} from \code{Biobase} package
#' containing a matrix of 16873 rows (features) and 963 columns (samples) 
#' see (\href{http://bit.ly/2rjy29l}{ExpressionSet documentation})
#'
#' @references Kartha VK, Kern JG, Sebastiani P, Zhang L,
#' Varelas X, Monti S (2017) CaDrA: A computational framework for performing 
#' candidate driver analyses using binary genomic features. 
#' (\href{https://www.biorxiv.org/content/early/2017/11/23/221846}{bioRxiv})
#'
"BRCA_GISTIC_MUT_SIG"