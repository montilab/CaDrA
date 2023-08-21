#'
#' knnmi Scoring Method
#'
#' Compute conditional mutual information of \code{x} and \code{y}
#' given \code{z} for each row of a given binary feature matrix
#' @param FS a matrix of binary features where rows represent features of 
#' interest (e.g. genes, transcripts, exons, etc...) and columns represent 
#' the samples.
#' @param input_score a vector of continuous scores representing a phenotypic
#' readout of interest such as protein expression, pathway activity, etc.
#' The \code{input_score} object must have names or labels that match the column
#' names of FS object.
#' @param meta_feature a vector of one or more features representing known 
#' causes of activation or features associated with a response of interest, 
#' \code{e.g. input_score}. Default is NULL.
#' @param assoc_metric an association metric: \code{"IC"} for information
#' coefficient or \code{"COR"} for correlation. Default is \code{IC}.
#' 
#' @noRd
#' 
#' @examples 
#'  
#' mat <- matrix(c(1,0,1,0,0,0,0,0,1,0, 
#'                 0,0,1,0,1,0,1,0,0,0,
#'                 0,0,0,0,1,0,1,0,1,0), nrow=3)
#' 
#' colnames(mat) <- 1:10
#' row.names(mat) <- c("TP_1", "TP_2", "TP_3")
#' 
#' set.seed(42)
#' input_score = rnorm(n = ncol(mat))
#' names(input_score) <- colnames(mat)
#' 
#' knnmi_rs <- knnmi_rowscore(
#'   FS = mat, 
#'   input_score = input_score, 
#'   assoc_metric = "IC"
#' )
#' 
#' @return return a vector of row-wise scores where its labels or names 
#' must match the row names of \code{FS} object
#' 
#' @import knnmi
knnmi_rowscore <- function
(
  FS,
  input_score,
  meta_feature = NULL,
  assoc_metric = c("IC", "COR")
)
{

  assoc_metric <- match.arg(assoc_metric)
  
  # Check if meta_feature is provided
  if(is.null(meta_feature)){
    meta_vector <- as.vector(rep(0, ncol(FS)))
    
    cmi <- mutual_inf_cd(input_score, FS)
    
  }else{
    # Getting the position of the known meta features
    locs <- match(meta_feature, row.names(FS))
    # Taking the union across the known meta features
    if(length(locs) > 1) {
      meta_vector <- as.numeric(ifelse(colSums(FS[locs,]) == 0, 0, 1))
    }else{
      meta_vector <- as.numeric(FS[locs,])
    }
    # Remove the seeds from the binary feature matrix
    FS <- FS[-locs, , drop=FALSE]
    
    cmi <- cond_mutual_inf(input_score, FS, meta_vector )
  }
  names(cmi) <- rownames(FS)
  
  return(cmi)

}


