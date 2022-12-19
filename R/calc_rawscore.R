
#' Calculate Raw Score using the appropriate method to compute scores based on 
#' skewness of a given binary matrix
#'
#' @param x a numeric matrix ot an expression set of binary features.
#' @param method a character string specifies a scoring method that is 
#' used in the search. There are 4 options: (\code{"ks"} or \code{"wilcox"} or 
#' \code{"revealer"} (conditional mutual information from REVEALER) or 
#' \code{"custom"} (a user customized scoring method)). Default is \code{ks}.
#' @param alternative a character string specifies an alternative hypothesis 
#' testing (\code{"two.sided"} or \code{"greater"} or \code{"less"}). 
#' Default is \code{less} for left-skewed significance testing.
#' @param metric a character string specifies a metric to search 
#' for best features. \code{"pval"} or \code{"stat"} may be used which 
#' corresponding to p-value or score statistic. Default is \code{pval}. 
#' NOTE: \code{Revealer} method only utilized score statistics values 
#' (no p-values).
#' @param weights if method is \code{ks}, specifies a vector of weights 
#' will perform a weighted-KS testing. Default is \code{NULL}.   
#' @param input_score a vector of continuous values of a response of 
#' interest (required). The \code{input_score} object must have names or 
#' labels that match the colnames of the expression matrix.
#' @param seed_names one or more features representing known “causes” 
#' of activation or features associated with a response of interest. 
#' It is used for revealer method only.
#' @param custom_function if method is \code{"custom"}, specifies 
#' the customized function here. Default is \code{NULL}.
#' @param custom_parameters if method is \code{"custom"}, specifies a list of 
#' arguments to be passed to the custom_function(). Default is \code{NULL}.
#' @param warning a logical value indicates whether or not to print the 
#' warning message messages. Default is \code{FALSE}. 
#'
#' @return Score
#' @noRd
calc_rawscore <- function(x,
                          method = c("ks","wilcox","revealer", "custom"),
                          alternative = c("less", "greater", "two.sided"), 
                          metric,
                          weights, 
                          input_score,
                          seed_names=NULL,
                          custom_function,
                          custom_parameters, 
                          warning=FALSE){
  
  
  if (class(x)[1] == "ExpressionSet") x <- exprs(x)
  
  score <- switch(
    method,
    ks = ks_genescore_mat(
      mat = x,
      alternative = alternative, 
      weights = weights
    ),
    wilcox = wilcox_genescore_mat(
      mat = x,
      alternative = alternative,
      ranks = NULL
    ),
    revealer = revealer_genescore_mat(
      mat = x,                                   
      input_score = input_score,      
      seed_names = seed_names,
      target_match = "positive",
      assoc_metric = "IC"
    ),
    custom = custom_genescore_mat(
      mat = x,
      input_score = input_score,
      custom_function = custom_function,
      custom_parameters = custom_parameters
    )
  )
  
  if (warning){
    # Check if the returning result has one or two columns: 
    # score or p_value or both
    if(ncol(score) == 1){
      if(colnames(score) == "score" & metric == "pval"){
        warning("metric = 'pval' is provided but the method ",
                "function only return score values. ",
                "Thus, using 'stat' as metric to search for best features.")
        metric <- "stat"
      }else if(colnames(score) == "p_value" & metric == "stat"){
        warning("metric provided is 'stat' but the method function only ",
                "return p-values. Thus, using 'pval' as metric to ",
                "search for best features.")
        metric <- "pval"
      }
    }
    
  }
  score
  
}