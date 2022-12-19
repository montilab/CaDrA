
#' \code{REVEALER} Scoring Method
#' 
#' Compute conditional mutual information of \code{x} and \code{y} 
#' given \code{z} for each row of a given binary feature matrix
#' @param mat a matrix of binary features (required)
#' @param input_score a vector of continuous scores of a response of interest 
#' (required). \code{input_score} must have labels or names that associated 
#' with the colnames of the feature matrix. 
#' @param target_match a direction of target matching (\code{"negative"} or 
#' \code{"positive"}). Use \code{"positive"} to match higher values of 
#' \code{input_score}, \code{"negative"} to match lower values of 
#' \code{input_score}. Default is \code{positive}. 
#' @param seed_names one or more features representing known “causes” 
#' of activation or features associated with a response of interest
#' @param assoc_metric an association metric: \code{"IC"} for information 
#' coefficient or \code{"COR"} for correlation. Default is \code{IC}. 
#'
#' @return a data frame with one column: \code{score}
#' @examples
#' 
#' # Load R library
#' library(Biobase)
#'
#' # Load pre-computed expression set
#' data(sim.ES)
#' 
#' # set seed
#' set.seed(123)
#' 
#' # Provide a vector of continuous scores for a target profile
#' input_score <- rnorm(n = ncol(sim.ES))
#' names(input_score) <- colnames(sim.ES)
#' 
#' # Define additional parameters and start the function
#' result <- revealer_rowscore(
#'   mat = exprs(sim.ES), input_score = input_score, 
#'   target_match = "positive", assoc_metric = "IC"
#' )
#'  
#' @export
revealer_rowscore <- function
(
  mat,                                   
  input_score, 
  target_match = c("positive", "negative"),             
  seed_names = NULL,
  assoc_metric = c("IC", "COR")
)
{
  
  verbose("Using Revealer's Mutually Exclusive method for features scoring")
  
  assoc_metric <- match.arg(assoc_metric)
  target_match <- match.arg(target_match)
  
  ## Make sure mat variable is a matrix
  mat <- as.matrix(mat)
  
  # If mat has only one column, it must be converted to a row-wise matrix 
  # form as it is needed for backward_forward_search() computation
  # mat must have rownames to track features and columns to track samples
  # for n = 1 case, it is only in backward_forward_search(), thus we can 
  # assign a random labels to it
  if(ncol(mat) == 1){
    mat <- matrix(t(mat), nrow=1, byrow=TRUE, 
                  dimnames = list("my_label", rownames(mat))) 
  }
  
  # Check if the matrix has only binary values and no empty values
  if(length(mat) == 0 || !is.matrix(mat) || 
     any(!mat %in% c(0,1)) || any(is.na(mat)))
    stop("mat variable must be a matrix of binary values (no empty values).")
  
  # Check if input_score is provided and no empty values
  if(length(input_score) == 0 || 
     any(!is.numeric(input_score)) || any(is.na(input_score)))
    stop("input_score variable must be provided and are numeric ",
         "with no empty values.\n")
  
  # Make sure the mat variable has rownames for features tracking
  if(is.null(rownames(mat)))
    stop("The mat object does not have rownames or featureData to ",
         "track the features by. Please provide unique features or ",
         "rownames for the expression matrix.\n")
  
  # Make sure the input_score variable has names as the colnames in mat
  if(is.null(names(input_score)))
    stop("The input_score object must have names or labels to track ",
         "the samples by. Please provide the sample names or labels that ",
         "matches the colnames of the expression matrix.\n")
  
  # Make sure the input_score has the same length as number of samples in mat
  if(length(input_score) != ncol(mat)){
    stop("The input_score variable must have the same length as ",
         "the number of columns in mat.\n")
  }else{
    # check if input_score has any labels or names
    if(length(names(input_score)) == 0){
      stop("The input_score object must have names or labels that ",
           "match the colnames of the expression matrix.\n")
    }
    
    # check if input_score has labels or names that matches the 
    # colnames of the expression matrix
    if(any(!names(input_score) %in% colnames(mat))){
      stop("The input_score object have names or labels that do not match ",
           "the colnames of the expression matrix.\n")
    }
    
    # match colnames of expression matrix with names of 
    # provided input_score values
    # if nrow = 1, if it is, convert to matrix form as it is needed for 
    # backward_forward_search with one dimension matrix computation
    if(nrow(mat) == 1){
      mat <- matrix(t(mat[,names(input_score)]), nrow=1, byrow=TRUE, 
                    dimnames = list(rownames(mat), colnames(mat))) 
    }else{
      mat <- mat[,names(input_score)]
    }
  }
  
  # Check if seed_names is provided
  if(length(seed_names) > 0){
    if(any(!seed_names %in% rownames(mat))){
      stop(paste0("The provided seed_names, ", 
                  paste0(seed_names, collapse = ","), ", 
                  does not exist among the rownames of expression matrix."))
    } 
  }
  
  # Check if the dataset has any all 0 or 1 features 
  # (these are to be removed since they are not informative)
  if(any(rowSums(mat) == 0) || any(rowSums(mat) == ncol(mat))){
    warning("The provided matrix has some features that are either all 0s or 1s.",
            "These features will be removed from downsteam computation.\n")
    mat <- mat[!(rowSums(mat) == 0 | rowSums(mat) == ncol(mat)),]
  }
  
  # Make sure matrix is not empty after removing uninformative features
  if(nrow(mat) == 0){
    stop("After removing features that are either all 0s or 1s. ",
         "There are no more features remained for downsteam computation.\n")
  }
  
  # Define seed from given seed_names
  if(length(seed_names) == 0){
    seed <- as.vector(rep(0, ncol(mat)))      
  } else {
    if (length(seed_names) > 1) {
      # Consolidate or summarize one or more seeds to one vector of values
      seed <- as.numeric(ifelse(colSums(mat[seed_names,]) == 0, 0, 1))
    } else {
      seed <- mat[seed_names,]
    }
    locs <- match(seed_names, row.names(mat))
    mat <- mat[-locs,]
  }
  
  # Compute MI and % explained with original seed(s)
  cmi <- seq_len(nrow(mat)) %>% 
    purrr::map_dbl(
      function(r){
        revealer_score(x=input_score, y=mat[r,], 
                       z=seed, assoc_metric=assoc_metric, 
                       target_match=target_match) 
      }
    )
  
  # Convert list to data.frame
  # Only score value from revealer is returned
  cmi <- data.frame(score=cmi)
  rownames(cmi) <- rownames(mat)
  
  return(cmi)
  
}





#' Compute Conditional Mutual Information of x and y given z from \code{REVEALER}
#' 
#' @param x a vector of continuous values of 
#' a given functional response of interest
#' @param y a binary present/absent feature typically representing 
#' genome-wide alterations (mutations, cnvs, amplifications/deletions)
#' @param z a consolidated or summarized vector of values obtained from 
#' one or more binary features(s) which representing known “causes” 
#' of activation or features associated with a given response of interest
#' @param assoc_metric an association metric: information coefficient 
#' (\code{"IC"} by default) or correlation (\code{"COR"}) from \code{REVEALER}
#' @param target_match a direction of target matching (\code{"negative"} or 
#' \code{"positive"}). Use \code{"positive"} to match the higher values of 
#' the target, \code{"negative"} to match the lower values. 
#' Default is \code{positive}.
#'
#' @return a data frame with two columns: \code{score} and \code{p_value}
#' 
#' @examples 
#'
#' # load R library
#' library(purrr)
#' library(Biobase)
#' 
#' # Load pre-computed expression set
#' data(sim.ES)
#' 
#' # set seed
#' set.seed(123)
#' 
#' # Extract the expression matrix
#' mat <- exprs(sim.ES)
#' 
#' # Provide a vector of continuous scores for a target profile
#' input_score <- rnorm(n = ncol(sim.ES))
#' names(input_score) <- colnames(sim.ES)
#' 
#' cmi <- seq_len(nrow(mat)) %>% 
#' purrr::map_dbl(
#'  function(r){
#'    revealer_score(
#'      x=input_score, 
#'      y=mat[r,], 
#'      z=NULL, 
#'      assoc_metric="IC", 
#'      target_match="positive"
#'    ) 
#'  }
#' )
#'   
#' @export
#' @importFrom MASS kde2d bcv
#' @importFrom misc3d kde3d
#' @importFrom stats cor median sd
#' @importFrom ppcor pcor.test
revealer_score <- function
(
  x,
  y, 
  z = NULL,
  assoc_metric = c("IC", "COR"),
  target_match = c("positive", "negative")
)
{
  
  assoc_metric <- match.arg(assoc_metric)
  target_match <- match.arg(target_match)
  
  # Check if x is provided 
  if(length(x) == 0 || any(is.na(x)))
    stop("x must be provided (no empty values).\n")
  
  # Check if y is provided 
  if(length(y) == 0 || any(!y %in% c(0,1)) || any(is.na(y)))
    stop("y must be provided (no empty values).\n")
  
  # Check if the x has the same length as the number of samples in the y
  if(length(x) != length(y))
    stop("The provided x must have the same length ",
         "as the number of samples in y.\n")
  
  # Check if z is provided 
  if(length(z) == 0){
    z <- as.vector(rep(0, length(x)))      
  }else if(length(z) != length(x)){
    stop("The provided z must have the same length as ",
         " the number of samples in x and y",
         " (no empty values).\n")
  }
  
  # Ordering x by target_match direction
  if (target_match == "negative") {
    ind <- order(x, decreasing=FALSE)
  } else {
    ind <- order(x, decreasing=TRUE)
  }
  
  x <- x[ind]
  y <- y[ind]
  z <- z[ind]
  
  # Compute CMI and % explained with or without the provided z
  cmi <- suppressWarnings(
    cond_assoc(x=x, y=y, z=z, metric=assoc_metric)
  )
  
  # Only return score value for revealer
  return(c(score=cmi))
  
}


# Compute the conditional mutual information of x and y given z
cond_assoc <-  function(x, y, z, metric) { 
  
  # Association of x and y given z
  #
  # Conditional mutual information I(x, y | z)
  
  if (length(unique(x)) == 1 || length(unique(y)) == 1) return(0)
  
  if (length(unique(z)) == 1) {  # e.g. for NULLSEED
    if (metric == "IC") {
      return(mutual_inf_v2(x = x, y = y, n.grid = 25)$IC)
    } else if (metric == "COR") {
      return(cor(x, y))
    }
  } else {
    if (metric == "IC") {
      return(cond_mutual_inf(x = x, y = y, z = z, n.grid = 25)$CIC)
    } else if (metric == "COR") {
      return(pcor.test(x, y, z)$estimate)
    }
  }
  
}

# Computes the conditional mutual information x, y | z
cond_mutual_inf <- function(x, y, z, 
                            n.grid=25, 
                            delta = 0.25*c(bcv(x), bcv(y), bcv(z))) {
  
  # Computes the Conditional mutual information: 
  # I(X, Y | X) = H(X, Z) + H(Y, Z) - H(X, Y, Z) - H(Z)
  # The 0.25 in front of the bandwidth is because different conventions 
  # between bcv and kde3d
  
  # Compute correlation-dependent bandwidth
  
  rho <- cor(x, y)
  rho2 <- ifelse(rho < 0, 0, rho)
  delta <- delta*(1 + (-0.75)*rho2)
  
  # Kernel-based prob. density
  
  kde3d.xyz <- kde3d(x=x, y=y, z=z, h=delta, n = n.grid)
  X <- kde3d.xyz$x
  Y <- kde3d.xyz$y
  Z <- kde3d.xyz$z
  PXYZ <- kde3d.xyz$d + .Machine$double.eps
  dx <- X[2] - X[1]
  dy <- Y[2] - Y[1]
  dz <- Z[2] - Z[1]
  
  # Normalize density and calculate marginal densities and entropies
  
  PXYZ <- PXYZ/(sum(PXYZ)*dx*dy*dz)
  PXZ <- colSums(aperm(PXYZ, c(2,1,3)))*dy
  PYZ <- colSums(PXYZ)*dx
  PZ <- rowSums(aperm(PXYZ, c(3,1,2)))*dx*dy
  PXY <- colSums(aperm(PXYZ, c(3,1,2)))*dz
  PX <- rowSums(PXYZ)*dy*dz
  PY <- rowSums(aperm(PXYZ, c(2,1,3)))*dx*dz
  
  HXYZ <- - sum(PXYZ * log(PXYZ))*dx*dy*dz
  HXZ <- - sum(PXZ * log(PXZ))*dx*dz
  HYZ <- - sum(PYZ * log(PYZ))*dy*dz
  HZ <-  - sum(PZ * log(PZ))*dz
  HXY <- - sum(PXY * log(PXY))*dx*dy   
  HX <-  - sum(PX * log(PX))*dx
  HY <-  - sum(PY * log(PY))*dy
  
  MI <- HX + HY - HXY   
  CMI <- HXZ + HYZ - HXYZ - HZ
  
  SMI <- sign(rho) * MI
  SCMI <- sign(rho) * CMI
  
  IC <- sign(rho) * sqrt(1 - exp(- 2 * MI))
  CIC <- sign(rho) * sqrt(1 - exp(- 2 * CMI))
  
  return(list(CMI=CMI, MI=MI, 
              SCMI=SCMI, SMI=SMI, 
              HXY=HXY, HXYZ=HXYZ, 
              IC=IC, CIC=CIC))
  
}

# Computes the Mutual Information/Information Coefficient IC(x, y)
mutual_inf_v2 <- function(x, y, n.grid=25, delta = c(bcv(x), bcv(y))) {
  
  # Computes the Mutual Information/Information Coefficient IC(x, y)
  #
  # Compute correlation-dependent bandwidth
  
  rho <- cor(x, y)
  rho2 <- abs(rho)
  delta <- delta*(1 + (-0.75)*rho2)
  
  # Kernel-based prob. density
  
  kde2d.xy <- kde2d(x, y, n = n.grid, h = delta)
  
  FXY <- kde2d.xy$z + .Machine$double.eps
  dx <- kde2d.xy$x[2] - kde2d.xy$x[1]
  dy <- kde2d.xy$y[2] - kde2d.xy$y[1]
  PXY <- FXY/(sum(FXY)*dx*dy)
  PX <- rowSums(PXY)*dy
  PY <- colSums(PXY)*dx
  HXY <- -sum(PXY * log(PXY))*dx*dy
  HX <- -sum(PX * log(PX))*dx
  HY <- -sum(PY * log(PY))*dy
  
  PX <- matrix(PX, nrow=n.grid, ncol=n.grid)
  PY <- matrix(PY, byrow = TRUE, nrow=n.grid, ncol=n.grid)
  
  MI <- sum(PXY * log(PXY/(PX*PY)))*dx*dy
  rho <- cor(x, y)
  SMI <- sign(rho) * MI
  
  # Use pearson correlation the get the sign (directionality)   
  
  IC <- sign(rho) * sqrt(1 - exp(- 2 * MI)) 
  
  NMI <- sign(rho) * ((HX + HY)/HXY - 1)  
  
  return(list(MI=MI, SMI=SMI, HXY=HXY, HX=HX, HY=HY, NMI=NMI, IC=IC))
  
}


