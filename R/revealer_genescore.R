
#' Compute Conditional Mutual Information of x and y given z from \code{REVEALER}
#' 
#' @param x a vector of continuous values of 
#' a given functional response of interest
#' @param y a binary feature for a given response of interest
#' @param z a binary feature often known as causes of activation that 
#' associated with a given response of interest
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
#'    revealer_genescore(
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
revealer_genescore <- function
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
  
  # Convert x, y, z as vector and numeric  
  x <- as.numeric(x); y <- as.numeric(y); z <- as.numeric(z); 
  
  # Check if x is provided 
  if(length(x) == 0 || any(is.na(x)))
    stop("x must be provided and are numeric (no empty values).\n")
  
  # Check if y is provided 
  if(length(y) == 0 || any(!y %in% c(0,1)) || any(is.na(y)))
    stop("y must be provided and are binary (no empty values).\n")
  
  # Check if the x has the same length as the number of samples in the y
  if(length(x) != length(y))
    stop("The provided x variable must have the same length ",
         "as the number of samples in y.\n")
  
  # Check if z is provided 
  if(length(z) > 0){
    if(length(z) != length(x) || any(!z %in% c(0,1)) || any(is.na(z)))
      stop("The provided z variable must have the same length as ",
           " the number of samples in x and y",
           " and contains only binary values (no empty values).\n")
  }else{
    z <- as.vector(rep(0, length(x)))      
  }
  
  # If target_match variable is not specified, use "positive" as default.
  if(length(target_match) == 0 || nchar(target_match) == 0){
    target_match <- "positive"
    warning("The target_match variable is not specified.",
            " Using 'positive' by default...\n")
  }else if(length(target_match) == 1 && 
           !target_match %in% c("positive", "negative")){
    stop(paste0(target_match, collapse=", "), 
         " is not a valid target_match value.",
         " The target_match variable must be 'positive' or 'negative'.\n")
  }else if(length(target_match) > 1 && 
           all(!target_match %in% c("positive", "negative"))){
    stop(paste0(target_match, collapse=", "), 
         " is not a valid target_match value.",
         " The target_match variable must be 'positive' or 'negative'.\n")
  }else if(length(target_match) > 1 && 
           any(target_match %in% c("positive", "negative"))){
    target_match <- target_match[which(target_match %in% 
                                         c("positive", "negative"))][1]
    warning("More than one target_match values were specified.", 
            "Only the first valid target_match value, '", 
            target_match, "', is used.\n")
  }
  
  # If assoc_metric variable is not specified, use "IC" as default.
  if(length(assoc_metric) == 0 || nchar(assoc_metric) == 0){
    assoc_metric <- "IC"
    warning("The assoc_metric variable is not specified.", 
                  "Using 'IC' by default ..\n")
  }else if(length(assoc_metric) == 1 && !assoc_metric %in% c("IC", "COR")){
    stop(paste0(assoc_metric, collapse=", "), 
         " is not a valid assoc_metric value.",
         " The assoc_metric variable must be 'IC' or 'COR'.\n")
  }else if(length(assoc_metric) > 1 && all(!assoc_metric %in% c("IC", "COR"))){
    stop(paste0(assoc_metric, collapse=", "), 
         " is not a valid assoc_metric value.",
         " The assoc_metric variable must be 'IC' or 'COR'.\n")
  }else if(length(assoc_metric) > 1 && any(assoc_metric %in% c("IC", "COR"))){
    assoc_metric <- assoc_metric[which(assoc_metric %in% c("IC", "COR"))][1]
    warning("More than one assoc_metric values were specified.",
            " Only the first valid assoc_metric value, '", 
            assoc_metric, "', is used.\n")
  }
  
  # reordering x by target_match direction
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

# Computes the Conditional mutual information
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





