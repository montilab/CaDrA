
#' Compute conditional mutual information scores of x and y given z
#' 
#' @param x A continuous functional response of interest
#' @param y A binary feature for response of interest
#' @param z A binary feature which is often known as "causes" of activation
#' @param assoc_metric Association Metric: "IC" information coefficient (default) or "COR" correlation
#' @param target_match Direction of the match (negative or positive). Use "positive" to match the higher values of the target, "negative" to match the lower values. Default is positive. 
#' @export
revealer_genescore <- function(
  x,
  y, 
  z = NULL,
  assoc_metric = "IC",
  target_match = "positive"
)
{
  
  # #Required R packages
  # library(NMF)
  # library(MASS)
  # library(ppcor)
  # library(misc3d)
  # library(smacof)
  # library(tidyverse)
  # library(maptools)
  # library(CaDrA)
  # library(Biobase)
  # library(tidyverse)
  # 
  # ## Read in the simulated dataset from CaDrA
  # data("sim.ES")
  # mat <- exprs(sim.ES)
  #
  # ## Define the parameters for testing
  # x = rnorm(n=ncol(mat), mean=0, sd=2)
  # y = mat[1,]
  # z = NULL
  # assoc_metric = "IC"
  # target_match = "positive"
   
#   # Convert x, y, z as vector and numeric input
#   x=as.vector(as.numeric(x)); y=as.vector(as.numeric(y)); z=as.vector(as.numeric(z));
# 
#   # Check if x is provided
#   if(length(x) == 0)
#     stop("x must be provided and are numeric values.\n")
# 
#   # Check if y is provided
#   if(length(y) == 0 || any(!y %in% c(0,1)))
#     stop("y must be provided and contains only binary values.\n")
# 
#   # Check if the x has the same length as the number of samples in the y
#   if(length(x) < length(y))
#     stop("'The provided x variable must have the same length as the number of samples in y.\n")
# 
#   # Check if z is provided
#   if(length(z) > 0){
#     if(length(z) < length(x) || any(!z %in% c(0,1)))
#       stop("'The provided z variable must have the same length as the number of samples in x and y and contains only binary values.\n")
#   }else{
#     z <- as.vector(rep(0, length(x)))
#   }
# 
#   # exclude samples with target == NA
#   x_locs <- seq(1, length(x))[!is.na(x)]
#   y_locs <- seq(1, length(y))[!is.na(y)]
#   z_locs <- seq(1, length(z))[!is.na(z)]
# 
#   # get overlap btw x, y and z
#   overlap <- intersect(x_locs, y_locs) %>% intersect(z_locs)
#   x <- x[overlap]
#   y <- y[overlap]
#   z <- z[overlap]
# 
#   # If target_match variable is not specified, use "positive" as default.
#   if(length(target_match) == 0 || nchar(target_match) == 0){
#     warning("The target_match variable is not specified. Using 'positive' by default ..\n")
#     target_match <- "positive"
#   }else if(length(target_match) == 1 && !target_match %in% c("positive", "negative")){
#     warning(paste0(target_match, collapse=", "), " is not a valid target_match value. The target_match variable must be 'positive' or 'negative'. Using 'positive' by default.\n")
#     target_match <- "positive"
#   }else if(length(target_match) > 1 && all(!target_match %in% c("positive", "negative"))){
#     warning(paste0(target_match, collapse=", "), " is not a valid target_match value. The target_match variable must be 'positive' or 'negative'. Using 'positive' by default.\n")
#     target_match <- "positive"
#   }else if(length(target_match) > 1 && any(target_match %in% c("positive", "negative"))){
#     target_match <- target_match[which(target_match %in% c("positive", "negative"))][1]
#     warning("More than one target_match values were specified. Only the first valid target_match value, '", target_match, "', is used.\n")
#   }
# 
#   # If assoc_metric variable is not specified, use "IC" as default.
#   if(length(assoc_metric) == 0 || nchar(assoc_metric) == 0){
#     warning("The assoc_metric variable is not specified. Using 'IC' by default ..\n")
#     assoc_metric <- "IC"
#   }else if(length(assoc_metric) == 1 && !assoc_metric %in% c("IC", "COR")){
#     warning(paste0(assoc_metric, collapse=", "), " is not a valid assoc_metric value. The assoc_metric variable must be 'IC' or 'COR'. Using 'IC' by default.\n")
#     assoc_metric <- "IC"
#   }else if(length(assoc_metric) > 1 && all(!assoc_metric %in% c("IC", "COR"))){
#     warning(paste0(assoc_metric, collapse=", "), " is not a valid assoc_metric value. The assoc_metric variable must be 'IC' or 'COR'. Using 'IC' by default.\n")
#     assoc_metric <- "IC"
#   }else if(length(assoc_metric) > 1 && any(assoc_metric %in% c("IC", "COR"))){
#     assoc_metric <- assoc_metric[which(assoc_metric %in% c("IC", "COR"))][1]
#     warning("More than one assoc_metric values were specified. Only the first valid assoc_metric value, '", assoc_metric, "', is used.\n")
#   }
# 
#   # reordering x by target_match direction
#   if (target_match == "negative") {
#     ind <- order(x, decreasing=F)
#   } else {
#     ind <- order(x, decreasing=T)
#   }
# 
#   x <- x[ind]
#   y <- y[ind]
#   z <- z[ind]
# 
#   # Compute CMI and % explained with/without the provided seed
#   median_target <- median(x)
# 
#   if (target_match == "negative") {
#     target_locs <- seq(1, length(x))[x <= median_target]
#   } else {
#     target_locs <- seq(1, length(x))[x > median_target]
#   }
# 
#   cmi <- cond_assoc(x=x, y=y, z=z, metric=assoc_metric)
#   pct_explained <- sum(y[target_locs])/length(target_locs)
# 
#   return(data.frame(score=cmi, p_value=pct_explained))
# 
# }
# 
# # Compute the Conditional mutual information
# cond_assoc <-  function(x, y, z, metric) {
# 
#   # Association of x and y given z
#   #
#   # Conditional mutual information I(x, y | z)
# 
#   if (length(unique(x)) == 1 || length(unique(y)) == 1) return(0)
# 
#   if (length(unique(z)) == 1) {  # e.g. for NULLSEED
#     if (metric == "IC") {
#       return(mutual_inf_v2(x = x, y = y, n.grid = 25)$IC)
#     } else if (metric == "COR") {
#       return(cor(x, y))
#     }
#   } else {
#     if (metric == "IC") {
#       return(cond_mutual_inf(x = x, y = y, z = z, n.grid = 25)$CIC)
#     } else if (metric == "COR") {
#       return(pcor.test(x, y, z)$estimate)
#     }
#   }
# 
# }
# 
# mutual_inf_v2 <- function(x, y, n.grid=25, delta = c(bcv(x), bcv(y))) {
# 
#   # Computes the Mutual Information/Information Coefficient IC(x, y)
#   #
#   # Compute correlation-dependent bandwidth
# 
#   rho <- cor(x, y)
#   rho2 <- abs(rho)
#   delta <- delta*(1 + (-0.75)*rho2)
# 
#   # Kernel-based prob. density
# 
#   kde2d.xy <- kde2d(x, y, n = n.grid, h = delta)
#   FXY <- kde2d.xy$z + .Machine$double.eps
#   dx <- kde2d.xy$x[2] - kde2d.xy$x[1]
#   dy <- kde2d.xy$y[2] - kde2d.xy$y[1]
#   PXY <- FXY/(sum(FXY)*dx*dy)
#   PX <- rowSums(PXY)*dy
#   PY <- colSums(PXY)*dx
#   HXY <- -sum(PXY * log(PXY))*dx*dy
#   HX <- -sum(PX * log(PX))*dx
#   HY <- -sum(PY * log(PY))*dy
# 
#   PX <- matrix(PX, nrow=n.grid, ncol=n.grid)
#   PY <- matrix(PY, byrow = TRUE, nrow=n.grid, ncol=n.grid)
# 
#   MI <- sum(PXY * log(PXY/(PX*PY)))*dx*dy
#   rho <- cor(x, y)
#   SMI <- sign(rho) * MI
# 
#   # Use pearson correlation the get the sign (directionality)
# 
#   IC <- sign(rho) * sqrt(1 - exp(- 2 * MI))
# 
#   NMI <- sign(rho) * ((HX + HY)/HXY - 1)
# 
#   return(list(MI=MI, SMI=SMI, HXY=HXY, HX=HX, HY=HY, NMI=NMI, IC=IC))
# 
# }
# 
# 
# cond_mutual_inf <- function(x, y, z, n.grid=25, delta = 0.25*c(bcv(x), bcv(y), bcv(z))) {
# 
#   # Computes the Conditional mutual information:
#   # I(X, Y | X) = H(X, Z) + H(Y, Z) - H(X, Y, Z) - H(Z)
#   # The 0.25 in front of the bandwidth is because different conventions between bcv and kde3d
# 
#   # Compute correlation-dependent bandwidth
# 
#   rho <- cor(x, y)
#   rho2 <- ifelse(rho < 0, 0, rho)
#   delta <- delta*(1 + (-0.75)*rho2)
# 
#   # Kernel-based prob. density
# 
#   kde3d.xyz <- kde3d(x=x, y=y, z=z, h=delta, n = n.grid)
#   X <- kde3d.xyz$x
#   Y <- kde3d.xyz$y
#   Z <- kde3d.xyz$z
#   PXYZ <- kde3d.xyz$d + .Machine$double.eps
#   dx <- X[2] - X[1]
#   dy <- Y[2] - Y[1]
#   dz <- Z[2] - Z[1]
# 
#   # Normalize density and calculate marginal densities and entropies
# 
#   PXYZ <- PXYZ/(sum(PXYZ)*dx*dy*dz)
#   PXZ <- colSums(aperm(PXYZ, c(2,1,3)))*dy
#   PYZ <- colSums(PXYZ)*dx
#   PZ <- rowSums(aperm(PXYZ, c(3,1,2)))*dx*dy
#   PXY <- colSums(aperm(PXYZ, c(3,1,2)))*dz
#   PX <- rowSums(PXYZ)*dy*dz
#   PY <- rowSums(aperm(PXYZ, c(2,1,3)))*dx*dz
# 
#   HXYZ <- - sum(PXYZ * log(PXYZ))*dx*dy*dz
#   HXZ <- - sum(PXZ * log(PXZ))*dx*dz
#   HYZ <- - sum(PYZ * log(PYZ))*dy*dz
#   HZ <-  - sum(PZ * log(PZ))*dz
#   HXY <- - sum(PXY * log(PXY))*dx*dy
#   HX <-  - sum(PX * log(PX))*dx
#   HY <-  - sum(PY * log(PY))*dy
# 
#   MI <- HX + HY - HXY
#   CMI <- HXZ + HYZ - HXYZ - HZ
# 
#   SMI <- sign(rho) * MI
#   SCMI <- sign(rho) * CMI
# 
#   IC <- sign(rho) * sqrt(1 - exp(- 2 * MI))
#   CIC <- sign(rho) * sqrt(1 - exp(- 2 * CMI))
#   
#   return(list(CMI=CMI, MI=MI, SCMI=SCMI, SMI=SMI, HXY=HXY, HXYZ=HXYZ, IC=IC, CIC=CIC))

}


