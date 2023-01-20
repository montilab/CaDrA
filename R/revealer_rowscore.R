#'
#' REVEALER Scoring Method
#'
#' Compute conditional mutual information of \code{x} and \code{y}
#' given \code{z} for each row of a given binary feature matrix
#' @param FS a feature set of binary features. It can be a matrix or
#' a \code{SummarizedExperiment} class object from SummarizedExperiment package.
#' 
#' @param input_score a vector of continuous scores representing a phenotypic
#' readout of interest such as protein expression, pathway activity, etc.
#' The \code{input_score} object must have names or labels that match the column
#' names of FS object.
#' @param seed_names one or more features representing known “causes”
#' of activation or features associated with a response of interest.
#' Default is NULL.
#' @param assoc_metric an association metric: \code{"IC"} for information
#' coefficient or \code{"COR"} for correlation. Default is \code{IC}.
#' @param warning a logical value indicates whether or not to print the 
#' diagnostic messages. Default is \code{TRUE}
#' 
#' @noRd
#'
#' @return a matrix with one column: \code{score}
#' @import SummarizedExperiment
revealer_rowscore <- function
(
  FS,
  input_score,
  seed_names = NULL,
  assoc_metric = c("IC", "COR"),
  warning = TRUE
)
{

  assoc_metric <- match.arg(assoc_metric)
  
  # Check of FS and input_score are valid inputs
  if(warning == TRUE) check_data_input(FS = FS, input_score = input_score, warning=warning)

  # Get the feature names
  feature_names <- rownames(FS)
  
  # Check if seed_names is provided
  if(length(seed_names) == 0){
    seed_vector <- as.vector(rep(0, ncol(FS)))
  }else{
    if(any(!seed_names %in% rownames(FS)))
      stop(paste0("The provided seed_names, ",
                  paste0(seed_names[which(!seed_names %in% rownames(FS))], collapse = ","), ",
                  do not exist among the row names of the FS object"))

    # Consolidate or summarize one or more seeds into one vector of values
    if(length(seed_names) > 1) {
      seed_vector <- as.numeric(ifelse(colSums(assay(FS)[seed_names,]) == 0, 0, 1))
    }else{
      seed_vector <- as.numeric(assay(FS)[seed_names,])
    }
    
    # Remove the seeds from the binary feature matrix
    locs <- match(seed_names, row.names(FS))
    FS <- FS[-locs,]
  }

  # Extract the feature binary matrix
  if(is(FS, "SummarizedExperiment")){
    mat <- as.matrix(SummarizedExperiment::assay(FS))
  }else if(is(FS, "matrix")){
    mat <- as.matrix(FS)
  }else{
    mat <- matrix(t(FS), nrow=1, byrow=TRUE,
                  dimnames=list(feature_names, names(FS)))
  }
  
  # Compute CMI
  cmi <- apply(X=mat, MARGIN=1, function(x){
    revealer_score(
      x = input_score,
      y = x,
      z = seed_vector,
      assoc_metric = assoc_metric
    )
  })

  # Convert results as matrix
  # Revealer method only returns  score statistics (no p-values)
  revealer_mat <- matrix(NA, nrow=nrow(mat), 
                         ncol=1, byrow=TRUE, 
                         dimnames=list(rownames(mat), "score"))
  
  revealer_mat[,1] <- cmi
  
  return(revealer_mat)

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
#' (\code{"IC"} by default) or correlation (\code{"COR"}) from \code{REVEALER}.
#'
#' @noRd
#' 
#' @return a score statistics value
#'
#' @importFrom MASS kde2d bcv
#' @importFrom misc3d kde3d
#' @importFrom stats cor median sd
#' @importFrom ppcor pcor.test
revealer_score <- function
(
  x,
  y,
  z,
  assoc_metric = c("IC", "COR")
)
{

  assoc_metric <- match.arg(assoc_metric)

  # Compute CMI
  cmi <- suppressWarnings(
    cond_assoc(x=x, y=y, z=z, metric=assoc_metric)
  )

  # Revealer only returns score statistics value
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

