#'
#' REVEALER Scoring Method
#'
#' Compute conditional mutual information of \code{x} and \code{y}
#' given \code{z} for each row of a given binary feature matrix
#' @param FS_mat a matrix of binary features where 
#' rows represent features of interest (e.g. genes, transcripts, exons, etc...)
#' and columns represent the samples.
#' @param input_score a vector of continuous scores representing a phenotypic
#' readout of interest such as protein expression, pathway activity, etc.
#' The \code{input_score} object must have names or labels that match the column
#' names of FS_mat object.
#' @param seed_names a vector of one or more features representing known causes
#' of activation or features associated with a response of interest.
#' Default is NULL.
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
#' revealer_rs <- revealer_rowscore(
#'   FS_mat = mat, 
#'   input_score = input_score, 
#'   assoc_metric = "IC"
#' )
#' 
#' @return return a vector of row-wise scores where its labels or names 
#' must match the row names of \code{FS_mat} object
#' 
revealer_rowscore <- function
(
  FS_mat,
  input_score,
  seed_names = NULL,
  assoc_metric = c("IC", "COR")
)
{

  assoc_metric <- match.arg(assoc_metric)
  
  # Check if seed_names is provided
  if(length(seed_names) == 0){
    seed_vector <- as.vector(rep(0, ncol(FS_mat)))
  }else{
    
    # Check if seed_names exit among the row names of the FS_mat object 
    if(any(!seed_names %in% rownames(FS_mat)))
      stop(paste0(
        "The provided seed_names, ",
        paste0(seed_names[which(!seed_names %in% rownames(FS_mat))], collapse = ","), 
        ", do not exist among the row names of the FS_mat object")
      )
    
    # Consolidate or summarize one or more seeds into one vector of values
    if(length(seed_names) > 1) {
      seed_vector <- as.numeric(ifelse(colSums(FS_mat[seed_names,]) == 0, 0, 1))
    }else{
      seed_vector <- as.numeric(FS_mat[seed_names,])
    }
    
    # Remove the seeds from the binary feature matrix
    locs <- match(seed_names, row.names(FS_mat))
    FS_mat <- FS_mat[-locs,]
    
  }
  
  # Compute CMI
  cmi <- apply(X=FS_mat, MARGIN=1, function(x){
    revealer_score(
      x = input_score,
      y = x,
      z = seed_vector,
      assoc_metric = assoc_metric
    )
  })

  names(cmi) <- rownames(FS_mat)
  
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
                            delta = 0.25*c(MASS::bcv(x), MASS::bcv(y), MASS::bcv(z))) {

  # Computes the Conditional mutual information:
  # I(X, Y | X) = H(X, Z) + H(Y, Z) - H(X, Y, Z) - H(Z)
  # The 0.25 in front of the bandwidth is because different conventions
  # between bcv and kde3d

  # Compute correlation-dependent bandwidth

  rho <- cor(x, y)
  rho2 <- ifelse(rho < 0, 0, rho)
  delta <- delta*(1 + (-0.75)*rho2)

  # Kernel-based prob. density

  kde3d.xyz <- misc3d::kde3d(x=x, y=y, z=z, h=delta, n = n.grid)
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
mutual_inf_v2 <- function(x, y, n.grid=25, delta = c(MASS::bcv(x), MASS::bcv(y))) {

  # Computes the Mutual Information/Information Coefficient IC(x, y)
  #
  # Compute correlation-dependent bandwidth

  rho <- cor(x, y)
  rho2 <- abs(rho)
  delta <- delta*(1 + (-0.75)*rho2)

  # Kernel-based prob. density

  kde2d.xy <- MASS::kde2d(x, y, n = n.grid, h = delta)

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


