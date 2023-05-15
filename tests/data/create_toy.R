## TESTING THE FOLLOWING FUNCTIONS
##
## mutual_inf_cc_1d
## mutual_inf_cc_2d
## mutual_inf_cd_1d
## mutual_inf_cd_2d
## cond_mutual_inf_ccc_1d
## cond_mutual_inf_cdd_1d
## cond_mutual_inf_ccc_2d
## cond_mutual_inf_cdd_2d
##  
setwd(file.path(Sys.getenv("CBMGIT"),"CaDrA","tests","data"))
N <- 100 # number of observations
SD <- 1  # noise
MU <- 0  # mean of independent variables
PATH <- file.path(".")
do_save <- TRUE

## Xc --> Zc <-- Yc (to test ccc_1d)
## all continuous
set.seed(123)
Xc <- rnorm(N, mean = MU, sd = SD)
Yc <- rnorm(N, mean = MU, sd = SD)
Zc_XcYc <- rnorm(N, mean = Xc + Yc, sd = SD)

## Xc --> Zc <-- Yc,Wc (to test ccc_2d)
## all continuous
Wc <- rnorm(N, mean = MU, sd = SD)
Zc_XcYcWc <- rnorm(N, mean = Xc + Yc + Wc, sd = SD)

## Xd --> Zc <-- Yd
## continuous child, discrete parents (to test cdd_1d)
set.seed(456)
Xd <- sample(c(0, 1), size = N, replace = TRUE)
Yd <- sample(c(0, 1), size = N, replace = TRUE)
Zc_XdYd <- rnorm(N, mean = Xd + Yd, sd = SD)

## Xd --> Zc <-- Yd,Wd
## continuous child, discrete parents (to test cdd_2d)
Wd <- sample(c(0, 1), size = N, replace = TRUE)
Zc_XdYdWd <- rnorm(N, mean = Xd + Yd + Wd, sd = SD)

toyDF <- data.frame(
  Xc, Yc, Wc, Zc_XcYc, Zc_XcYcWc, Xd, Yd, Wd, Zc_XdYd, Zc_XdYdWd
)
## marginal correlations
print(as.dist(round(cor(toyDF),3)))
## partial correlations for Xc --> Zc <-- Yc
print(as.dist(ppcor::pcor(toyDF |> dplyr::select(Xc,Yc,Zc_XcYc))$estimate))
## partial correlations for Xc --> Zc <-- Yc,Wc
print(as.dist(ppcor::pcor(toyDF |> dplyr::select(Xc,Yc,Wc,Zc_XcYcWc))$estimate))
## partial correlations for Xd --> Zc <-- Yd
print(as.dist(ppcor::pcor(toyDF |> dplyr::select(Xd,Yd,Zc_XdYd))$estimate))
## partial correlations for Xd --> Zc <-- Yd,Wd
print(as.dist(ppcor::pcor(toyDF |> dplyr::select(Xd,Yd,Wd,Zc_XdYdWd))$estimate))
## partial correlations of any two variables conditioned on all the others
print(as.dist(ppcor::pcor(toyDF)$estimate))

## save toy dataset
if (do_save) {
  saveRDS(toyDF, file = file.path(PATH,"toy_dataset.rds"))
  write.csv(toyDF, file = file.path(PATH,"toy_dataset.csv"))
}
