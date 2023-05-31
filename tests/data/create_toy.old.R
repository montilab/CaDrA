N <- 100 # number of observations
SD <- 1  # noise
MU <- 0  # mean of independent variables
         # conditional probability table
         # P(Z=1| X,Y == {00,01,10,11})
CPT <- c("00"=.1,"01"=.5,"10"=.3,"11"=.7) 
PATH <- file.path(".")

## Xc --> Zc <-- Yc
## all continuous
set.seed(123)
Xc <- rnorm(N, mean = MU, sd = SD)
Yc <- rnorm(N, mean = MU, sd = SD)
Zc_XcYc <- rnorm(N, mean = Xc+Yc, sd = SD)

## Xd --> Wc <-- Yc
## mixed discrete/continuous
set.seed(456)
Xd <- sample(c(0,1),size = N, replace = TRUE)
Wc_XdYc <- rnorm(N, mean = Xd + Yc, sd = SD)

## Xd --> Zd <-- Yd
## all discrete

## create conditional probability table P(Zd | Xd, Yd)
##
## generate data
set.seed(789)
Yd <- sample(c(0,1),size = N, replace = TRUE)
probs <- data.frame(Xd, Yd, XY = sprintf("%d%d", Xd, Yd)) |>
  dplyr::mutate(P1 = CPT[XY]) |>
  dplyr::mutate(P0 = 1 - P1)
head(probs)

set.seed(987)
Zd_XdYd <- apply(probs |> dplyr::select(P0,P1),1,function(P) sample(c(0,1),size=1,prob=P))

## Xd --> Zc <-- Yd
## continuous child, discrete parents
set.seed(654)
Zc_XdYd <- rnorm(N,mean = Xd + Yd, sd = SD)

toyDF <- data.frame(
  Xc, Yc, Zc_XcYc, Xd, Wc_XdYc, Yd, Zd_XdYd, Zc_XdYd
)
## marginal correlations
print(as.dist(round(cor(toyDF),3)))
## partial correlations for Xc --> Zc <-- Yc
print(as.dist(ppcor::pcor(toyDF |> dplyr::select(Xc,Yc,Zc_XcYc))$estimate))
## partial correlations for Xd --> Wc <-- Yc
print(as.dist(ppcor::pcor(toyDF |> dplyr::select(Xd,Yc,Wc_XdYc))$estimate))
## partial correlations for Xd --> Zd <-- Yd
print(as.dist(ppcor::pcor(toyDF |> dplyr::select(Xd,Yd,Zd_XdYd))$estimate))
## partial correlations for Xd --> Zc <-- Yd
print(as.dist(ppcor::pcor(toyDF |> dplyr::select(Xd,Yd,Zc_XdYd))$estimate))
## partial correlations of any two variables conditioned on all the others
print(as.dist(ppcor::pcor(toyDF)$estimate))

## save toy dataset
saveRDS(toyDF, file = file.path(PATH,"toy_dataset.rds"))
write.csv(toyDF, file = file.path(PATH,"toy_dataset.csv"))

