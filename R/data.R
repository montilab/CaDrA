#' Toy Dataset for knnmi package
#'
#' @docType data
#'
#' @usage data(mutual_info_df)
#'
#' @format A data frame with 100 rows and 10 columns
#'
"mutual_info_df"
## Some parameters
######################################
# N <- 100 # number of observations
# SD <- 1  # noise
# MU <- 0  # mean of independent variables
# # conditional probability table
# CPT <- c("00"=.1,"01"=.5,"10"=.3,"11"=.7)
# PATH <- file.path(".")
#
# ######################################
# ## Case 1: all continuous
# ## Xc --> Zc <-- Yc
# ######################################
# set.seed(123)
# Xc <- rnorm(N, mean = MU, sd = SD)
# Yc <- rnorm(N, mean = MU, sd = SD)
# Zc_XcYc <- rnorm(N, mean = Xc+Yc, sd = SD)
#
# ######################################
# # Case 2: mixed discrete/continuous
# ## Xd --> Wc <-- Yc
# ######################################
# set.seed(456)
# Xd <- sample(c(0L,1L),size = N, replace = TRUE)
# Wc_XdYc <- rnorm(N, mean = Xd + Yc, sd = SD)
#
# ######################################
# ## Case 3: all discrete
# ## Xd --> Zd <-- Yd
# ######################################
#
# ## create conditional probability table P(Zd | Xd, Yd)
# ##
# ## generate data
# set.seed(789)
# Yd <- sample(c(0L,1L),size = N, replace = TRUE)
# probs <- data.frame(Xd, Yd, XY = sprintf("%d%d", Xd, Yd)) |>
#   dplyr::mutate(P1 = CPT[XY]) |>
#   dplyr::mutate(P0 = 1 - P1)
# head(probs)
#
# set.seed(987)
# Zd_XdYd <- apply(probs |> dplyr::select(P0,P1),1,function(P) sample(c(0L,1L),size=1,prob=P))
#
# ######################################
# ## Create a data frame
# ######################################
# mutual_info_df <- data.frame(
#   Xc,Yc,Zc_XcYc,Xd,Wc_XdYc,Yd,Zd_XdYd
# )
#
# # usethis::use_data(mutual_info_df )
#
# ## save toy dataset
# #saveRDS(mutual_info_df, file = file.path(PATH,"toy_dataset.rds"))
# #write.csv(mutual_info_df, file = file.path(PATH,"toy_dataset.csv"))
#
# ######################################
# ## Check marginal/partial correlations
# ######################################
# ## marginal correlations
# as.dist(round(cor(mutual_info_df),3))
# #             Xc     Yc Zc_XcYc     Xd Wc_XdYc     Yd
# # Yc      -0.050
# # Zc_XcYc  0.478  0.612
# # Xd      -0.043  0.042  -0.061
# # Wc_XdYc  0.003  0.700   0.412  0.365
# # Yd       0.006  0.068   0.049  0.090   0.102
# # Zd_XdYd -0.060 -0.076  -0.139  0.257   0.026  0.410
#
# ## partial correlations for Xc --> Zc <-- Yc
# as.dist(ppcor::pcor(mutual_info_df |> dplyr::select(Xc,Yc,Zc_XcYc))$estimate)
# #                 Xc         Yc
# # Yc      -0.4917019
# # Zc_XcYc  0.6428611  0.7241640
#
# ## partial correlations for Xd --> Wc <-- Yc
# as.dist(ppcor::pcor(mutual_info_df |> dplyr::select(Xd,Yc,Wc_XdYc))$estimate)
# #                 Xd         Yc
# # Yc      -0.3214877
# # Wc_XdYc  0.4706841  0.7359997
#
# ## partial correlations for Xd --> Zd <-- Yd
# as.dist(ppcor::pcor(mutual_info_df |> dplyr::select(Xd,Yd,Zd_XdYd))$estimate)
# #                 Xd         Yd
# # Yd      -0.0178907
# # Zd_XdYd  0.2428992  0.4022989
