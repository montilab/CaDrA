## generating a binary feature, and a set of 100 scores with increasing association with the feature.

setwd(file.path(Sys.getenv("CBMGIT"),"CaDrA","tests","data"))
PATH <- file.path(".")
do_save <- TRUE

N <- 100 # number of observations
M <- 50  # number of scores
SD <- 0.25

set.seed(123)
feature <- sample(seq(0,1),size = N, replace = TRUE)
scores <- sapply(
  seq(0, 5, length.out = M),
  function(delta) rnorm(n = N, mean = feature, sd = SD + delta)
)
toy <- data.frame(
  feature = feature,
  scores[,order(cor(feature,scores), decreasing = TRUE)]
)
## show correlation between feature and each of the 50 scores
COR <- cor(toy$feature, toy |> dplyr::select(-feature)) |> drop()
stopifnot(all(diff(COR)<0)) # check that correlations are in descending order
plot(COR, ylab = "correlation")

## show few entries
toy[1:10,1:5]

if (do_save) {
  saveRDS(toy, file = file.path(PATH,"toy_dataset2.rds"))
  write.csv(toy, file = file.path(PATH,"toy_dataset2.csv"), row.names = FALSE)
}