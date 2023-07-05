test_that("cond_mutual_inf returns expected result - mat/mat", {
  set.seed(654321)
  
  data(mutual_info_df)
  M <- cbind(mutual_info_df$Xc, mutual_info_df$Yc)
  ZM <- cbind(mutual_info_df$Yc, mutual_info_df$Wc)
  tM <- t(M)
  tZM <- t(ZM)
  result <- cond_mutual_inf(mutual_info_df$Zc_XcYcWc, tM, tZM)
  
  expect_length(result, 2L)
  expect_type(result, "double")
  expect_equal(result, c(0.1171533, 0.2192397), tolerance = 0.000001)
  
})

test_that("cond_mutual_inf returns expected result - vec/mat", {
  set.seed(654321)
  
  data(mutual_info_df)
  M <- cbind(mutual_info_df$Xc, mutual_info_df$Yc)
  ZM <- cbind(mutual_info_df$Yc, mutual_info_df$Wc)
  tM <- t(M)
  tZM <- t(ZM)
  
  # Use a row of tM as a vector
  tM1 <- c(tM[1,1:ncol(tM)])
  result <- cond_mutual_inf(mutual_info_df$Zc_XcYcWc, tM1, tZM)
  
  expect_length(result, 2L)
  expect_type(result, "double")
  expect_equal(result, c(0.1171533, 0.2029773), tolerance = 0.000001)
  
  tM2 <- c(tM[2,1:ncol(tM)])
  result <- cond_mutual_inf(mutual_info_df$Zc_XcYcWc, tM2, tZM)
  
  expect_length(result, 2L)
  expect_type(result, "double")
  expect_equal(result, c( 0.035512432, 0.219239735), tolerance = 0.000001)  
  
})


test_that("cond_mutual_inf returns expected result - mat/vec", {
  set.seed(654321)
  
  data(mutual_info_df)
  M <- cbind(mutual_info_df$Xc, mutual_info_df$Yc)
  ZM <- cbind(mutual_info_df$Yc, mutual_info_df$Wc)
  tM <- t(M)
  tZM <- t(ZM)
  
  # Use a row of tM as a vector
  tZM1 <- c(tZM[1,1:ncol(tZM)])
  result <- cond_mutual_inf(mutual_info_df$Zc_XcYcWc, tM, tZM1)
  
  expect_length(result, 2L)
  expect_type(result, "double")
  expect_equal(result, c(0.11715329, 0.02990735), tolerance = 0.000001)
  
  tZM2 <- c(tZM[2,1:ncol(tZM)])
  result <- cond_mutual_inf(mutual_info_df$Zc_XcYcWc, tM, tZM2)
  
  expect_length(result, 2L)
  expect_type(result, "double")
  expect_equal(result, c(0.2029773, 0.2192397), tolerance = 0.000001)  
  
})

test_that("cond_mutual_inf  returns expected result - vec/vec", {
  set.seed(654321)
  
  data(mutual_info_df)
  
  result <- cond_mutual_inf (mutual_info_df$Zc_XcYc,
                             mutual_info_df$Xc, mutual_info_df$Yc)
  expect_length(result, 1L)
  expect_type(result, "double")
  expect_equal(result, 0.2936858, tolerance = 0.00001)
  
})

test_that("cond_mutual_inf  issues error messages when vector and matrix have different sizes", {
  
  data(mutual_info_df)
  M <- cbind(mutual_info_df$Xc, mutual_info_df$Yc)
  ZM <- cbind(mutual_info_df$Yc, mutual_info_df$Wc)
  expect_error( cond_mutual_inf(mutual_info_df$Zc_XcYcWc[-1],
                                    M,
                                    ZM))
  expect_error( cond_mutual_inf(mutual_info_df$Zc_XcYcWc,
                                    M, ZM[-1,]))
  expect_error( cond_mutual_inf(mutual_info_df$Zc_XcYcWc,
                                    M, ZM[, 1]))
  expect_error( cond_mutual_inf(mutual_info_df$Zc_XcYcWc,
                                    M[-1,], ZM))
  expect_error( cond_mutual_inf(mutual_info_df$Zc_XcYcWc,
                                    M[, 1], ZM))
})


test_that("cond_mutual_inf issues error messages when the value of k is too large", {
  
  data(mutual_info_df)
  M <- cbind(mutual_info_df$Xc, mutual_info_df$Yc)
  ZM <- cbind(mutual_info_df$Yc, mutual_info_df$Wc)
  expect_error( cond_mutual_inf(mutual_info_df$Zc_XcYcWc,
                                    M,
                                    ZM, k=101))
})


test_that("cond_mutual_inf returns expected result - integer mat/mat", {
  set.seed(654321)
  
  data(mutual_info_df)
  M <- cbind(mutual_info_df$Xd, mutual_info_df$Yd)
  ZM <- cbind(mutual_info_df$Yd, mutual_info_df$Wd)
  tM <- t(M)
  tZM <- t(ZM)
  result <- cond_mutual_inf(mutual_info_df$Zc_XdYdWd, tM, tZM)
  
  expect_length(result, 2L)
  expect_type(result, "double")
  expect_equal(result, c(0.1757598, 0.1086227), tolerance = 0.000001)
  
})

test_that("cond_mutual_inf returns expected result - integer vec/mat", {
  set.seed(654321)
  
  data(mutual_info_df)
  M <- cbind(mutual_info_df$Xd, mutual_info_df$Yd)
  ZM <- cbind(mutual_info_df$Yd, mutual_info_df$Wd)
  tM <- t(M)
  tZM <- t(ZM)
  
  tM1 <- c(tM[1,1:ncol(tM)])
  result <- cond_mutual_inf(mutual_info_df$Zc_XdYdWd, tM1, tZM)
  
  expect_length(result, 2L)
  expect_type(result, "double")
  expect_equal(result, c(0.17575981, 0.19006222), tolerance = 0.000001)
  
  tM2 <- c(tM[2,1:ncol(tM)])
  result <- cond_mutual_inf(mutual_info_df$Zc_XdYdWd, tM2, tZM)
  
  expect_length(result, 2L)
  expect_type(result, "double")
  expect_equal(result, c(0.0000000, 0.1086227), tolerance = 0.000001)  

})

test_that("cond_mutual_inf returns expected result - integer mat/vec", {
  set.seed(654321)
  
  data(mutual_info_df)
  M <- cbind(mutual_info_df$Xd, mutual_info_df$Yd)
  ZM <- cbind(mutual_info_df$Yd, mutual_info_df$Wd)
  tM <- t(M)
  tZM <- t(ZM)
  
  tZM1 <- c(tZM[1,1:ncol(tZM)])
  result <- cond_mutual_inf(mutual_info_df$Zc_XdYdWd, tM, tZM1)
  
  expect_length(result, 2L)
  expect_type(result, "double")
  expect_equal(result, c(0.1757598, 0.0000000), tolerance = 0.000001)
  
  tZM2 <- c(tZM[2,1:ncol(tZM)])
  result <- cond_mutual_inf(mutual_info_df$Zc_XdYdWd, tM, tZM2)
  
  expect_length(result, 2L)
  expect_type(result, "double")
  expect_equal(result, c(0.1904468, 0.1086227), tolerance = 0.000001)  
  
})

test_that("cond_mutual_inf returns expected result - integer vec/vec", {
  set.seed(654321)
  
  data(mutual_info_df)
  M <- cbind(mutual_info_df$Xd, mutual_info_df$Yd)
  ZM <- cbind(mutual_info_df$Yd, mutual_info_df$Wd)
  tM <- t(M)
  tZM <- t(ZM)
  
  tM1 <- c(tM[1,1:ncol(tM)])
  tZM1 <- c(tZM[1,1:ncol(tZM)])
  result <- cond_mutual_inf(mutual_info_df$Zc_XdYdWd, tM1, tZM1)
  
  expect_length(result, 1L)
  expect_type(result, "double")
  expect_equal(result, 0.1757598, tolerance = 0.00001)
  
  tM2 <- c(tM[2,1:ncol(tM)])
  tZM2 <- c(tZM[2,1:ncol(tZM)])
  result <- cond_mutual_inf(mutual_info_df$Zc_XdYdWd, tM2, tZM2)
  
  expect_length(result, 1L)
  expect_type(result, "double")
  expect_equal(result, 0.1086227, tolerance = 0.00001) 
  
})


test_that("cond_mutual_inf issues error messages when vector sizes are different", {
  
  data(mutual_info_df)
  
  expect_error(   cond_mutual_inf(mutual_info_df$Zc_XcYc[-1],
                                      mutual_info_df$Xc, mutual_info_df$Yc) )
  expect_error(   cond_mutual_in(mutual_info_df$Zc_XcYc,
                                      mutual_info_df$Xc[-1], mutual_info_df$Yc) )
  expect_error(   cond_mutual_inf(mutual_info_df$Zc_XcYc,
                                      mutual_info_df$Xc, mutual_info_df$Yc[-1]) )
  
})

 