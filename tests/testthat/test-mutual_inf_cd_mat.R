test_that("mutual_inf_cd_mat returns expected result", {

  #' data(mutual_info_df)
  #'
  #' M <- cbind(mutual_info_df$Xd, mutual_info_df$Yd)
  #' mutual_inf_cd_mat(mutual_info_df$Zc_XdYdWd, M)
  #' ## 0.1070804 0.1041177


  data(mutual_info_df)

  M <- cbind(mutual_info_df$Xd, mutual_info_df$Yd)
  result <- mutual_inf_cd_mat(mutual_info_df$Zc_XdYdWd, M)

  expect_length(result, 2L)
  expect_type(result, "double")
  expect_equal(result, c(0.1070804, 0.1041177), tolerance = 0.00001)

})

test_that("mutual_inf_cd_mat issues error messages when needed", {

  data(mutual_info_df)

  M <- cbind(mutual_info_df$Xd, mutual_info_df$Yd)
  expect_error( mutual_inf_cd_mat(mutual_info_df$Zd_XdYd, M, k=3) )

  M <- cbind(mutual_info_df$Xd, mutual_info_df$Yd)
  expect_error( mutual_inf_cd_mat(mutual_info_df$Zc_XcYc[-1], M, k=3) )
})
