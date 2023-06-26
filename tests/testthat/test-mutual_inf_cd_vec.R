test_that("mutual_inf_cd_vec returns expected result", {

  data(mutual_info_df)

  result <- mutual_inf_cd_vec(mutual_info_df$Zc_XdYd, mutual_info_df$Xd)

  expect_length(result, 1L)
  expect_type(result, "double")
  expect_equal(result, 0.128029, tolerance = 0.00001)
})


test_that("mutual_inf_cc_vec issues error messages when vector sizes are different", {

  data(mutual_info_df)

  expect_error(  mutual_inf_cd_vec(mutual_info_df$Zc_XdYd[-1],
                                  mutual_info_df$Xd) )
})


test_that("mutual_inf_cc_vec issues error messages when the value of k is too large", {

  data(mutual_info_df)


  expect_error(  mutual_inf_cd_vec(mutual_info_df$Zc_XdYd,
                                  mutual_info_df$Xd, k=101) )
})
