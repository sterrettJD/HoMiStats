test_that("check_proportional", {
  df <- data.frame(feature_1=c(1, 0.2, 0.4),
                   feature_2=c(0, 0.6, 0.5),
                   feature_3=c(0, 0.199999, 0))
  expect_error(check_proportional(df,
                                  "The following sample does not sum to 1: 3"))
})


