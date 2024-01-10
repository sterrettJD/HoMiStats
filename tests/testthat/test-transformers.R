test_that("check_proportional", {
  df <- data.frame(feature_1=c(1, 0.2, 0.4),
                   feature_2=c(0, 0.6, 0.5),
                   feature_3=c(0, 0.199999, 0))
  expect_error(check_proportional(df,
                                  "The following sample does not sum to 1: 3"))
})

test_that("arcsinh works",{
    vec <- seq(from=0, to=1, by=0.2)
    exp.scale.1 <- c(0, 0.198, 0.390, 0.568, 0.732, 0.881)
    exp.scale.100 <- c(0, 3.689, 4.382, 4.787, 5.075, 5.298)
    exp.scale.100.norm <- c(0, 0.696, 0.827, 0.903, 0.957, 1)

    expect_equal(arcsinh(vec, scaling.factor=1, norm.by.scaling=F),
                 exp.scale.1,
                 tolerance=1e-3)
    expect_equal(arcsinh(vec, scaling.factor=100, norm.by.scaling=F),
                 exp.scale.100,
                 tolerance=1e-3)
    expect_equal(arcsinh(vec, scaling.factor=100, norm.by.scaling=T),
                 exp.scale.100.norm,
                 tolerance=1e-3)
})

test_that("multi arcsinh works",{
    df <- data.frame(feature_1=c(1, 0.2),
                     feature_2=c(0, 0.8))
    exp.scale.1 <- data.frame(feature_1=c(0.881, 0.198),
                              feature_2=c(0, 0.732))

    t.func <- function(x) arcsinh(x, scaling.factor=1, norm.by.scaling=F)
    expect_equal(transform_feature_table(df, transformation=t.func),
                 exp.scale.1,
                 tolerance=1e-3)

})
