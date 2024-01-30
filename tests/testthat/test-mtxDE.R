test_that("tidy zibr", {
    dirty.mod <- list(
        logistic_est_table=data.frame(Estimate=c(0, 0),
                                      Pvalue=c(1, 1),
                                      row.names=c("intersept", "phenotype")),
        logistic_s1_est=1,
        beta_est_table=data.frame(Estimate=c(0, 0),
                                  Pvalue=c(1, 1),
                                  row.names=c("intersept", "phenotype")),
        beta_s2_est=1,
        beta_v_est=1,
        loglikelihood=500,
        joint_p=c(phenotype=1)
    )

    expected.clean.mod <- data.frame(
        parameter=c("logistic", "logistic",
                    "beta", "beta"),
        term=c("intersept", "phenotype",
               "intersept", "phenotype"),
        estimate=c(0, 0, 0, 0),
        p.value=c(1, 1, 1, 1),
        joint.p=c(NA, 1, NA, 1)
    )

    clean.mod <- tidy_zibr_results(dirty.mod)
    expect_equal(clean.mod, expected.clean.mod)
})


test_that("rand effects getter", {
    expect_null(get_random_fx(y ~ x))
    expect_null(get_random_fx(~ x))

    expect_equal(get_random_fx(y ~ x + (1|z)),
                 c("z"))
    expect_equal(get_random_fx(y ~ x*z + (1|a) + (1|b)),
                 c("a", "b"))
})
