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

test_that("check_for_ones works", {
    should.error.df <- data.frame(a=c(0, 1),
                                  b=c(1, 0))

    expect_error(check_for_ones(should.error.df),
                 "The following rows contains a value of one which the zero-inflated beta regression cannot handle: 1")
    expect_error(check_for_ones(should.error.df),
                 "The following rows contains a value of one which the zero-inflated beta regression cannot handle: 2")

    no.error.df <- data.frame(a=c(0.1, 0.2),
                              b=c(0.9, 0.8))
    expect_no_error(check_for_ones(no.error.df))
})


test_that("run_mtxDE works", {
    feature.table <- data.frame(a=c(0.1, 0.0, 0.0, 0.0),
                                b=c(0.5, 0.5, 0.5, 0.4),
                                c=c(0.4, 0.5, 0.0, 0.0),
                                d=c(0.0, 0.0, 0.5, 0.6))
    metadata <- data.frame(SampleID=1:4,
                           phenotype=c(0,0,1,1),
                           participant=c(0,1,0,1),
                           timepoint=c(0,0,1,1))

    expected.zibr.nolong <- read.csv("data/expected_mtxDE_results_zibrnolong.csv")
    expected.zibr.long <- read.csv("data/expected_mtxDE_results_zibrlong.csv")
    expected.gamlss <- read.csv("data/expected_mtxDE_results_gamlss.csv")

    expect_equal(suppressWarnings(run_mtxDE("phenotype", feature.table,
                                           metadata, sampleID="SampleID")),
                 expected.zibr.nolong,
                 tolerance=1e-3)

    expect_equal(suppressWarnings(
                    run_mtxDE("phenotype + (1|participant)", feature.table,
                              metadata, sampleID="SampleID", zibr_time_ind="timepoint")
                    ),
                 expected.zibr.long,
                 tolerance=1e-3)

    expect_equal(suppressWarnings(run_mtxDE("phenotype", feature.table,
                                            metadata, sampleID="SampleID",
                                            reg.method="gamlss")),
                 expected.gamlss,
                 tolerance=1e-3)

})
