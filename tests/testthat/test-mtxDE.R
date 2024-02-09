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
    row.names(feature.table) <- paste0("sample_", 1:4)
    metadata <- data.frame(SampleID=paste0("sample_", 1:4),
                           phenotype=c(0,0,1,1),
                           participant=c(0,1,0,1),
                           timepoint=c(0,0,1,1))

    expected.zibr.nolong <- read.csv("data/expected_mtxDE_results_zibrnolong.csv")
    expected.zibr.long <- read.csv("data/expected_mtxDE_results_zibrlong.csv")
    expected.gamlss <- read.csv("data/expected_mtxDE_results_gamlss.csv")

    zibr.nolong.actual <-suppressWarnings(run_mtxDE("phenotype",
                                                    feature.table,
                                                    metadata,
                                                    sampleID="SampleID"))

    expect_equal(zibr.nolong.actual$p.value,
                 expected.zibr.nolong$p.value,
                 tolerance=1e-3)
    expect_equal(zibr.nolong.actual$estimate,
                 expected.zibr.nolong$estimate,
                 tolerance=5 # some of these values are big and variable
                 # and the p value/CI matters more
                 )
    expect_equal(colnames(zibr.nolong.actual),
                 colnames(expected.zibr.nolong))

    zibr.long.actual <- suppressWarnings(
                                run_mtxDE("phenotype + (1|participant)",
                                          feature.table,
                                          metadata, sampleID="SampleID",
                                          zibr_time_ind="timepoint")
                            )

    expect_equal(zibr.long.actual$p.value,
                 expected.zibr.long$p.value,
                 tolerance=1e-2)
    expect_equal(zibr.long.actual$estimate,
                 expected.zibr.long$estimate,
                 tolerance=5 # some of these values are big and variable
                 # and the p value/CI matters more
                 )
    expect_equal(colnames(zibr.long.actual),
                 colnames(expected.zibr.long))


    # gamlss
    gamlss.actual <- suppressWarnings(run_mtxDE("phenotype", feature.table,
                                      metadata, sampleID="SampleID",
                                      reg.method="gamlss")
                                     )
    expect_equal(gamlss.actual$p.value,
                 expected.gamlss$p.value,
                 tolerance=1e-2)
    expect_equal(gamlss.actual$p.value,
                 expected.gamlss$p.value,
                 tolerance=1e-3)
    expect_equal(colnames(gamlss.actual),
                 colnames(expected.gamlss))

})


test_that("run_mtxDE works linear models", {
    feature.table <- data.frame(a=c(0.1, 0.0, 0.0, 0.0),
                                b=c(0.5, 0.5, 0.5, 0.4),
                                c=c(0.4, 0.5, 0.0, 0.0),
                                d=c(0.0, 0.0, 0.5, 0.6))
    row.names(feature.table) <- paste0("sample_", 1:4)
    metadata <- data.frame(SampleID=paste0("sample_", 1:4),
                           phenotype=c(0,0,1,1),
                           participant=c(0,1,0,1),
                           timepoint=c(0,0,1,1))

    expected.lm <- read.csv("data/expected_mtxDE_results_lm.csv")
    expected.lmer <- read.csv("data/expected_mtxDE_results_lmer.csv")

    lm.actual <- run_mtxDE("phenotype",
                           feature.table,
                           metadata,
                           reg.method="lm",
                           sampleID="SampleID")

    expect_equal(lm.actual$p.value,
                 expected.lm$p.value,
                 tolerance=1e-3)
    expect_equal(lm.actual$estimate,
                 expected.lm$estimate,
                 tolerance=1 # some of these values are big and variable
                 # and the p value/CI matters more
    )
    expect_equal(colnames(lm.actual),
                 colnames(expected.lm))

    lmer.actual <- run_mtxDE("phenotype + (1|participant)",
                            feature.table,
                            metadata,
                            reg.method="lmer",
                            sampleID="SampleID")

    expect_equal(lmer.actual$p.value,
                 expected.lmer$p.value,
                 tolerance=1e-3)
    expect_equal(lmer.actual$estimate,
                 expected.lmer$estimate,
                 tolerance=1 # some of these values are big and variable
                 # and the p value/CI matters more
    )
    expect_equal(colnames(lmer.actual),
                 colnames(expected.lmer))

    # Make sure there aren't errors from values >=1 that don't fit the beta reg
    feature.table.bigvals <- feature.table*100
    expect_no_error(run_mtxDE("phenotype + (1|participant)",
                             feature.table.bigvals,
                             metadata,
                             reg.method="lmer",
                             sampleID="SampleID"))
})

