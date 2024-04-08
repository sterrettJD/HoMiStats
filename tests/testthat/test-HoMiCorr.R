test_that("HoMiCorr works", {
    mtx <- data.frame(a1=c(0.1, 0.0, 0.0, 0.0),
                      b1=c(0.5, 0.5, 0.5, 0.4),
                      c1=c(0.4, 0.5, 0.0, 0.0),
                      d1=c(0.0, 0.0, 0.5, 0.6))
    host <- data.frame(a2=c(0.5, 0.0, 0.0, 0.0),
                      b2=c(0.5, 0.5, 0.5, 0.4),
                      c2=c(0.4, 0.5, 0.0, 0.0),
                      d2=c(0.0, 0.0, 0.5, 0.6))
    row.names(mtx) <- paste0("sample_", 1:4)
    row.names(host) <- paste0("sample_", 1:4)
    metadata <- data.frame(SampleID=paste0("sample_", 1:4),
                           phenotype=c(0,0,1,1),
                           participant=c(0,1,0,1),
                           timepoint=c(0,0,1,1))



    expect_no_error(actual <- run_HoMiCorr(mtx, host,
                                 reg.method="zibr",
                                 show_progress=FALSE))
    expected <- read.csv("data/expected_homicorr_zibr.csv")
    expect_equal(actual$term,
                 expected$term)
    expect_equal(actual$p.value < 0.01,
                 expected$p.value < 0.01)
    expect_equal(actual$estimate,
                 expected$estimate,
                 tolerance=5 # some of these values are big and variable
                 # and the p value/CI matters more
    )

    expect_no_error(actual <- run_HoMiCorr(mtx, host,
                                 reg.method="zibr",
                                 covariates="(1|participant)",
                                 metadata=metadata,
                                 sampleID="SampleID",
                                 zibr_time_ind="timepoint",
                                 show_progress=FALSE))
    expected <- read.csv("data/expected_homicorr_zibr_rand.csv")
    expect_equal(actual$term,
                 expected$term)
    expect_equal(actual$p.value < 0.01,
                 expected$p.value < 0.01)
    expect_equal(actual$estimate,
                 expected$estimate,
                 tolerance=5 # some of these values are big and variable
                 # and the p value/CI matters more
    )


    expect_no_error(actual <- run_HoMiCorr(mtx, host,
                                 reg.method="gamlss",
                                 show_progress=FALSE))
    expected <- read.csv("data/expected_homicorr_gamlss.csv")
    expect_equal(actual$term,
                 expected$term)
    expect_equal(actual$p.value,
                 expected$p.value,
                 tolerance=1e-3)
    expect_equal(actual$estimate,
                 expected$estimate,
                 tolerance=5 # some of these values are big and variable
                 # and the p value/CI matters more
    )


    expect_no_error(actual <- run_HoMiCorr(mtx, host,
                                 reg.method="lm",
                                 show_progress=FALSE))
    expected <- read.csv("data/expected_homicorr_lm.csv")
    expect_equal(actual$term,
                 expected$term)
    expect_equal(actual$p.value,
                 expected$p.value,
                 tolerance=1e-3)
    expect_equal(actual$estimate,
                 expected$estimate,
                 tolerance=5 # some of these values are big and variable
                 # and the p value/CI matters more
    )



    expect_no_error(actual <- run_HoMiCorr(mtx, host,
                                 covariates="(1|participant)",
                                 metadata=metadata,
                                 sampleID="SampleID",
                                 reg.method="lmer",
                                 show_progress=FALSE))
    expected <- read.csv("data/expected_homicorr_lmer.csv")
    expect_equal(actual$term,
                 expected$term)
    expect_equal(actual$p.value,
                 expected$p.value,
                 tolerance=1e-3)
    expect_equal(actual$estimate,
                 expected$estimate,
                 tolerance=5 # some of these values are big and variable
                 # and the p value/CI matters more
    )

})


test_that("HoMiCorr errors with duplicated column names", {
    mtx <- data.frame(a=c(0.1, 0.0, 0.0, 0.0),
                      b=c(0.5, 0.5, 0.5, 0.4),
                      c=c(0.4, 0.5, 0.0, 0.0),
                      d=c(0.0, 0.0, 0.5, 0.6))
    host <- data.frame(a=c(0.5, 0.0, 0.0, 0.0),
                       b=c(0.5, 0.5, 0.5, 0.4),
                       c=c(0.4, 0.5, 0.0, 0.0),
                       d=c(0.0, 0.0, 0.5, 0.6))
    row.names(mtx) <- paste0("sample_", 1:4)
    row.names(host) <- paste0("sample_", 1:4)

    expect_error(actual <- run_HoMiCorr(mtx, host,
                                           reg.method="zibr",
                                           show_progress=FALSE),
                 "a is found in both datasets.")

})
