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



    expect_no_error(run_HoMiCorr(mtx, host,
                                 reg.method="zibr",
                                 sampleID="SampleID",
                                 show_progress=FALSE))
    expect_no_error(run_HoMiCorr(mtx, host,
                                 reg.method="gamlss",
                                 sampleID="SampleID",
                                 show_progress=FALSE))
    expect_no_error(run_HoMiCorr(mtx, host,
                                 sampleID="lm",
                                 show_progress=FALSE))
    expect_no_error(run_HoMiCorr(mtx, host,
                                 covariates="(1|participant)",
                                 metadata=metadata,
                                 sampleID="SampleID",
                                 reg.method="lmer",
                                 show_progress=FALSE))
})
