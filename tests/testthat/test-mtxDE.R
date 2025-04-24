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

test_that("map zibr term names", {
    unmapped <- rep(c("intersept", "var1", "var2"), 10)

    mapped <- map_zibr_termnames(unmapped, c("phenotype", "treatment"))
    expect_equal(mapped,
                 rep(c("intercept", "phenotype", "treatment"), 10)
    )

    mapped <- map_zibr_termnames(unmapped, c("treatment", "phenotype"))
    expect_equal(mapped,
                 rep(c("intercept", "treatment", "phenotype"), 10)
    )
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
                 paste0("The following rows contains a value of one ",
                 "which the zero-inflated beta regression cannot handle: ",
                 "(1, 2|2, 1)"))

    no.error.df <- data.frame(a=c(0.1, 0.2),
                              b=c(0.9, 0.8))
    expect_no_error(check_for_ones(no.error.df))
})

test_that("check_data_mtxDE works", {
  feature.table <- data.frame(a=c(0.1, 0.2),
                              b=c(0.1, 0))

  unmatched.dna.table <- data.frame(c=c(0.2, 0.4),
                                    d=c(0.1, 0))

  metadata <- data.frame(SampleID = paste0("sample_", seq_len(2)),
                         phenotype = c(0,1))

  row.names(feature.table) <- paste0("sample_", seq_len(2))
  row.names(unmatched.dna.table) <- paste0("sample_", seq_len(2))

  expect_error(check_data_mtxDE(feature.table = feature.table,
                                dna.table = unmatched.dna.table,
                                metadata = metadata,
                                sampleID = "SampleID"),
               "None of the colnames of `dna.table` match `feature.table`")

  good.dna.table <- data.frame(a=c(0.2, 0.4),
                               b=c(0.1, 0))
  row.names(good.dna.table) <- paste0("sample_", seq_len(2))

  expect_no_error(check_data_mtxDE(feature.table = feature.table,
                                   dna.table = good.dna.table,
                                   metadata = metadata,
                                   sampleID = "SampleID"))
})

test_that("filter_tables_by_shared_columns correctly filters columns", {

  dna.table <- data.frame(a = c(1, 2), b = c(3, 4), c = c(5, 6))
  feature.table <- data.frame(a = c(7, 8), b = c(9, 10))  # Missing column 'c'
  expect_warning(
    result <- filter_tables_by_shared_columns(dna.table, feature.table,
                                              "dna.table", "feature.table"),
    "The following feature\\(s\\) were removed from `dna.table`"
  )
  filtered.dna.table <- result$dna.table
  filtered.feature.table <- result$feature.table

  expected.dna.table <- dna.table[, c("a", "b"), drop = FALSE]
  expected.feature.table <- feature.table[, c("a", "b"), drop = FALSE]

  expect_equal(filtered.dna.table, expected.dna.table)
  expect_equal(filtered.feature.table, expected.feature.table)

})

test_that(".prepare_data_mtxDE works", {
  feature.table <- data.frame(a=c(0.1, 0.0, 0.0, 0.0),
                              b=c(0.5, 0.5, 0.5, 0.4),
                              c=c(0.4, 0.5, 0.0, 0.0),
                              d=c(0.0, 0.0, 0.5, 0.6))
  row.names(feature.table) <- paste0("sample_", seq_len(4))
  metadata <- data.frame(SampleID=paste0("sample_", seq_len(4)),
                         phenotype=c(0,0,1,1),
                         participant=c(0,1,0,1),
                         timepoint=c(0,0,1,1))
  dna.table <- feature.table

  expected_data <- data.frame(
    Row.names=paste0("sample_", seq_len(4)),
    a=c(0.1, 0.0, 0.0, 0.0),
    b=c(0.5, 0.5, 0.5, 0.4),
    c=c(0.4, 0.5, 0.0, 0.0),
    d=c(0.0, 0.0, 0.5, 0.6),
    phenotype=c(0,0,1,1),
    a_mgx=c(0.1, 0.0, 0.0, 0.0),
    b_mgx=c(0.5, 0.5, 0.5, 0.4),
    c_mgx=c(0.4, 0.5, 0.0, 0.0),
    d_mgx=c(0.0, 0.0, 0.5, 0.6)
    )

  data <- .prepare_data_mtxDE(feature.table, metadata, " ~ phenotype",
                      "SampleID", NULL, dna.table)
  data$Row.names <- as.character(data$Row.names)
  expect_equal(data, expected_data)

  # testing the case where mismatch in dna.table
  dna.table$b <- NULL
  expected_data2 <- expected_data
  expected_data2$b <- NULL
  expected_data2$b_mgx <- NULL
  expect_warning(
    data2 <- .prepare_data_mtxDE(feature.table, metadata, " ~ phenotype",
                              "SampleID", NULL, dna.table))
  data2$Row.names <- as.character(data2$Row.names)
  expect_equal(data2, expected_data2)

  # testing the case where there's no dna.table
  expected_data2 <- expected_data[, !grepl("mgx$", colnames(expected_data))]
  data3 <- .prepare_data_mtxDE(feature.table, metadata, " ~ phenotype",
                              "SampleID", NULL)
  data3$Row.names <- as.character(data3$Row.names)
  expect_equal(data3, expected_data2)
})

test_that(".add_dna_to_formula correctly adds dna col", {
  feature.table <- data.frame(a=c(0.1, 0.0, 0.0, 0.0),
                              b=c(0.5, 0.5, 0.5, 0.4),
                              c=c(0.4, 0.5, 0.0, 0.0),
                              d=c(0.0, 0.0, 0.5, 0.6))
  row.names(feature.table) <- paste0("sample_", seq_len(4))
  metadata <- data.frame(SampleID=paste0("sample_", seq_len(4)),
                         phenotype=c(0,0,1,1),
                         participant=c(0,1,0,1),
                         timepoint=c(0,0,1,1))
  dna.table <- feature.table
  formula <- " ~ phenotype"
  data <- .prepare_data_mtxDE(feature.table, metadata, formula,
                              "SampleID", NULL, dna.table)
  fixed.vars <- all.vars(as.formula(formula))
  col <- "a"
  formula <- paste0(col, formula)

  # successfully add to fixed.vars
  expected.fixed.vars <- c(fixed.vars, "a_mgx")

  result <- .add_dna_to_formula(data, col, formula, fixed.vars,
                                reg.method = "zibr")
  actual.fixed.vars <- result$fixed.vars
  expect_equal(actual.fixed.vars, expected.fixed.vars)
  expect_equal(result$formula, formula) # shouldn't add to formula

  # successfully add to formula
  expected.formula <- as.formula("a ~ phenotype + a_mgx")
  result <- .add_dna_to_formula(data, col, formula, fixed.vars,
                                reg.method = "lm")
  actual.formula <- as.formula(result$formula)
  expect_equal(actual.formula, expected.formula)
  expect_equal(result$fixed.vars, fixed.vars) # shouldn't add to fixed.vars

  # has an error if feature doesn't exist
  col <- "nonexistent_feature"
  expect_error(
    .add_dna_to_formula(data, col, formula, fixed.vars, reg.method = "zibr"),
    "The following DNA feature was not found in metadata: nonexistent_feature_mgx"
  )

  expect_no_error(
      # lm can run with the formula
      mod <- lm(as.formula(actual.formula),
                data=data)
  )

})

test_that("run_mtxDE works", {
    feature.table <- data.frame(a=c(0.1, 0.0, 0.0, 0.0),
                                b=c(0.5, 0.5, 0.5, 0.4),
                                c=c(0.4, 0.5, 0.0, 0.0),
                                d=c(0.0, 0.0, 0.5, 0.6))
    row.names(feature.table) <- paste0("sample_", seq_len(4))
    metadata <- data.frame(SampleID=paste0("sample_", seq_len(4)),
                           phenotype=c(0,0,1,1),
                           participant=c(0,1,0,1),
                           timepoint=c(0,0,1,1))

    expected.zibr.nolong <- read.csv(
                            "data/expected_mtxDE_results_zibrnolong.csv")
    expected.zibr.long <- read.csv("data/expected_mtxDE_results_zibrlong.csv")
    expected.gamlss <- read.csv("data/expected_mtxDE_results_gamlss.csv")

    zibr.nolong.actual <-suppressWarnings(run_mtxDE("phenotype",
                                                    feature.table,
                                                    metadata,
                                                    sampleID="SampleID",
                                                    show_progress=FALSE))
    expect_equal(zibr.nolong.actual$term,
                 rep(c("intercept", "phenotype"), ncol(feature.table)*2))
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
                                          zibr_time_ind="timepoint",
                                          show_progress=FALSE)
                            )

    expect_equal(zibr.long.actual$term,
                 rep(c("intercept", "phenotype"), ncol(feature.table)*2))
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
                                      reg.method="gamlss",
                                      show_progress=FALSE)
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


test_that("ZIBR timepoint errors", {
    feature.table <- data.frame(a=c(0.1, 0.0, 0.0, 0.0),
                                b=c(0.5, 0.5, 0.5, 0.4),
                                c=c(0.4, 0.5, 0.0, 0.0),
                                d=c(0.0, 0.0, 0.5, 0.6))
    row.names(feature.table) <- paste0("sample_", seq_len(4))
    metadata <- data.frame(SampleID=paste0("sample_", seq_len(4)),
                           phenotype=c(0,0,1,1),
                           participant=c(0,1,0,1),
                           unique_participants=seq_len(4),
                           timepoint=c(0,0,1,1))
    # Longitudinal but no timepoint column
    expect_error(run_mtxDE("phenotype + (1|participant)",
                           feature.table,
                           metadata, sampleID="SampleID",
                           zibr_time_ind=NULL,
                           show_progress=FALSE),
                 paste0("A timepoint column is necessary ",
                        "if there are longitudinal samples.")
                )
    # Not longitudinal with no timepoint column
    expect_no_error(suppressWarnings(
                run_mtxDE("phenotype + (1|unique_participants)",
                           feature.table,
                           metadata, sampleID="SampleID",
                           zibr_time_ind=NULL,
                          show_progress=FALSE)
                        )
                    )
})

test_that("run_mtxDE works linear models", {
    feature.table <- data.frame(a=c(0.1, 0.0, 0.0, 0.0),
                                b=c(0.5, 0.5, 0.5, 0.4),
                                c=c(0.4, 0.5, 0.0, 0.0),
                                d=c(0.0, 0.0, 0.5, 0.6))
    row.names(feature.table) <- paste0("sample_", seq_len(4))
    metadata <- data.frame(SampleID=paste0("sample_", seq_len(4)),
                           phenotype=c(0,0,1,1),
                           participant=c(0,1,0,1),
                           timepoint=c(0,0,1,1))

    expected.lm <- read.csv("data/expected_mtxDE_results_lm.csv")
    expected.lmer <- read.csv("data/expected_mtxDE_results_lmer.csv")

    lm.actual <- run_mtxDE("phenotype",
                           feature.table,
                           metadata,
                           reg.method="lm",
                           sampleID="SampleID",
                           show_progress=FALSE)

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
                            sampleID="SampleID",
                            show_progress=FALSE)

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
                             sampleID="SampleID",
                             show_progress=FALSE))
})

test_that("run_mtxDE works with multiple cores", {
    feature.table <- data.frame(a=c(0.1, 0.0, 0.0, 0.0),
                                b=c(0.5, 0.5, 0.5, 0.4),
                                c=c(0.4, 0.5, 0.0, 0.0),
                                d=c(0.0, 0.0, 0.5, 0.6))
    row.names(feature.table) <- paste0("sample_", seq_len(4))
    metadata <- data.frame(SampleID=paste0("sample_", seq_len(4)),
                           phenotype=c(0,0,1,1),
                           participant=c(0,1,0,1),
                           timepoint=c(0,0,1,1))


    expect_no_error(run_mtxDE("phenotype + (1|participant)",
                              feature.table,
                              metadata,
                              reg.method="lmer",
                              sampleID="SampleID",
                              ncores=2,
                              show_progress=FALSE))
})



test_that("run_mtxDE warns about undetected features", {
    feature.table <- data.frame(a=c(0.0, 0.0, 0.0, 0.0),
                                b=c(0.0, 0.0, 0.0, 0.0),
                                c=c(0.0, 0.0, 0.0, 0.0),
                                d=c(0.0, 0.0, 0.0, 0.0))
    row.names(feature.table) <- paste0("sample_", seq_len(4))
    metadata <- data.frame(SampleID=paste0("sample_", seq_len(4)),
                           phenotype=c(0,0,1,1),
                           participant=c(0,1,0,1),
                           timepoint=c(0,0,1,1))


    expect_warning(run_mtxDE("phenotype",
                              feature.table,
                              metadata,
                              reg.method="gamlss",
                              sampleID="SampleID",
                              ncores=2,
                              show_progress=FALSE),
                   "Undetected features",
                   perl=TRUE)

    expect_warning(run_mtxDE("phenotype",
                             feature.table,
                             metadata,
                             reg.method="zibr",
                             sampleID="SampleID",
                             ncores=2,
                             show_progress=FALSE),
                   "Undetected features",
                   perl=TRUE)
})


test_that("run_mtxDE works with DNA", {
    feature.table <- data.frame(a=c(0.1, 0.0, 0.0, 0.0),
                                b=c(0.5, 0.5, 0.5, 0.4),
                                c=c(0.4, 0.5, 0.0, 0.0),
                                d=c(0.0, 0.0, 0.5, 0.6))
    dna.table <- feature.table
    row.names(feature.table) <- paste0("sample_", seq_len(4))
    row.names(dna.table) <- paste0("sample_", seq_len(4))
    metadata <- data.frame(SampleID=paste0("sample_", seq_len(4)),
                           phenotype=c(0,0,1,1),
                           participant=c(0,1,0,1),
                           timepoint=c(0,0,0,0))


    # lm based
    expect_no_error(res <- run_mtxDE("phenotype",
                             feature.table,
                             metadata,
                             dna.table=dna.table,
                             reg.method="lm",
                             sampleID="SampleID",
                             ncores=2,
                             show_progress=FALSE))

    expect_equal(nrow(dplyr::filter(res, term=="a_mgx")), 1)
    expect_true(res[res$term=="d_mgx", "q"] < 0.05)
    expect_true(res[res$term=="phenotype" & res$feature=="d", "q"] > 0.05)

    # gamlss based results
    set.seed(1234)
    N <- 25
    a <- runif(N, max=1-(1E-4))
    b <- runif(N, max=1-(1E-4))
    feature.table <- data.frame(a=a,
                                b=b)
    dna.table <- data.frame(a=a,
                            b=b)
    row.names(feature.table) <- paste0("sample_", seq_len(N))
    row.names(dna.table) <- paste0("sample_", seq_len(N))
    metadata <- data.frame(SampleID=paste0("sample_", seq_len(N)),
                           phenotype=rbinom(N, 1, 0.5))


    # gamlss based
    res <- HoMiStats::run_mtxDE("phenotype",
                                feature.table,
                                metadata,
                                dna.table=dna.table,
                                reg.method="gamlss",
                                sampleID="SampleID",
                                ncores=2,
                                show_progress=FALSE)
    expect_no_error(res)
    expect_equal(nrow(
                dplyr::filter(res,
                              term %in% c("a_mgx", "phenotype", "b_mgx"))),
                 4)
    # mgx should be very predictive, but phenotype shouldn't have an effect
    expect_true(res[res$term=="a_mgx", "q"] < 0.05)
    expect_true(res[res$term=="b_mgx", "q"] < 0.05)
    expect_true(all(res[res$term=="phenotype", "q"] > 0.05))

})
