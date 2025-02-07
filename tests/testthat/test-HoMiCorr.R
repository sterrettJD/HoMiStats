test_that("HoMiCorr works", {
    mtx <- data.frame(a1=c(0.1, 0.0, 0.0, 0.0),
                      b1=c(0.5, 0.5, 0.5, 0.4),
                      c1=c(0.4, 0.5, 0.0, 0.0),
                      d1=c(0.0, 0.0, 0.5, 0.6))
    host <- data.frame(a2=c(0.5, 0.0, 0.0, 0.0),
                      b2=c(0.5, 0.5, 0.5, 0.4),
                      c2=c(0.4, 0.5, 0.0, 0.0),
                      d2=c(0.0, 0.0, 0.5, 0.6))
    row.names(mtx) <- paste0("sample_", seq_len(4))
    row.names(host) <- paste0("sample_", seq_len(4))
    metadata <- data.frame(SampleID=paste0("sample_", seq_len(4)),
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
    row.names(mtx) <- paste0("sample_", seq_len(4))
    row.names(host) <- paste0("sample_", seq_len(4))

    expect_error(actual <- run_HoMiCorr(mtx, host,
                                           reg.method="zibr",
                                           show_progress=FALSE),
                 "Columns found in both datasets: a, b, c, d")

})

test_that("prepare_data processes inputs correctly", {
    mtx <- data.frame(a1=c(0.1, 0.0, 0.2, 0.3),
                      b1=c(0.5, 0.5, 0.5, 0.4))
    row.names(mtx) <- paste0("sample_", seq_len(4))
    host <- data.frame(a2 = c(0.5, 0.0, 0.0, 0.0),
                       b2 = c(0.5, 0.5, 0.5, 0.4))
    rownames(host) <- paste0("sample_", seq_len(4))

    metadata <- data.frame(SampleID=paste0("sample_", seq_len(4)),
                           phenotype=c(0,0,1,1))

    # Case: No metadata, no sampleID
    expect_no_error(prepped <- prepare_data(mtx, host))
    expect_equal(ncol(prepped), ncol(mtx) + ncol(host))
    expect_equal(rownames(prepped), rownames(mtx))

    # Case: Metadata with valid sampleID
    expect_no_error(prepped <- prepare_data(mtx, host, covariates="phenotype",
                                            metadata, sampleID="SampleID"))
    expect_equal(ncol(prepped), ncol(mtx) + ncol(host) + 1)
    expect_equal(rownames(prepped), metadata$SampleID) # Ensure rownames align

    # Case: SampleID column missing
    expect_error(prepare_data(mtx, host, covariates="phenotype",
                              metadata=metadata, sampleID="MissingID"),
                 "sampleID column 'MissingID' not found in metadata.")

    # Case: Mismatched row names
    bad_host <- host
    rownames(bad_host) <- paste0("other_", seq_len(4))
    expect_error(prepare_data(mtx, bad_host, covariates="phenotype",
                              metadata=metadata, "SampleID"),
                 paste0("mtx and host rownames do not match. ",
                        "Please only make sure all samples in one dataset ",
                        "exist in the other."))

    bad_metadata <- metadata
    bad_metadata$SampleID <- paste0("other_", seq_len(4))
    expect_error(prepare_data(mtx, host, covariates="phenotype",
                              metadata=bad_metadata, "SampleID"),
                 "Row names of `mtx` do not match `metadata\\$SampleID`.")
})


test_that("generate_feature_combinations works as expected", {
    mtx <- data.frame(a1 = c(0.1, 0.0, 0.0, 0.0),
                       b1 = c(0.5, 0.5, 0.5, 0.4))
    host <- data.frame(a2 = c(0.5, 0.0, 0.0, 0.0),
                        b2 = c(0.5, 0.5, 0.5, 0.4))

    # Generate feature combinations
    feature_combos <- generate_feature_combinations(mtx, host)

    expect_equal(length(feature_combos), 6)
    expect_equal(length(feature_combos[[1]]), 2)
    expect_equal(feature_combos[[1]], c(1, 2))
})

test_that("initialize_model_summaries creates correct structure", {
    feature.combos <- list(c("feature1", "feature2"), c("feature3", "feature4"))

    # Test for gamlss
    gamlss_summary <- initialize_model_summaries(feature.combos, "gamlss")
    expect_equal(nrow(gamlss_summary), length(feature.combos))
    expect_equal(colnames(gamlss_summary),
                 c("parameter", "term", "estimate",
                   "std.error", "statistic", "p.value", "feature"))

    # Test for zibr
    zibr_summary <- initialize_model_summaries(feature.combos, "zibr")
    expect_equal(nrow(zibr_summary), length(feature.combos))
    expect_equal(colnames(zibr_summary),
                 c("parameter", "term", "estimate",
                   "p.value", "joint.p", "feature"))

    # Test for lm
    lm_summary <- initialize_model_summaries(feature.combos, "lm")
    expect_equal(nrow(lm_summary), length(feature.combos))
    expect_equal(colnames(lm_summary),
                 c("term", "estimate", "std.error",
                   "statistic", "p.value", "feature"))

    # Test for lmer
    lmer_summary <- initialize_model_summaries(feature.combos, "lmer")
    expect_equal(nrow(lmer_summary), length(feature.combos))
    expect_equal(colnames(lmer_summary),
                 c("effect", "group", "term", "estimate",
                   "std.error", "statistic", "df", "p.value", "feature"))
})


test_that("stop_parallel_processing works as expected", {
    # Create a mock parallel cluster object
    cl <- parallel::makeCluster(2)

    # Check that the cluster is running by calling a simple function on the cluster
    cluster_alive <- tryCatch({
        parallel::clusterCall(cl, function() { return("alive") })
        TRUE
    }, error = function(e) {
        FALSE
    })
    expect_true(cluster_alive)

    # Test case where show_progress is FALSE
    expect_silent(stop_parallel_processing(cl, show_progress = FALSE))

    # Verify that the cluster is stopped by trying to interact with it
    cluster_alive_after <- tryCatch({
        parallel::clusterCall(cl, function() { return("alive") })
        TRUE
    }, error = function(e) {
        FALSE
    })
    expect_false(cluster_alive_after)

    # Test case where show_progress is TRUE
    # Create a mock progress bar object to test
    cl <- parallel::makeCluster(2)
    pb <- utils::txtProgressBar(min = 0, max = 100, style = 3)
    expect_silent(stop_parallel_processing(cl, show_progress = TRUE, pb))
})

test_that("adjust_p_values works for ZIBR with zero_prop_from_formula", {
    mod.summaries <- data.frame(term = c("(Intercept)", "GeneX"),
                                parameter = c("beta", "beta"),
                                p.value = c(0.01, 0.002),
                                joint.p = c(0.02, 0.005))

    # Test with ZIBR and zero_prop_from_formula = TRUE, using FDR adjustment
    adjusted <- adjust_p_values(mod.summaries, "zibr",
                                zero_prop_from_formula = TRUE,
                                padj = "fdr")

    # Ensure that only the non-intercept terms are considered for adjusted p-values
    expect_equal(sum(!is.na(adjusted$q)), 1)  # Exclude the intercept term
    expect_true(adjusted[2, "q"] <= 0.05)  # Check that adjusted p-values are <= 0.05
})

test_that("adjust_p_values works for ZIBR without zero_prop_from_formula", {
    mod.summaries <- data.frame(term = c("(Intercept)", "GeneX"),
                                parameter = c("beta", "beta"),
                                p.value = c(0.01, 0.002),
                                joint.p = c(0.02, 0.005))

    # Test with ZIBR but zero_prop_from_formula = FALSE, using FDR adjustment
    adjusted <- adjust_p_values(mod.summaries, "zibr",
                                zero_prop_from_formula = FALSE,
                                padj = "fdr")

    # Ensure that only the non-intercept terms are considered for adjusted p-values
    expect_equal(sum(!is.na(adjusted$q)), 1)  # Exclude the intercept term
    expect_true(adjusted[2, "q"] <= 0.05)  # Check that adjusted p-values are <= 0.05
})

test_that("adjust_p_values works for other regression methods", {
    mod.summaries <- data.frame(term = c("(Intercept)", "GeneX"),
                                parameter = c("beta", "beta"),
                                p.value = c(0.01, 0.002),
                                joint.p = c(0.02, 0.005))

    # Test for "lm" regression method using Bonferroni adjustment
    adjusted <- adjust_p_values(mod.summaries, "lm",
                                zero_prop_from_formula = FALSE,
                                padj = "bonferroni")

    # Ensure that only the non-intercept terms are considered for adjusted p-values
    expect_equal(sum(!is.na(adjusted$q)), 1)  # Exclude the intercept term
    expect_true(adjusted[2, "q"] <= 0.05)  # Check that adjusted p-values are <= 0.05
})
