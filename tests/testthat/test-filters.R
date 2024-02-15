test_that("filter_low_abundance_by_mean function filters low abundance features", {
    df <- data.frame(
        sample1 = c(0.1, 0.2, 0.3, 0.4),
        sample2 = c(0.2, 0.3, 0.4, 0.5),
        sample3 = c(0.3, 0.4, 0.5, 0.6)
    )
    df <- data.frame(t(df)) #originally wrote this transposed and would rather just have R flip it
    colnames(df) <- c("feature1", "feature2", "feature3", "feature4")

    threshold <- 0.3
    filtered_df <- filter_low_abundance_by_mean(df, threshold)

    # Check that the filtered dataframe has the correct number of rows
    expect_equal(ncol(filtered_df), 3)

    # Check that the filtered dataframe has the correct row names
    expect_equal(colnames(filtered_df), c("feature2", "feature3", "feature4"))

    # Check that the filtered dataframe has the correct values
    expect_null(filtered_df$feature1)
    expect_equal(filtered_df$feature2, c(0.2, 0.3, 0.4))
    expect_equal(filtered_df$feature3, c(0.3, 0.4, 0.5))
})



test_that("filter_low_variance function filters low variance features", {
    df <- data.frame(
        feature1 = c(0.1, 0.1, 0.0),
        feature2 = c(0.2, 0.3, 0.5),
        feature3 = c(0.3, 0.4, 0.8),
        feature4 = c(0.4, 0.5, 0.7)
    )
    rownames(df) <- c("sample1", "sample2", "sample3")

    threshold <- 0.01
    filtered_df <- filter_low_variance(df, threshold)

    # Check that the filtered dataframe has the correct number of rows and columns
    expect_equal(ncol(filtered_df), 3)

    # Check that the filtered dataframe has the correct values
    expect_null(filtered_df$feature1)
    expect_equal(filtered_df$feature2, c(0.2, 0.3, 0.5))
    expect_equal(filtered_df$feature3, c(0.3, 0.4, 0.8))
})
