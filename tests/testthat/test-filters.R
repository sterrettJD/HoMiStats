test_that("filter_low_abundance function filters low abundance features", {
    # Create a sample dataframe with relative abundance data
    df <- data.frame(
        sample1 = c(0.1, 0.2, 0.3, 0.4),
        sample2 = c(0.2, 0.3, 0.4, 0.5),
        sample3 = c(0.3, 0.4, 0.5, 0.6)
    )
    rownames(df) <- c("feature1", "feature2", "feature3", "feature4")

    # Set the threshold for filtering low abundance features
    threshold <- 0.3

    # Call the filter_low_abundance function
    filtered_df <- filter_low_abundance_by_mean(df, threshold)

    # Check that the filtered dataframe has the correct number of rows
    expect_equal(nrow(filtered_df), 3)

    # Check that the filtered dataframe has the correct row names
    expect_equal(rownames(filtered_df), c("feature2", "feature3", "feature4"))

    # Check that the filtered dataframe has the correct values
    expect_equal(filtered_df$sample1, c(0.2, 0.3, 0.4))
    expect_equal(filtered_df$sample2, c(0.3, 0.4, 0.5))
    expect_equal(filtered_df$sample3, c(0.4, 0.5, 0.6))
})
