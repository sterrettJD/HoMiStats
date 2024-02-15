#' Filter low abundance features from relative abundance data
#' @description Filters low abundance features with means below a certain threshold.
#' @param df A dataframe containing relative abundance data.
#' @param threshold The threshold for filtering low abundance features.
#' @return A dataframe with low abundance features filtered out.
#' @examples
#' df <- data.frame(
#'   feature1 = c(0.1, 0.2, 0.3, 0.4),
#'   feature2 = c(0.2, 0.3, 0.4, 0.5),
#'   feature3 = c(0.3, 0.4, 0.5, 0.6)
#' )
#' filtered_df <- filter_low_abundance(df, 0.5)
#' print(filtered_df)
#' @export
filter_low_abundance_by_mean <- function(df, threshold) {
    # Calculate the mean abundance for each feature
    mean_abundance <- rowMeans(df)

    # Filter out features with mean abundance below the threshold
    filtered_df <- dplyr::filter(df, mean_abundance >= threshold)

    return(filtered_df)
}


