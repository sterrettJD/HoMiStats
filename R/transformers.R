#' Transform feature table for analysis
#' @description transforms feature table for differential expression or co-occurrence analysis
#' @param feature.table A dataframe, where rows are samples, and columns are genes/features. Row names should be sample IDs.
#' @return The transformed table
#' @export
#' @importFrom dplyr mutate_all
#'
transform_feature_table <- function(feature.table, method){
    mutate_all(feature.table, method)
}

#' Make sure data are proportional within a row
#' @description Some transformations are based on proportional data, where samples sum to 1. This checks that that's the case.
#' @param feature.table A dataframe, where rows are samples, and columns are genes/features. Row names should be sample IDs.
#' @return Nothing
#' @export
#'
check_proportional <- function(feature.table, tolerance=1e-3){
    sample.sums <- rowSums(feature.table)
    rows.not.prop <- which((sample.sums < (1-tolerance)) | (sample.sums > (1+tolerance)))
    rownames.not.prop <- rownames(feature.table)[rows.not.prop]

    if(length(rownames.not.prop) > 0){
        stop(paste("The following sample does not sum to 1:",
                   rownames.not.prop))
    }
}

