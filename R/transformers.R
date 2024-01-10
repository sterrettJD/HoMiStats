#' Transform feature table for analysis
#' @description transforms feature table for differential expression or co-occurrence analysis
#' @param feature.table A dataframe, where rows are samples, and columns are genes/features. Row names should be sample IDs.
#' @param transformation Method for transforming the feature.table. String options include NOT YET IMPLEMENTED, or transformation can be a function used for transforming data
#' @return The transformed table
#' @export
#' @importFrom dplyr mutate_all
#'
transform_feature_table <- function(feature.table, transformation){
    check_proportional(feature.table, soft=TRUE)
    dplyr::mutate_all(feature.table, transformation)
}

#' Make sure data are proportional within a row
#' @description Some transformations are based on proportional data, where samples sum to 1. This checks that that's the case.
#' @param feature.table A dataframe, where rows are samples, and columns are genes/features. Row names should be sample IDs.
#' @param tolerance Numerical value for tolerance in checking if samples sum to 1.
#' @param soft Boolean. If true, this function only warns.
#' @return Nothing
#' @export
#'
check_proportional <- function(feature.table, tolerance=1e-3, soft=FALSE){
    sample.sums <- rowSums(feature.table)
    rows.not.prop <- which((sample.sums < (1-tolerance)) | (sample.sums > (1+tolerance)))
    rownames.not.prop <- rownames(feature.table)[rows.not.prop]

    if(length(rownames.not.prop) > 0){
        if(soft){
            warning(paste("The following sample does not sum to 1:",
                       rownames.not.prop,
                       "Please make sure this is intentional..."))
        } else {
            stop(paste("The following sample does not sum to 1:",
                       rownames.not.prop))
        }

    }
}

#' Perform a inverse hyperbolic sin transformation
#' @description Performs an arcsinh transformation, with scaling and adjustment. The arcsinh transformation behaves like a log transformation for large values but linearly for small values, and it can handle 0s.
#' @param vec The numerical vector to be transformed
#' @param scaling.factor A scaling factor by which to multiply the values prior to transforming. Increasing this scaling factor increases how much of the lower end of the data is treated linearly by the transformation. This is the opposite of the standard "cofactor" used with the arcsinh transformation in flow cyotmetry, as we're working with input data < 1, so we likely need to scale our values up.
#' @param norm.by.scaling Boolean defining whether transformed data should be normalized (divided) by the arcsinh of the scaling factor. This is useful if you're using a model such as the beta regression implemented in mtxDE, which cannot handle values >= 1.
#' @return The transformed vector.
#' @export
#' @importFrom base asinh
#'
arcsinh <- function(vec, scaling.factor=100, norm.by.scaling=TRUE){
    scaled <- vec * scaling.factor

    if(norm.by.scaling){
        return( asinh(scaled)/asinh(scaling.factor) )
    } else {
        return( asinh(scaled) )
    }
}
