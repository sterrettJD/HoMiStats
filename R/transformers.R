#' Transform feature table for analysis
#' @description Transforms feature table for differential expression or
#' co-occurrence analysis.
#' @param feature.table A dataframe where rows are samples
#' and columns are genes/features. Row names should be sample IDs.
#' @param transformation Method for transforming the feature.table.
#' String options include "arcsinh_<scaling factor>_norm" and
#' "arcsinh_<scaling factor>_nonorm".
#' The arcsinh transformation applies a linear transformation
#' for small values and a log transformation for larger values.
#' The scaling factor controls how much of the smaller values
#' are transformed linearly.
#' "norm" arcsinh transformations normalize values such that
#' the maximum after transformation is 1.
#' Additionally, a function can be provided for custom transformations.
#' @return The transformed table.
#' @export
#' @importFrom dplyr mutate_all
#' @importFrom base strsplit
#' @examples
#' feature.table <- data.frame(gene1 = c(0.1, 0.5, 0.2),
#'                             gene2 = c(0.3, 0.7, 0.1))
#' rownames(feature.table) <- c("sample1", "sample2", "sample3")
#' transform_feature_table(feature.table, "arcsinh_100_norm")
#' transform_feature_table(feature.table, function(x) log(x + 1))
#'
transform_feature_table <- function(feature.table, transformation){
    check_proportional(feature.table, soft=TRUE)

    # SUPPORTED CHARACTER TRANSFORMATIONS
    if(mode(transformation)=="character"){
        # Arcsinh transform with renormalization (max return value is 1)
        # Scaling factor is passed between underscores
        if(grepl(pattern="^arcsinh_[0-9]+_norm$", x=transformation)){
            arcsinh.args <- unlist(strsplit(x=transformation, split="_"))
            scaling.factor <- as.numeric(arcsinh.args[2])

            return(dplyr::mutate_all(feature.table,
                                     function(x) arcsinh(x,
                                                         scaling.factor,
                                                         norm.by.scaling=TRUE)))
        }
        # Arcsinh transform with NO normalization (max return value can be > 1)
        # Scaling factor is passed between underscores
        if(grepl(pattern="^arcsinh_[0-9]+_nonorm$", x=transformation)){
            arcsinh.args <- unlist(strsplit(x=transformation, split="_"))
            scaling.factor <- as.numeric(arcsinh.args[2])

            return(dplyr::mutate_all(feature.table,
                                     function(x) arcsinh(x,
                                                         scaling.factor,
                                                         norm.by.scaling=FALSE)
                                                         )
                                    )
        }

        stop("The requested transformation is not yet supported. ",
             "You are welcome to implement it via passing a function ",
             "for the transformation argument.")
    }

    # If not a character, try to implement the function provided
    dplyr::mutate_all(feature.table, transformation)
}

#' Check if feature table rows sum to 1
#' @description Some transformations require proportional data
#' where sample rows sum to 1. This function checks that condition.
#' @param feature.table A dataframe where rows are samples
#' and columns are features.
#' @param tolerance Numerical value for tolerance in checking if rows sum to 1.
#' @param soft Boolean.
#' If TRUE, the function only warns instead of stopping execution.
#' @return Nothing; throws a warning or error if the check fails.
#' @export
#' @examples
#' feature.table <- data.frame(gene1 = c(0.4, 0.5, 0.6),
#'                             gene2 = c(0.6, 0.5, 0.4))
#' rownames(feature.table) <- c("sample1", "sample2", "sample3")
#' check_proportional(feature.table, soft=TRUE)
#'
check_proportional <- function(feature.table, tolerance=1e-3, soft=FALSE){
    sample.sums <- rowSums(feature.table)
    rows.not.prop <- which((sample.sums < (1-tolerance)) |
                           (sample.sums > (1+tolerance)))
    rownames.not.prop <- rownames(feature.table)[rows.not.prop]

    if(length(rownames.not.prop) > 0){
        rownames.str <- paste(rownames.not.prop, collapse=", ")
        if(soft){
            warning("The following samples do not sum to 1: ",
                    rownames.str,
                    "\nPlease make sure this is intentional...")
        } else {
            stop("The following samples do not sum to 1: ",
                 rownames.str)
        }

    }
}

#' Apply arcsinh transformation to a numeric vector
#' @description Performs an arcsinh transformation with scaling.
#' The arcsinh function behaves like a log transformation for large values
#' but remains linear for small values and can handle zeros.
#' @param vec A numeric vector to be transformed.
#' @param scaling.factor A scaling factor applied before transformation.
#' Increasing this factor increases how much of the lower end of the data
#' is treated linearly.
#' @param norm.by.scaling Boolean indicating whether transformed data
#' should be normalized by the arcsinh of the scaling factor.
#' @return A transformed numeric vector.
#' @export
#' @importFrom base asinh
#' @examples
#' vec <- c(0.01, 0.1, 1, 10, 100)
#' arcsinh(vec, scaling.factor=100, norm.by.scaling=TRUE)
#' arcsinh(vec, scaling.factor=50, norm.by.scaling=FALSE)
#'
arcsinh <- function(vec, scaling.factor=100, norm.by.scaling=TRUE){
    scaled <- vec * scaling.factor

    if(norm.by.scaling){
        return( asinh(scaled)/asinh(scaling.factor) )
    } else {
        return( asinh(scaled) )
    }
}
