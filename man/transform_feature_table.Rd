% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transformers.R
\name{transform_feature_table}
\alias{transform_feature_table}
\title{Transform feature table for analysis}
\usage{
transform_feature_table(feature.table, transformation)
}
\arguments{
\item{feature.table}{A dataframe where rows are samples
and columns are genes/features. Row names should be sample IDs.}

\item{transformation}{Method for transforming the feature.table.
String options include "arcsinh_<scaling factor>_norm" and
"arcsinh_<scaling factor>_nonorm".
The arcsinh transformation applies a linear transformation
for small values and a log transformation for larger values.
The scaling factor controls how much of the smaller values
are transformed linearly.
"norm" arcsinh transformations normalize values such that
the maximum after transformation is 1.
Additionally, a function can be provided for custom transformations.}
}
\value{
The transformed table.
}
\description{
Transforms feature table for differential expression or
co-occurrence analysis.
}
\examples{
feature.table <- data.frame(gene1 = c(0.1, 0.5, 0.2),
                            gene2 = c(0.3, 0.7, 0.1))
rownames(feature.table) <- c("sample1", "sample2", "sample3")
transform_feature_table(feature.table, "arcsinh_100_norm")
transform_feature_table(feature.table, function(x) log(x + 1))

}
