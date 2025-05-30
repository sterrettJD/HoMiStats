% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filters.R
\name{filter_low_abundance_by_mean}
\alias{filter_low_abundance_by_mean}
\title{Filter low abundance features from relative abundance data}
\usage{
filter_low_abundance_by_mean(df, threshold)
}
\arguments{
\item{df}{A dataframe containing relative abundance data,
where rows are samples and columns are features (genes).}

\item{threshold}{The threshold for filtering low abundance features.}
}
\value{
A dataframe with low abundance features filtered out.
}
\description{
Filters low abundance features with
means below a certain threshold.
}
\examples{
df <- data.frame(
  feature1 = c(0.1, 0.2, 0.3, 0.4),
  feature2 = c(0.2, 0.3, 0.4, 0.5),
  feature3 = c(0.3, 0.4, 0.5, 0.6)
)
filtered_df <- filter_low_abundance_by_mean(df, 0.5)
print(filtered_df)
}
