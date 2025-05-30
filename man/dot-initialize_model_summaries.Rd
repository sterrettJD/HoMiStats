% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HoMiCorr.R
\name{.initialize_model_summaries}
\alias{.initialize_model_summaries}
\title{Initialize Data Frame for Model Summaries (internal)}
\usage{
.initialize_model_summaries(feature.combos, reg.method)
}
\arguments{
\item{feature.combos}{A list of feature combinations to be analyzed.}

\item{reg.method}{A string specifying the regression method to be used.
Supported options: "gamlss", "zibr", "lm", "lmer".}
}
\value{
A data frame initialized with appropriate columns for the chosen
regression method, with the number of rows equal to the number of feature
combinations.
}
\description{
Creates an empty data frame to store regression model summaries.
The columns included in the data frame depend on the selected
regression method.
}
\examples{
# Simulated feature combinations
feature_combos <- list(c("a1", "b2"), c("b1", "a2"))

# Initialize model summaries for linear regression
summaries <- HoMiStats:::.initialize_model_summaries(feature_combos, "lm")

}
\keyword{internal⁠}
