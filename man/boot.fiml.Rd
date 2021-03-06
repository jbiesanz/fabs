% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/missing_data_lm_functions.R
\name{boot.fiml}
\alias{boot.fiml}
\title{Bootstraps the full information maximum likelihood regression model under lavaan.}
\usage{
boot.fiml(data, model, z.function = FALSE, R = 9999, conf = 0.95)
}
\arguments{
\item{data}{The dataset for the analysis.}

\item{model}{A regression model from lm()}

\item{z.function}{A user supplied function to standardize variables and create nonlinear terms. Called internally.}

\item{R}{The number of resamples to be conducted.}

\item{conf}{The confidence level for intervals.}
}
\value{
Returns (1) maximum likelihood regression coefficient estimates
(s) percentile and BCa confidence intervals using case-wise resampling
}
\description{
Bootstraps the full information maximum likelihood regression model under lavaan.
}
\examples{
\dontrun{
x <- c(1,2,3,4,5,6)
y <- c(2.1,1.9,2.1,1.8,10,2.2)
temp_data <- as.data.frame(cbind(y, x))
lm_model <- lm(y~x, data=temp_data)
robust_lm_inferences(data=temp_data, model=lm_model)
}
}
