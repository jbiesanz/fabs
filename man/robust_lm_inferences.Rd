% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/robust_lm_inferences.R
\name{robust_lm_inferences}
\alias{robust_lm_inferences}
\title{Estimates a series of robust regressions for a specified regression model.}
\usage{
robust_lm_inferences(data, model, R = 9999, conf = 0.95)
}
\arguments{
\item{data}{The dataset for the analysis.}

\item{model}{A regression model from lm()}

\item{R}{The number of resamples to be conducted.}

\item{conf}{The confidence level for intervals.}
}
\value{
Returns (1) the lm() regression model with HC4 heteroscedastic covariance consistent standard errors
(2) p-values based on the wild bootstrap (and HC4 standard errors) under the null hypothesis.
(3) percentile and BCa confidence intervals using case-wise resampling
(4) Robust regression analysis from rlm() under the default Huber loss function.
(5) Robust regression with resampled wild bootstrap standard errors.
}
\description{
Estimates a series of robust regressions for a specified regression model.
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
