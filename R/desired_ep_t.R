#' Determines the required sample size in a future study to achieve the desired expected power (ep)
#' based on the uncertainty associated with an existing study.
#' Uses the t-test of the previous study as input in the sample size planning process. Converts
#' the t-test to the correlational metric, determines sample size for the prospective study,
#' and provides the median confidence interval width for the correlation in the prospective study
#' based on that sample size.
#'
#' @param t The t-test of the previous study.
#' @param df The degrees of freedom associated with the previous study
#' @param desired_ep The desired expected power for the future study.
#' @param alpha The signficance level. Default is .05.
#' @param filter The filter value reflects the probability of nonsignificant results being filtered.
#' filter = 0 means that there is no filtering and you would have observed nonsignificant results.
#' filter = 1 means that only significant results are observed and you would never have seen nonsigificant results if they had occurred.
#' Filtering is based on alpha = .05 and assumes that are have observed a significant result.
#' Filtering is conducted by weighting (actually filtering) the posterior distribution.
#' For instance, if filter = 1, then the posterior of the null (i.e., the noncentrality parameter is 0)
#' is up to 20 times more likely than when the noncentrality parameter is very large. Setting filter > 0 slows estimation.
#' @param upper_null Specifies the upper value of the composite null hypothesis in units of the correlation coefficient eqn{r}.
#' The default value of upper_null = 0 keeps the point null hypothesis.
#' A value of, for instance, upper_null = .05 would remove all posterior values between -.05 and .05
#' from consideration when calculating expected power.
#' @param estimate_fixed Specifies whether the predictor in the regression model is either fixed or random.
#' The default is FALSE for random predictors.
#' @param future_fixed Specifies whether the future study will have fixed predictors.
#' @return Returns (1) the sample size required for the future study to achieve the specified level of expected power.
#' This reflects the uncertainty associated with the previous study and (2) the median 95% confidence interval width
#' for the correlation in the prospective study.
#' @export
#' @importFrom graphics axis par plot segments text
#' @examples
#' \dontrun{
#' desired_ep_t(t = 3.0, df=40, desired_ep = 0.80)
#' }
desired_ep_t <- function(t, df, desired_ep = 0.80, alpha=0.05, filter=0, upper_null=0, estimate_fixed=FALSE, future_fixed=FALSE){
	#convert t to r and call the desired_ep_r function.
	r <- t/sqrt(t^2 + df)
	desired_ep_r(r=r, df=df, desired_ep = desired_ep, alpha=alpha, filter=filter, upper_null= upper_null, estimate_fixed = estimate_fixed, future_fixed = future_fixed)
}
