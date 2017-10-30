#' Determines the expected power a future study with a specified sample size
#' based on the uncertainty associated with an existing study.
#' Uses the estimated correlation of the previous study as input in the sample size planning process.
#'
#' @param r The correlation from the previous study.
#' @param df The degrees of freedom associated with the previous study
#' @param dfnew The degrees of freedom associated with the prospective study
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
#' A value of, for instance, upper_null = .05 would remove all posterior values between -.05 and .05 from consideration when calculating expected power.
#' @param estimate_fixed Specifies whether the predictor in the regression model is either fixed or random. The default is FALSE for random predictors.
#' @param future_fixed Specifies whether the future study will have fixed predictors.
#' @return Returns (1) the sample size required for the future study to achieve the specified level of expected power.
#' This reflects the uncertainty associated with the previous study and (2) the median 95% confidence interval width for the correlation in the prospective study.
#' @export
#' @examples
#' \dontrun{
#' ep_r(r = .59, df=13, dfnew=50)
#'
#' ep_r(r = .59, df=13, dfnew=50, filter = 1)
#' }
ep_r <- function(r, df, dfnew, alpha=0.05, filter=0, upper_null=0, estimate_fixed=FALSE, future_fixed=FALSE){
	exit <- 0
	if (abs(r) > 1){
		cat("The correlation (r) must be between -1 and 1.\n")
		exit<-1
	}
	if (df < 3) {
		cat("Degrees of freedom must be greater than 3.\n")
		exit<-1
		}
	if (filter < 0) {cat("filter must be between or equal to 0 and 1.\n")
		exit<-1
		}
	if (filter > 1) {cat("filter must be between or equal to 0 and 1.\n")
		exit<-1
		}
	if (upper_null < 0) {cat("The upper_null value must be greater than 0 (and close to 0).\n")
		exit<-1
		}
	if (upper_null > 1) {cat("The upper_null value must be greater than 0 (and close to 0). Excluding extremely large effect sizes does not make sense.\n")
		exit<-1
		}
	if (alpha > 1) {cat("alpha must be between 0 and 1. Default value is 0.05 if not specified.\n")
		exit<-1
		}
	if (alpha < 0) {cat("alpha must be between 0 and 1. Default value is 0.05 if not specified.\n")
		exit<-1
		}
	if (alpha == 0) {cat("alpha must be between 0 and 1. Default value is 0.05 if not specified.\n")
		exit<-1
		}

	if (exit == 0){
  	#Initial model computing expected power across the starting values ranging from n=3 to 50,000
  	posterior <- posterior_r(r=r, df=df, filter=filter, upper_null=upper_null,
  	                         fixed=estimate_fixed, ndraws=5000000)  #Single static posterior distribution with large number of draws
  	ndraws<- length(posterior)

		n_new <- dfnew + 2

		z <- rnorm(n=ndraws,mean=0,sd=1)
		cdf <-  sqrt(rchisq(n= ndraws,df= dfnew,ncp=0))
		cdf1 <-  sqrt(rchisq(n= ndraws,df= dfnew + 1,ncp=0))

		#Posterior Predictive Distribution
		if (future_fixed == TRUE){
			posterior_lambda <- posterior*sqrt(n_new/(1-posterior^2))
			ppt <-(z + posterior_lambda)/(cdf/sqrt(dfnew))
		} else {
			posterior_lambda <- cdf1*posterior/(1-posterior^2)
			ppt <-(z + posterior_lambda)/(cdf/sqrt(dfnew))
		}

		abs_ppt <- abs(ppt)
		dist <- ecdf(abs_ppt)
		EP <- 1 - dist(qt(alpha/2,df=dfnew,ncp=0,lower.tail=FALSE))

	  return(EP)
	}
}


