#' Determines the expected power in a future study with specified df for the future study
#' based on the uncertainty associated with an existing study.
#' Uses the F-test of the previous study as input in the sample size planning process.
#'
#' @param Ftest The F-test of the previous study.
#' @param df1 The numerator degrees of freedom of the previous study. This must be the same in the previous and new study.
#' @param df2 The denominator degrees of freedom of the previous study.
#' @param df2new The denominator degrees of freedom for the new study.
#' @param alpha The significance level. Default is \eqn{\alpha = .05}.
#' @param filter The filter value reflects the probability of nonsignificant results being filtered.
#' filter = 0 means that there is no filtering and you would have observed nonsignificant results.
#' filter = 1 means that only significant results are observed and you would never have seen nonsigificant results if they had occurred.
#' Filtering is based on alpha = .05 and assumes that are have observed a significant result.
#' Filtering is conducted by weighting (actually filtering) the posterior distribution.
#' For instance, if filter = 1, then the posterior of the null (i.e., the noncentrality parameter is 0)
#' is up to 20 times more likely than when the noncentrality parameter is very large. Setting filter > 0 slows estimation.
#' @param upper_null Specifies the upper value of the composite null hypothesis in units of Cohen's f.
#' The default value of upper_null = 0 keeps the point null hypothesis.
#' A value of, for instance, upper_null = .05 would remove all posterior values between -.05 and .05 from consideration when calculating expected power.
#' @param estimate_fixed Specifies whether the predictor in the regression model is either fixed or random. The default is FALSE for random predictors.
#' @param future_fixed Specifies whether the future study will have fixed predictors.
#' @return Returns expected power for the prospective study.
#' @export
#' @examples
#' \dontrun{
#' ep_Ftest(Ftest = 5.0, df1=2, df2=50, df2new=100, alpha = .05)
#' }
ep_Ftest <- function(Ftest, df1, df2, df2new, alpha=0.05, filter=0, upper_null=0, estimate_fixed=TRUE, future_fixed=TRUE){
	# upper_null is in units of Cohen's f.

	exit <- 0
	if (Ftest < 0){
		cat("The F-test must be > 0.\n")
		exit<-1
	}
	if (df1 < 1) {
		cat("df1 must be greater than 1.\n")
		exit<-1
		}
	if (df2 < 1) {
		cat("df2 must be greater than 1.\n")
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

	  posterior <- posterior_Cohen_f(Ftest=Ftest, df1=df1, df2=df2, filter=filter,
	                                 upper_null=upper_null, fixed=estimate_fixed, ndraws=5000000)
	  ndraws<- length(posterior)

		N_new <- df2new + df1 + 1

		z <- rnorm(n= ndraws,mean=0,sd=1)
		c1 <- sqrt(rchisq(n= ndraws,df= df1-1,ncp=0)) #This is indeed df1 - 1 and is a constant 0 when df1=1.
		c2 <- sqrt(rchisq(n= ndraws,df= df2new,ncp=0))
		cN <- sqrt(rchisq(n= ndraws,df= N_new-1,ncp=0))

		#Posterior Predictive Distribution
		if (future_fixed == TRUE){
			ppF <- (df2new/df1)*((z+sqrt(N_new)*posterior)^2 + c1^2)/c2^2
		} else {
			ppF <- (df2new/df1)*((z+ cN*posterior)^2 + c1^2)/c2^2
		}

		dist <- ecdf(ppF)
		EP <- 1 - dist(qf(alpha, df1=df1, df2= df2new, ncp=0, lower.tail=FALSE))

	  return(EP)
	}
}


