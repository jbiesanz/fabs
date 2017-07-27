#' Samples from the posterior distribution of Cohen's f. Requires the F-test as input
#' and calculates Cohen's \eqn{f} internally. The posterior distribution is generated using
#' Fisher's fiducial approach.
#'
#' @param Ftest The F-test of the previous study.
#' @param df1 The numerator degrees of freedom.
#' @param df2 The denominator degrees of freedom.
#' @param filter The filter value reflects the probability of nonsignificant results being filtered.
#' filter = 0 means that there is no filtering and you would have observed nonsignificant results.
#' filter = 1 means that only significant results are observed and you would never have seen nonsigificant results if they had occurred.
#' Filtering is based on alpha = .05 and assumes that are have observed a significant result.
#' Filtering is conducted by weighting (actually filtering) the posterior distribution.
#' For instance, if filter = 1, then the posterior of the null (i.e., the noncentrality parameter is 0)
#' is up to 20 times more likely than when the noncentrality parameter is very large. Setting filter > 0 slows estimation.
#' @param upper_null Specifies the upper value of the composite null hypothesis in units of Cohen's f.
#' The default value of upper_null = 0 keeps the point null hypothesis.
#' A value of, for instance, upper_null = .05 would remove all posterior values between -.05 and .05.
#' @param fixed Specifies whether the predictors in the model are either fixed or random. The default is FALSE for random predictors.
#' @param ndraws Specifies the number of initial samples from the posterior. For small effect sizes and when filter or upper_null > 0
#' the number of returned samples from the posterior distribution will be lower than ndraws.
#' @return Returns a vector of samples from the posterior distribution in units of Cohen's \eqn{f}.
#' @export
#' @examples
#' posterior_Cohen_f(Ftest = 5.0, df1=2, df2=50, ndraws = 10)
posterior_Cohen_f <- function(Ftest, df1, df2, filter=0, upper_null=0, fixed=FALSE, ndraws = 200000){
	#upper_null is in units of Cohen's f
	f <- sqrt((df1/df2)*Ftest)
	options(warn=-1)  #setting warning option to silent since NAs are produced when we generate the posterior of f. These can be ignored as they refer to negative variance estimates (and are set to 0)


	z <- rnorm(n= ndraws,mean=0,sd=1)
	c1 <- sqrt(rchisq(n= ndraws,df=df1-1,ncp=0))  #This is indeed df1 - 1 and is a constant 0 when df1=1.
	c2 <- sqrt(rchisq(n= ndraws,df=df2,ncp=0))
	cN <- sqrt(rchisq(n= ndraws,df=df2+df1,ncp=0))

	if (fixed == FALSE){
		posterior_f <- sqrt((sqrt(c2^2*f^2 - c1^2)+ z)^2/cN^2)
	} else if (fixed == TRUE){
		posterior_f <- sqrt((sqrt(c2^2*f^2 - c1^2)+ z)^2/(df2+1))
	} else cat("Fixed refers to whether the predictor is an observed variable and must be either TRUE or FALSE.\n")

	#Cleaning up negative variance estimates. These produce either NAs or negative values.
	#Setting these negative posterior parameter estimates to 0.
	posterior_f[is.na(posterior_f)] <- 0
	posterior_f[posterior_f < 0] <- 0
	posterior_f[posterior_f > upper_null]
	ndraws <- length(posterior_f)

	#Weighting using the publication bias filter
	if (filter > 0){
		F.crit <- qf(.95,df1=df1, df2=df2,ncp=0)
		posterior_lambda <- (df2+df1+1)*posterior_f^2
		lower.tail.proportion <- pf(F.crit,df1=df1, df2=df2, ncp= posterior_lambda, lower.tail=TRUE)

		W <- 1/(lower.tail.proportion*(1-filter) + (1-lower.tail.proportion))
		select <- runif(ndraws, min=0, max=1)
		posterior_f <- posterior_f[select < W/max(W)]
	}
	posterior_f
}

#' Samples from the posterior distribution of the correlation coefficient. The posterior distribution is
#' generated using Fisher's fiducial approach and corresponds exactly to the use of a Jeffrey's prior.
#'
#' @param r The correlation of the previous study.
#' @param df The degrees of freedom in the previous study.
#' @param filter The filter value reflects the probability of nonsignificant results being filtered.
#' filter = 0 means that there is no filtering and you would have observed nonsignificant results.
#' filter = 1 means that only significant results are observed and you would never have seen nonsigificant results if they had occurred.
#' Filtering is based on alpha = .05 and assumes that are have observed a significant result.
#' Filtering is conducted by weighting (actually filtering) the posterior distribution.
#' For instance, if filter = 1, then the posterior of the null (i.e., the noncentrality parameter is 0)
#' is up to 20 times more likely than when the noncentrality parameter is very large. Setting filter > 0 slows estimation.
#' @param upper_null Specifies the upper value of the composite null hypothesis in units of the correlation.
#' The default value of upper_null = 0 keeps the point null hypothesis.
#' A value of, for instance, upper_null = .05 would remove all posterior values between -.05 and .05.
#' @param fixed Specifies whether the predictors in the model are either fixed or random. The default is FALSE for random predictors.
#' @param ndraws Specifies the number of initial samples from the posterior. For small effect sizes and when filter or upper_null > 0
#' the number of returned samples from the posterior distribution will be lower than ndraws.
#' @return A vector of samples from the posterior distribution in units of the correlation parameter \eqn{\rho}.
#' @export
posterior_r <- function(r, df, filter=0, upper_null=0, fixed=FALSE, ndraws = 200000){
	#upper_null is in units of the correlation
	tobs <- r*sqrt(df)/sqrt(1 - r^2)
	N <- df + 2

	z <- rnorm(n=ndraws,mean=0,sd=1)
	cdf <-  sqrt(rchisq(n= ndraws,df= df, ncp=0))
	cdf1 <-  sqrt(rchisq(n= ndraws,df= df + 1,ncp=0))

	if (abs(r) > 1) cat("The correlation (r) must be between -1 and 1.\n")
	if (df < 3) cat("Degrees of freedom (df) must be greater than 3.\n")
	if (filter < 0) cat("filter must be between or equal to 0 and 1.\n")
	if (filter > 1) cat("filter must be between or equal to 0 and 1.\n")
	if (upper_null < 0) cat("The upper_null value must be between 0 and 1 (and closer to 0).\n")
	if (upper_null > 1) cat("The upper_null value must be between 0 and 1 (and closer to 0).\n")

	if (fixed == FALSE){
		posteriorL <- (z*sqrt(df) + tobs*cdf)/cdf1
		posterior <- posteriorL/sqrt(posteriorL*posteriorL + N)
	} else if (fixed == TRUE){
		posteriorL <- z + tobs*cdf/sqrt(df)
		posterior <- posteriorL/sqrt(posteriorL*posteriorL + N)
	} else 	cat("Fixed refers to whether the predictor is an observed variable and must be either TRUE or FALSE.\n")

	#applying the composite null filter
	posterior <- posterior[abs(posterior) > upper_null]

	#updating ndraws to reflect the current length of the posterior
	ndraws <- length(posterior)

	#converting rho back to lambda
	posterior_lambda <- posterior*sqrt(N)/sqrt(1 - posterior^2)

	#Weighting using the publication bias filter
	if (filter > 0){
		t.lower.crit <- qt(.025,df=df,ncp=0)
		t.upper.crit <- qt(.975,df=df,ncp=0)

		lower.tail.proportion <- pt(t.lower.crit, df=df, ncp= posterior_lambda, lower.tail=TRUE)
		upper.tail.proportion <- pt(t.upper.crit, df=df, ncp= posterior_lambda, lower.tail=FALSE)
		middle.proportion <- 1 - (lower.tail.proportion + upper.tail.proportion)

		W <- 1/(lower.tail.proportion + upper.tail.proportion + (1-filter)* middle.proportion)
		select <- runif(ndraws, min=0, max=1)
		posterior <- posterior[select < W/max(W)]
	}
	posterior
}

#' Samples from the posterior distribution of the standardized mean difference.
#' The posterior distribution is generated using Fisher's fiducial approach and
#' corresponds exactly to the use of a Jeffrey's prior. Assumes fixed predictor.
#'
#' @param d The observed standardized mean difference of the previous study based on
#' a pooled standard deviation.
#' @param n1 The number of observations in the first group.
#' @param n2 The number of observations in the second group.
#' @param filter The filter value reflects the probability of nonsignificant results being filtered.
#' filter = 0 means that there is no filtering and you would have observed nonsignificant results.
#' filter = 1 means that only significant results are observed and you would never have seen nonsigificant results if they had occurred.
#' Filtering is based on alpha = .05 and assumes that are have observed a significant result.
#' Filtering is conducted by weighting (actually filtering) the posterior distribution.
#' For instance, if filter = 1, then the posterior of the null (i.e., the noncentrality parameter is 0)
#' is up to 20 times more likely than when the noncentrality parameter is very large. Setting filter > 0 slows estimation.
#' @param upper_null Specifies the upper value of the composite null hypothesis in units of Cohen's f.
#' The default value of upper_null = 0 keeps the point null hypothesis.
#' A value of, for instance, upper_null = .05 would remove all posterior values between -.05 and .05.
#' @param ndraws Specifies the number of initial samples from the posterior. For small effect sizes and when filter or upper_null > 0
#' the number of returned samples from the posterior distribution will be lower than ndraws.
#' @return A vector of samples from the posterior distribution in units of Cohen's \eqn{f}.
#' @export
posterior_d <- function(d, n1, n2, filter=0, upper_null=0, 	ndraws = 200000){
	#upper_null is in units of Cohen's d -- the standardized mean difference
	df <- n1 + n2 - 2
	tobs <- d*sqrt((n1*n2)/(n1+n2))

	z <- rnorm(n=ndraws,mean=0,sd=1)
	cdf <-  sqrt(rchisq(n= ndraws,df= df, ncp=0))
	cdf1 <-  sqrt(rchisq(n= ndraws,df= df + 1,ncp=0))

	posterior <- z*sqrt(1/n1 + 1/n2) +d*cdf/sqrt(df)

	#applying the composite null filter
	posterior <- posterior[abs(posterior) > upper_null]

	#updating ndraws to reflect the current length of the posterior
	ndraws <- length(posterior)

	#converting delta back to lambda
	posterior_lambda <- posterior/sqrt(1/n1 + 1/n2)

	#Weighting using the publication bias filter
	if (filter > 0){
		t.lower.crit <- qt(.025,df=df,ncp=0)
		t.upper.crit <- qt(.975,df=df,ncp=0)

		lower.tail.proportion <- pt(t.lower.crit, df=df, ncp= posterior_lambda, lower.tail=TRUE)
		upper.tail.proportion <- pt(t.upper.crit, df=df, ncp= posterior_lambda, lower.tail=FALSE)
		middle.proportion <- 1 - (lower.tail.proportion + upper.tail.proportion)

		W <- 1/(lower.tail.proportion + upper.tail.proportion + (1-filter)* middle.proportion)
		select <- runif(ndraws, min=0, max=1)
		posterior <- posterior[select < W/max(W)]
	}
	posterior
}

#' Samples from the posterior distribution of the noncentral t-distribution.
#' The posterior distribution is generated using Fisher's fiducial approach.
#'
#' @param t The t-test of the previous study.
#' @param df The degrees of freedom associated with the t-test.
#' @param fixed Is the predictor fixed or random? Default is FALSE to indicate randomly sampled predictor.
#' @param ndraws Specifies the number of initial samples from the posterior. For small effect sizes and when filter or upper_null > 0
#' the number of returned samples from the posterior distribution will be lower than ndraws.
#' @return A vector of samples from the posterior distribution of the noncentral t-distribution.
#' @export
posterior_t <- function(t, df, fixed=FALSE, ndraws = 200000){
	ndraws<- 200000
	tobs <- t
	if (df < 1) cat("Degrees of freedom (df) must be greater than 1.\n")
	if (fixed == FALSE){
		posterior <- (rnorm(n= ndraws,mean=0,sd=1)*sqrt(df) + tobs*sqrt(rchisq(n= ndraws,df=df,ncp=0)))/sqrt(rchisq(n= ndraws,df=df+1,ncp=0))
	} else if (fixed == TRUE){
		posterior <- rnorm(n= ndraws,mean=0,sd=1) + tobs*sqrt(rchisq(n= ndraws,df=df,ncp=0)/(df))
	} else 	cat("Fixed refers to whether the predictor is an observed variable and must be either TRUE or FALSE.\n")
	posterior
}
