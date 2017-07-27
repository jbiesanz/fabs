#' Creates a confidence interval for the squared multiple correlation.
#'
#' @param r2 The estimated squared multiple correlation.
#' @param df1 The numerator degrees of freedom.
#' @param df2 The denominator degrees of freedom.
#' @param conf The desired (1-\eqn{\alpha}) confidence level. Default is .95.
#' @param fixed Specifies whether the predictors in the regression model are either fixed or random. The default is random predictors.
#' @return Returns the lower and upper endpoints of the confidence interval.
#' @export
#' @import stats
#' @import graphics
ci_r2 <- function(r2, df1, df2, conf=.95, fixed =FALSE){
  options(warn=-1)  #setting warning option to silent since NAs are produced.
  M <-5000000
  alphaL <- (1-conf)/2
  alphaU <- 1- alphaL
  predictor_warning <- 0

  if (fixed == TRUE){
    Cp1 <- rchisq(n= M, df=df1-1, ncp=0)
    Cnp1<- rchisq(n= M, df=df2, ncp=0)
    z <- rnorm(n=M, mean=0, sd=1)
    samplesize <- df1 + df2 + 1

    B <- Cnp1*r2/(1-r2) - Cp1
    D <- (sqrt(B) + z)^2/samplesize
    D[is.na(D)] <- 0
    CID <- D/(1 + D)
    conf_interval <- quantile(CID,probs=c(alphaL, alphaU))

  } else if (fixed == FALSE) {
    Cp1 <- rchisq(n= M, df=df1-1, ncp=0)
    Cnp1<- rchisq(n= M, df=df2, ncp=0)
    Cn1<- rchisq(n=M , df=df1+df2, ncp=0)
    z <- rnorm(n=M, mean=0, sd=1)

    B <- Cnp1*r2/(1-r2) - Cp1
    D <- (sqrt(B) + z)^2/(Cn1)
    D[is.na(D)] <- 0
    CID <- D/(1 + D)

    conf_interval <-  quantile(CID,probs=c(alphaL, alphaU))
  } else {
    predictor_warning <- 1
    cat("Please specify fixed_predictors as TRUE or FALSE")
  }
  options(warn=0)   #resetting warning back to original state
  if (predictor_warning == 0) conf_interval
}

#' Creates a confidence interval for Cohen's f.
#'
#' @param f Estimated value of Cohen's f.
#' @param df1 The degrees of freedom numerator for the F-test associated with f.
#' @param df2 The degrees of freedom denominator for the F-test associated with f.
#' @param conf The desired (1-\eqn{\alpha}) confidence level. Default is .95.
#' @param fixed Specifies whether the predictors in the model are either fixed or random. The default is FALSE for random predictors.
#' @return Returns the lower and upper endpoints of the confidence interval.
#' @export
ci_f <- function(f, df1, df2, fixed = FALSE, conf = 0.95){
	options(warn=-1)  #setting warning option to silent since NAs are produced.
	M <-2000000
	alphaL <- (1-conf)/2
	alphaU <- 1- alphaL

	z <- rnorm(n= M,mean=0,sd=1)
	c1 <- sqrt(rchisq(n= M,df=df1-1,ncp=0))  #This is indeed df1 - 1 and is a constant 0 when df1=1.
	c2 <- sqrt(rchisq(n= M,df=df2,ncp=0))
	cN <- sqrt(rchisq(n= M,df=df1+df2-1,ncp=0))

	if (fixed == TRUE){
    	CID <- (sqrt(f^2*c2^2 - c1^2) + z)/sqrt(df2 + df1 + 1)
	} else {
		CID <- (sqrt(f^2*c2^2 - c1^2) + z)/cN
	}

 	CID[is.na(CID)] <- 0
	CID[CID < 0] <- 0
 	options(warn=0)   #resetting warning back to original state

	initial <- quantile(CID,probs=c(alphaL, alphaU))
	c(initial[1], initial[2])
}

#' Creates a confidence interval for the standardized mean difference. Assumes predictors are fixed.
#'
#' @param d The estimated standardized mean difference.
#' @param n1 The sample size for group 1.
#' @param n2 The sample size for group 2.
#' @param conf The desired (1-\eqn{\alpha}) confidence level. Default is .95.
#' @param iter Number of iterations for the stochastic approximation optimization. Default is 50.
#' @return Returns the lower and upper endpoints of the confidence interval.
#' @export
ci_d <- function(d, n1, n2, conf=.95, iter=50){
  M <-200000
  alphaL <- (1-conf)/2
  alphaU <- 1- alphaL
  df <- n1+n2-2
  tobs<- d/sqrt(1/n1 + 1/n2)

  CID <- rnorm(M,mean=0,sd=1)*sqrt(1/n1 + 1/n2) + sqrt(1/n1 + 1/n2)*tobs*sqrt(rchisq(n= M,df=df,ncp=0)/(df))
  initial <- quantile(CID,probs=c(alphaL, alphaU))
  initialqL <- initial[1]
  initialqU <- initial[2]
  lastqL <- initialqL
  lastqU <- initialqU

  lastdenL <- density(CID, from= lastqL, to = lastqL +1)$y[1]
  initialdenL <- lastdenL
  lastdenU <- density(CID, from= lastqU, to = lastqU +1)$y[1]
  initialdenU <- lastdenU

  for (n in 1: iter){
    CID <- rnorm(M,mean=0,sd=1)*sqrt(1/n1 + 1/n2) +sqrt(1/n1 + 1/n2)*tobs*sqrt(rchisq(n= M,df=df,ncp=0)/(df))

    fnL <- (1-1/n)*lastdenL + (1/n)*(sum(abs(CID-lastqL)<=(1/sqrt(n))))/(2*M*(1/sqrt(n)))
    SnL <- lastqL + (1/(n*max(fnL,initialdenL/sqrt(n))))*(alphaL - sum(CID <= lastqL)/M)

    fnU <- (1-1/n)*lastdenU + (1/n)*(sum(abs(CID-lastqU)<=(1/sqrt(n))))/(2*M*(1/sqrt(n)))
    SnU <- lastqU + (1/(n*max(fnU,initialdenU/sqrt(n))))*(alphaU - sum(CID <= lastqU)/M)
    lastqL <- SnL
    lastdenL <- fnL
    lastqU <- SnU
    lastdenU <- fnU
  }
  c(SnL, SnU)
}

#' Creates a confidence interval for the standardized regression coefficient \eqn{\beta}.
#'
#' @param t The t-test for the standardized regression coefficient.
#' @param df1 The numerator degrees of freedom for the regression model.
#' @param df2 The denominator degrees of freedom for the regression model.
#' @param r2 The squared multiple correlation for the full regression model.
#' @param r2x The coefficient of determination for the predictor of interest. This is the squared multiple correlation when the predictor of interest is predicted by all of the other variables in the regression model.  1 - tolerance.
#' @param fixed Specifies whether the predictors in the regression model are either fixed or random. The default is random predictors.
#' @param conf The desired (1-\eqn{\alpha}) confidence level. The default is set to .95.
#' @param iter Number of iterations for the stochastic approximation optimization. Default is 50.
#' @return Returns the lower and upper endpoints of the confidence interval.
#' @export
ci_beta <- function(t, df1, df2, r2, r2x, conf=.95, fixed=FALSE, iter=50){
  options(warn=-1)  #setting warning option to silent since NAs are produced when we take the square root of B and this generates warnings that are not needed.
  predictor_warning <- 0

  M <-200000
  alphaL <- (1-conf)/2
  alphaU <- 1- alphaL
  r2change <- (t*t)*(1-r2)/df2
  r2reduced <- r2 - r2change
  samplesize <- df1 + df2 + 1

  if (fixed == TRUE){
    #posterior of the partial correlation (postpr)
    cv<- sqrt(rchisq(n= M,df=(df2),ncp=0))
    z <- rnorm(n=M,mean=0,sd=1)
    pr <- (t*cv/sqrt(df2) + z)/sqrt(samplesize)
    postpr <- pr/sqrt(1+ pr*pr)

    #confidence distribution for R2reduced (CIDR2)
    if (df1 > 2) Cp1 <- rchisq(n= M,df=df1-2,ncp=0)
    if (df1 == 2) Cp1 <- 0
    Cnp1<- rchisq(n= M,df=df2+1,ncp=0)
    z <- rnorm(n=M ,mean=0,sd=1)
    samplesize <- df1 + df2 + 1

    B <- Cnp1*r2/(1-r2) - Cp1
    D <- (sqrt(B) + z)^2/samplesize
    D[is.na(D)] <- 0
    CIDR2 <- D/(1 + D)

    if (df1 > 1) CID <- postpr*(sqrt(1-CIDR2)/sqrt(1-r2x))
    if (df1 == 1) CID <- postpr

    initial <- quantile(CID,probs=c(alphaL, alphaU))
    initialqL <- initial[1]
    initialqU <- initial[2]
    lastqL <- initialqL
    lastqU <- initialqU

    lastdenL <- density(CID, from= lastqL, to = lastqL +1)$y[1]
    initialdenL <- lastdenL
    lastdenU <- density(CID, from= lastqU, to = lastqU +1)$y[1]
    initialdenU <- lastdenU

    for (n in 1:iter){
      #posterior of the partial correlation (postpr)
      cv<- sqrt(rchisq(n= M,df=(df2),ncp=0))
      z <- rnorm(n=M,mean=0,sd=1)
      pr <- (t*cv/sqrt(df2) + z)/sqrt(samplesize)
      postpr <- pr/sqrt(1+ pr*pr)

      #confidence distribution for R2reduced (CIDR2)
      if (df1 > 2) Cp1 <- rchisq(n= M, df=df1-2, ncp=0)
      if (df1 == 2) Cp1 <- 0
      Cnp1<- rchisq(n= M,df=df2+1,ncp=0)
      z <- rnorm(n=M ,mean=0,sd=1)
      samplesize <- df1 + df2 + 1

      B <- Cnp1*r2/(1-r2) - Cp1
      D <- (sqrt(B) + z)^2/samplesize
      D[is.na(D)] <- 0
      CIDR2 <- D/(1 + D)

      if (df1 > 1) CID <- postpr*(sqrt(1-CIDR2)/sqrt(1-r2x))
      if (df1 == 1) CID <- postpr

      fnL <- (1-1/n)*lastdenL + (1/n)*(sum(abs(CID-lastqL)<=(1/sqrt(n))))/(2*M*(1/sqrt(n)))
      SnL <- lastqL + (1/(n*max(fnL,initialdenL/sqrt(n))))*(alphaL - sum(CID <= lastqL)/M)

      fnU <- (1-1/n)*lastdenU + (1/n)*(sum(abs(CID-lastqU)<=(1/sqrt(n))))/(2*M*(1/sqrt(n)))
      SnU <- lastqU + (1/(n*max(fnU,initialdenU/sqrt(n))))*(alphaU - sum(CID <= lastqU)/M)
      lastqL <- SnL
      lastdenL <- fnL
      lastqU <- SnU
      lastdenU <- fnU
    }
  } else if (fixed == FALSE){
    #posterior of the partial correlation (postpr)
    cv <- sqrt(rchisq(n= M,df=(df2),ncp=0))
    cv1 <- sqrt(rchisq(n= M,df=(df2+1),ncp=0))
    z <- rnorm(n=M,mean=0,sd=1)
    pr <- (t*cv/sqrt(df2) + z)/(cv1)
    postpr <- pr/sqrt(1+ pr*pr)

    #confidence distribution for R2reduced (CIDR2)
    if (df1 > 2) Cp2 <- rchisq(n= M,df=df1-2,ncp=0)
    if (df1 == 2) Cp2 <- 0
    Cnp1<- rchisq(n= M,df=df2+1,ncp=0)
    Cn1<- rchisq(n= M ,df=df1+df2,ncp=0)
    z <- rnorm(n= M ,mean=0,sd=1)
    B <- Cnp1* r2reduced/(1-r2reduced) - Cp2
    D <- (sqrt(B) + z)^2/(Cn1)
    D[is.na(D)] <- 0
    CIDR2 <- D/(1 + D)

    if (df1 > 1) CID <- postpr*(sqrt(1-CIDR2)/sqrt(1-r2x))
    if (df1 == 1) CID <- postpr

    initial <- quantile(CID,probs=c(alphaL, alphaU))
    initialqL <- initial[1]
    initialqU <- initial[2]
    lastqL <- initialqL
    lastqU <- initialqU

    lastdenL <- density(CID, from= lastqL, to = lastqL +1)$y[1]
    initialdenL <- lastdenL
    lastdenU <- density(CID, from= lastqU, to = lastqU +1)$y[1]
    initialdenU <- lastdenU

    for (n in 1: iter){
      #posterior of the partial correlation (postpr)
      cv<- sqrt(rchisq(n= M,df=(df2),ncp=0))
      cv1<- sqrt(rchisq(n= M,df=(df2+1),ncp=0))
      z <- rnorm(n=M,mean=0,sd=1)
      pr <- (t*cv/sqrt(df2) + z)/(cv1)
      postpr <- pr/sqrt(1+ pr*pr)

      #confidence distribution for R2reduced (CIDR2)
      if (df1 > 2) Cp2 <- rchisq(n= M,df=df1-2,ncp=0)
      if (df1 == 2) Cp2 <- 0
      Cnp1<- rchisq(n= M,df=df2+1,ncp=0)
      Cn1<- rchisq(n= M ,df=df1+df2,ncp=0)
      z <- rnorm(n= M ,mean=0,sd=1)
      B <- Cnp1* r2reduced/(1-r2reduced) - Cp2
      D <- (sqrt(B) + z)^2/(Cn1)
      D[is.na(D)] <- 0
      CIDR2 <- D/(1 + D)

      if (df1 > 1) CID <- postpr*(sqrt(1-CIDR2)/sqrt(1-r2x))
      if (df1 == 1) CID <- postpr

      fnL <- (1-1/n)*lastdenL + (1/n)*(sum(abs(CID-lastqL)<=(1/sqrt(n))))/(2*M*(1/sqrt(n)))
      SnL <- lastqL + (1/(n*max(fnL,initialdenL/sqrt(n))))*(alphaL - sum(CID <= lastqL)/M)

      fnU <- (1-1/n)*lastdenU + (1/n)*(sum(abs(CID-lastqU)<=(1/sqrt(n))))/(2*M*(1/sqrt(n)))
      SnU <- lastqU + (1/(n*max(fnU,initialdenU/sqrt(n))))*(alphaU - sum(CID <= lastqU)/M)
      lastqL <- SnL
      lastdenL <- fnL
      lastqU <- SnU
      lastdenU <- fnU
    }
  } else {
    predictor_warning <- 1
    cat("Please specify fixed as TRUE or FALSE")
  }
  options(warn=0)   #resetting warning back to original state
  if (predictor_warning == 0) c(SnL,SnU)
}

#' Creates a confidence interval for the correlation.
#'
#' @param r The estimated correlation or partial correlation.
#' @param df The degrees of freedom.
#' @param conf The desired (1-\eqn{\alpha}) confidence level. Default is .95.
#' @param fixed Specifies whether the variable in the regression model is fixed or random. The default is for random predictor.
#' @param  iter Number of iterations for the stochastic approximation optimization. Default is 50.
#' @return Returns the lower and upper endpoints of the confidence interval.
#' @export
ci_r <- function(r, df, fixed=FALSE, conf=.95, iter=50){
	M <-200000
	alphaL <- (1-conf)/2
	alphaU <- 1- alphaL

	cv<- sqrt(rchisq(n= M,df=(df),ncp=0))
	cvp<- sqrt(rchisq(n= M,df=(df+1),ncp=0))
	z <- rnorm(n=M,mean=0,sd=1)

	if (fixed == FALSE){
		y <- (r*cv/sqrt(1 - r*r) + z)/(cvp)
	} else {
		y <- (r*cv/sqrt(1 - r*r) + z)/sqrt(df + 2)
	}

    CID <- y/sqrt(1 + y*y)

	initial <- quantile(CID,probs=c(alphaL, alphaU))
	initialqL <- initial[1]
	initialqU <- initial[2]
	lastqL <- initialqL
	lastqU <- initialqU

	lastdenL <- density(CID, from= lastqL, to = lastqL +1)$y[1]
	initialdenL <- lastdenL
	lastdenU <- density(CID, from= lastqU, to = lastqU +1)$y[1]
	initialdenU <- lastdenU

	for (n in 1: iter){
		cv<- sqrt(rchisq(n= M,df=(df),ncp=0))
		cvp<- sqrt(rchisq(n= M,df=(df+1),ncp=0))
		z <- rnorm(n=M,mean=0,sd=1)

		if (fixed == FALSE){
			y <- (r*cv/sqrt(1 - r*r) + z)/(cvp)
		} else {
			y <- (r*cv/sqrt(1 - r*r) + z)/sqrt(df + 2)
		}

   		CID <- y/sqrt(1+y*y)

		fnL <- (1-1/n)*lastdenL + (1/n)*(sum(abs(CID-lastqL)<=(1/sqrt(n))))/(2*M*(1/sqrt(n)))
		SnL <- lastqL + (1/(n*max(fnL,initialdenL/sqrt(n))))*(alphaL - sum(CID <= lastqL)/M)

		fnU <- (1-1/n)*lastdenU + (1/n)*(sum(abs(CID-lastqU)<=(1/sqrt(n))))/(2*M*(1/sqrt(n)))
		SnU <- lastqU + (1/(n*max(fnU,initialdenU/sqrt(n))))*(alphaU - sum(CID <= lastqU)/M)
		lastqL <- SnL
		lastdenL <- fnL
		lastqU <- SnU
		lastdenU <- fnU
	}
	c(SnL,SnU)
}
