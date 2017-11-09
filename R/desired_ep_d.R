#' Determines the required sample size in a future study to achieve the desired expected power (ep)
#' based on the uncertainty associated with an existing study.
#' Uses the observed standardized mean difference from the previous study as input in the sample size planning process.
#' Assumes both that the previous and future study will involve fixed predictors. If this is not the case, use
#' the desired_ep_t() function for full flexibility.
#'
#' @param d The observed standardized mean difference of the previous study based on
#' a pooled standard deviation.
#' @param n1 The number of observations in the first group.
#' @param n2 The number of observations in the second group.
#' @param desired_ep The desired expected power for the future study.
#' @param alpha The significance level. Default is \eqn{\alpha = .05}.
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
#' @return Returns (1) the sample size per group required for the future study to achieve the specified level of expected power.
#' This reflects the uncertainty associated with the previous study and (2) the median 95% confidence interval width for the correlation in the prospective study.
#' @export
#' @examples
#' \dontrun{
#' desired_ep_d(d = 0.50, n1=25, n2=25, desired_ep = 0.80)
#' }
desired_ep_d <- function(d, n1, n2, desired_ep=0.80, alpha=0.05, filter=0, upper_null=0){
	exit <- 0
	if (n1 < 3) {
		cat("Degrees of freedom for each group (n1 and n2) must be greater than 3.\n")
		exit<-1
		}
	if (n2 < 3) {
		cat("Degrees of freedom for each group (n1 and n2) must be greater than 3.\n")
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
	if (desired_ep > 1) {cat("desired_ep must be between alpha and 1. Default value is 0.80 if not specified.\n")
		exit<-1
		}
	if (desired_ep < alpha) {cat("desired_ep must be between alpha and 1. Default value is 0.80 if not specified.\n")
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

	if(exit == 0){
	  start_n <- c(3,4,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,110,120,130,140,
	               150,160,175,185,200,225,250,275,300,350,400,450,500,
	               550,600,650,700,750,800,850,900,950,1000,
	               1200,1400,1600,1800,2000,2400,2600,2800,3000,
	               4000, 5000, 6000, 7000, 8000, 9000, 10000,
	               12000, 14000, 16000, 18000, 20000, 25000, 30000, 40000, 50000)

	#Initial model computing expected power across the starting values ranging from n=3 to 50,000
	posterior <- posterior_d(d=d, n1=n1, n2=n2, filter=filter, upper_null=upper_null)  #Single static posterior distribution for this initial model
	ndraws<- length(posterior)

	start_values <- matrix(NA, nrow=length(start_n), ncol=2)
	for (i in 1:length(start_n)){
			n <- start_n[i]				#Sample size per group
			dfnew <- 2*n - 2

			posterior_lambda <- posterior/sqrt(2/n)

			z <- rnorm(n=ndraws,mean=0,sd=1)
			cdf <-  sqrt(rchisq(n= ndraws,df= dfnew,ncp=0))
			cdf1 <-  sqrt(rchisq(n= ndraws,df= dfnew + 1,ncp=0))

			#Posterior Predictive Distribution
			ppt <- (z + posterior_lambda)/(cdf/sqrt(dfnew))
			ppt <- abs(ppt)

			dist <- ecdf(ppt)
			EP <- 1 - dist(qt(alpha/2,df=dfnew,ncp=0,lower.tail=FALSE))
			start_values[i, 1] <- start_n[i]
			start_values[i, 2] <- EP
		}

	#Estimate a smoothed curve across the range of sample sizes in start_n.
	#We then use the predicted values from this model and solve for DesirePower and a lower and higher value
	#of expected power around that. This provides the initial range of sample sizes to conduct a more
	#detailed and precise estimate of expected power and a new smoothed spline model.
	Starting_Expected_Power <- smooth.spline(x= start_values[,1],y= start_values[,2])

	Upper_Power_Value <- min(c((1 - desired_ep)*.25 + desired_ep,.99))
	Lower_Power_Value <- .93*desired_ep
	if (predict(Starting_Expected_Power, max(start_n))$y < Upper_Power_Value){
		cat("Required sample size is too high to estimate accurately.\n")
		cat("Sample size of ", max(start_n)," results in initial estimated statistical power of ", predict(Starting_Expected_Power, max(start_n))$y, ".\n")
	}

	f_start <- function(x) (predict(Starting_Expected_Power, x)$y - desired_ep)
	f_high  <- function(x) (predict(Starting_Expected_Power, x)$y - Upper_Power_Value)
	f_low   <- function(x) (predict(Starting_Expected_Power, x)$y - Lower_Power_Value)

	#Here we define a fairly tight range around the solution to do more refined modeling of df and expected power
	Starting_Solution <- uniroot(f_start,interval=c(start_values[1,1],  start_values[NROW(start_values),1]),tol = 0.001)$root
	Upper <- uniroot(f_high,interval=c(start_values[1,1],  start_values[NROW(start_values),1]),tol = 0.001)$root
	Lower <- uniroot(f_low,interval=c(start_values[1,1],  start_values[NROW(start_values),1]),tol = 0.001)$root

	#Evenly sample 80 df points from the lower to upper values
	final_n1 <- seq(from=Lower-2, to=Upper+1, by=(Upper - Lower)/80)

	#Add more concentration around the initial solution.
	#This initial solution is often quite decent so sampling intensively around this initial estimate
	#provides data to really help the precision of the estimate of expected statistical power.

	final_n2 <- rnorm(n = 40, mean = Starting_Solution, sd = (Upper-Starting_Solution)/20 )
	final_n <- c(final_n1, final_n2)
	final_n[final_n < 3] <- 3


	val <- matrix(NA, nrow=length(final_n), ncol=4)
	for (i in 1:length(final_n) ){
			n_new <- final_n[i]
			dfnew <- 2*n_new - 2
			posterior <- posterior_d(d=d, n1= n1, n2= n2, filter=filter, upper_null=upper_null) #obtaining a new posterior for every evaluated df
			posterior_lambda <- posterior/sqrt(2/n_new)

			ndraws <- length(posterior)
			z <- rnorm(n=ndraws,mean=0,sd=1)
			cdf <-  sqrt(rchisq(n= ndraws,df= dfnew, ncp=0))
			cdf1 <-  sqrt(rchisq(n= ndraws,df= dfnew + 1,ncp=0))

			#Posterior Predictive Distribution
			ppt <- (z + posterior_lambda)/(cdf/sqrt(dfnew))
			ppt <- abs(ppt)
			CI_Interval <- ci_d(d=quantile(abs(posterior),probs=c(.50)), n1=n_new, n2=n_new, conf=.95, iter=5)
			val[i, 2]<- CI_Interval[2] - CI_Interval[1]

			dist <- ecdf(ppt)
			EP <- 1 - dist(qt(alpha/2,df=dfnew,ncp=0,lower.tail=FALSE))
			val[i ,1]<- d
			val[i, 3]<- final_n[i]
			val[i, 4]<- EP

	}
	power <- smooth.spline(x=val[,3],y=val[,4])
	CIW <- smooth.spline(x=val[,3],y=val[,2])

	f_final <- function(x) (predict(power, x)$y - desired_ep)
	Desired_Power_n <- uniroot(f_final,interval=c(min(val[,3]),  max(val[,3])),tol = 0.001)$root
	lowern <- round(min(final_n))
	uppern <- round(max(final_n))

	upper_power <- 	predict(power, uppern)$y
	lower_power <- 	predict(power, lowern)$y
	upper_CIW <- 	predict(CIW, uppern)$y
	lower_CIW <- 	predict(CIW, lowern)$y

	plot(power, xlab="Sample Size per Group for Future Study",ylab="Expected Statistical Power",xlim=c(lowern, uppern),ylim=c(lower_power, upper_power), frame.plot=FALSE,type="l",lwd=2)
	segments(lowern, desired_ep, Desired_Power_n, desired_ep,col= "grey40", lwd=1,lty="dashed")
	segments(Desired_Power_n, min(val[,4]), Desired_Power_n, desired_ep,col= "grey40", lwd=1,lty="dashed")
	text(uppern,predict(power, uppern)$y,"Expected Power",adj=c(1,0),cex=.7)
	text(Desired_Power_n ,predict(power, lowern)$y,bquote("N ="~.(Desired_Power_n)), adj=c(0,1),cex=.7)

	par(new=TRUE)
	plot(CIW,type="l",lty="dashed",lwd=2,main=NA,col=c("grey60"),xlim=c(lowern, uppern),ylim=c(upper_CIW, lower_CIW),frame.plot=FALSE,axes=FALSE,xlab=NA, ylab=NA)
	axis(side=4,at=NULL,tick=TRUE,outer=FALSE)
	text(uppern,(max(val[,2])-min(val[,2]))/2+min(val[,2]),"Median Confidence Interval Width",srt=90)
	segments(Desired_Power_n,  predict(CIW, uppern)$y, Desired_Power_n, predict(CIW, Desired_Power_n)$y,col= "grey40", lwd=1,lty="dashed")
	segments(Desired_Power_n, predict(CIW, Desired_Power_n)$y, uppern, predict(CIW, Desired_Power_n)$y,col= "grey40", lwd=1,lty="dashed")
	text(lowern,predict(CIW, lowern)$y,"95% CI Width",adj=c(0,1),cex=.7)

	ExpectedPower <- as.data.frame(matrix(NA, nrow = 1, ncol = 2))
	colnames(ExpectedPower) <- c("GroupSampleSize", "Median.95CI.Width")
	ExpectedPower$GroupSampleSize <- Desired_Power_n
	ExpectedPower$Median.95CI.Width <- predict(CIW, Desired_Power_n)$y

	return(ExpectedPower)
	}
}


