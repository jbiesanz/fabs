#' Estimates a series of robust regressions for a specified regression model.
#'
#' @param data The dataset for the analysis.
#' @param model A regression model from lm()
#' @param R The number of resamples to be conducted.
#' @param conf The confidence level for intervals.
#' @return Returns (1) the lm() regression model with HC4 heteroscedastic covariance consistent standard errors
#' (2) p-values based on the wild bootstrap (and HC4 standard errors) under the null hypothesis.
#' (3) percentile and BCa confidence intervals using case-wise resampling
#' (4) Robust regression analysis from rlm() under the default Huber loss function.
#' (5) Robust regression with resampled wild bootstrap standard errors.
#' @export
#' @import MASS
#' @import boot
#' @examples
#' \dontrun{
#' x <- c(1,2,3,4,5,6)
#' y <- c(2.1,1.9,2.1,1.8,10,2.2)
#' temp_data <- as.data.frame(cbind(y, x))
#' lm_model <- lm(y~x, data=temp_data)
#' robust_lm_inferences(data=temp_data, model=lm_model)
#' }
robust_lm_inferences <- function(data, model, R=9999, conf=.95){
  formula <- formula(model)
  reg <- lm(formula, data=data, model=TRUE, x=TRUE, y=TRUE)
  x <- reg$x
  y <- reg$y
  p1 <- ncol(x)

  HC4 <- lm_HC4(x, y, alpha = (1-conf))
  Wild_Boot_pvalues <- HC4_wildboot_pvalues(x, y, R)
  resampling <- boot(data=x, statistic= lm_boot_casewise,
                     R=R,y=y, parallel ="multicore", ncpus=(detectCores()-1))
  boot_out <- matrix(NA, nrow=p1, ncol=5)
  for (i in 1:p1){
    perc <- boot.ci(resampling, index=i,conf=conf, type="perc")
    bca <-  boot.ci(resampling, index=i,conf=conf, type="bca")
    boot_out[i,1] <- coef(reg)[i]
    boot_out[i,2] <- perc$percent[4]
    boot_out[i,3] <- perc$percent[5]
    boot_out[i,4] <- bca$bca[4]
    boot_out[i,5] <- bca$bca[5]
  }
  colnames(boot_out)<-c("Estimate","perc.lower","perc.upper","bca.lower","bca.upper")
  rownames(boot_out) <- names(coef(reg))

  r_lm <- rlm(formula, data)
  r_output <- as.data.frame(summary(r_lm)$coefficients)
  colnames(r_output) <- c("Estimate","Std_Error","t_statistic")
  r_output$df <- length(y) - p1
  r_output$p_value <- 2*pt(abs(r_output[,3]), df= length(y) - p1, lower.tail=FALSE)

  r_wild_output <- r_output
  r_wild_output$Std_Error <- robust_wild_se(x, y, R=R)
  r_wild_output$t_statistic <- r_wild_output$Estimate/r_wild_output$Std_Error
  r_wild_output$p_value <- 2*pt(abs(r_wild_output$t_statistic), df= r_wild_output$df, lower.tail=FALSE)

  list(OLS_HC4_Results = HC4[[1]],
       wild_HC4_pvalues = Wild_Boot_pvalues[[1]],
       ResamplingPairwise=boot_out,
       Robust_Regression = r_output,
       Wild_Robust_Regression = r_wild_output)
}

lm_boot_casewise <- function(data, y, indices){
  x <- as.matrix(data[indices,]) #data here is the predictor matrix
  y <- y[indices]
  ols <- lsfit(x[,-1], y, intercept=TRUE) #Using lsfit to speed up the resampling
  coef(ols)
}

lm_HC4 <- function(x, y, alpha=.05){
  # Compute confidence intervals via least squares
  # regression using heteroscedastic method HC4
  # recommended by Cribari-Neto (2004).
  x <- as.matrix(x)
  p1 <- ncol(x)	#number of coefficients in the regression model
  if (p1 == 1) ols <- lm(y ~ 0 + x)
  if (p1 > 1) ols <- lsfit(x[,-1], y, intercept=TRUE) #faster than lm()
  xtx <- solve(t(x)%*%x)
  h <- diag(x%*%xtx%*%t(x))
  n <- length(h)
  d <- (n*h)/sum(h)
  d[d > 4] <- 4	#d is the minimum of {4, d}

  hc4 <- xtx%*%t(x)%*%diag(ols$res^2/(1-h)^d)%*%x%*%xtx
  df <- nrow(x) - p1
  crit <- qt(1-alpha/2,df)
  ci <- matrix(NA,nrow=ncol(x),ncol=7)

  for(j in 1:ncol(x)){
    ci[j,1]<-ols$coef[j]
    ci[,2]=sqrt(diag(hc4))
    ci[j,3]<-ols$coef[j]-crit*sqrt(hc4[j,j])
    ci[j,4]<-ols$coef[j]+crit*sqrt(hc4[j,j])
    ci[j,5]<-ci[j,1]/ci[j,2]
    ci[j,6]<-df
    ci[j,7]<-2*(1-pt(abs(ols$coef[j]/sqrt(hc4[j,j])),df))
  }
  dimnames(ci)<-list(colnames(x),c("Estimate","Std. Error","lower.ci","upper.ci","t-statistic","df","p-value"))
  list(OLS_HC4 = ci, cov = hc4)
}


HC4_wildboot_pvalues <- function(x, y, R=9999){
  x <- as.matrix(x)
  p1 <-ncol(x)	#number of coefficients in the regression model
  if (p1 == 1) temp <- lm(y ~ 0 + x)
  if (p1 > 1) temp <- lsfit(x[,-1], y, intercept=TRUE) #faster than lm()
  yhat <- 0         #set equal to 0 to provide inferences for the intercept.
  res <- temp$residuals
  s <- lm_HC4(x, y)$cov
  b <- temp$coef
  sample_ttest <- abs(b)/sqrt(diag(as.matrix(s)))
  vstar <- matrix(ifelse(rbinom(length(y)*R,1,0.5)==1,-1,1),nrow=R)
  wild_tstar <- apply(vstar,1,olswbtest.sub,yhat,res,x) #rvalb is a p by R matrix
  if(p1==1)  abs_wild_tstar <- t(as.matrix(abs(wild_tstar)))
  if(p1 > 1) abs_wild_tstar <- as.matrix(abs(wild_tstar))
  pvals=NA
  for(j in 1:p1) pvals[j] <- mean(abs_wild_tstar[j,] >= sample_ttest[j])
  wildHC4_pvalue <- cbind(b,pvals)
  colnames(wildHC4_pvalue) <- c("Estimate","p-value")
  list(wild_HC4_pvalue = wildHC4_pvalue)
}

olswbtest.sub<-function(vstar,yhat,res,x){
  ystar <- yhat + res*vstar
  p1 <- ncol(x)
  if (p1 == 1) vals <- t(as.matrix(lm(ystar ~ 1)$coef[1:p1]))
  if (p1 > 1) vals <- t(as.matrix(lsfit(x[,-1],ystar)$coef[1:p1]))
  sa <- lm_HC4(x, ystar)$cov
  test <- vals/sqrt(diag(as.matrix(sa)))
  test
}

robust_wild_se <- function(x, y, R=9999){
  r_lm <- rlm(y ~ 0 + x)
  resid <- r_lm$residuals
  yhat <- r_lm$fitted

  vstar <- matrix(ifelse(rbinom(length(y)*R,1,0.5)==1,-1, 1),nrow=R)  #R by n matrix
  r_wild_tstar <- apply(vstar, 1, rlm.sub, yhat, resid, x) #rvalb is a p1 by R matrix

  #Robust regression wild standard errors
  if (ncol(x) ==1 ) r_wild_se <- var(r_wild_tstar)
  if (ncol(x) > 1 ) r_wild_se <- sqrt(diag(cov(t(r_wild_tstar))))
  r_wild_se
}

rlm.sub <- function(vstar, yhat, resid, x){
  ystar <- yhat + resid*vstar
  r_lm_star <- rlm(ystar ~ 0 + x)
  r_lm_star$coef
}

#' Estimates the robust standardized mean difference based on trimmed means and pooled winsorized standard deviation.
#'
#' @param x Vector of observations from the first group.
#' @param y Vector of observations from the second group.
#' @param trim The proportion of observations to trim from each tail within each group.
#' @return Returns the robust standardized mean difference $d_r$ from Algina, Keselman, and Penfield (2005)
#' based on the trimmed means and the pooled winsorized standard deviation. The estimate for the robust standardized
#' mean difference is adjusted to be equal to that of the standardized mean difference based on normal distribution
#' with equal variances. The adjustment factor depends on the level of trimming.
#' @export
#' @import psych
#' @examples
#' x <- c(-2.40, -1.87, -0.60, -0.54, -0.12, -0.02, 0.12, 0.34, 0.40, 0.53, 0.55, 0.62, 0.92, 1.21,
#'         1.49, 1.55, 1.57, 1.57, 1.82, 1.87, 1.90, 1.91, 1.93, 2.34, 2.37)
#' y <- c(-1.32, -1.25, -0.91, -0.62, -0.55, -0.41,-0.40, -0.31, -0.28, -0.21, -0.18, -0.16, -0.03,
#'        -0.02, 0.04, 0.22, 0.38, 0.51, 0.53, 0.61, 1.09, 1.47, 1.59, 2.39, 2.47)
#' d_robust(x, y, trim=.20)
d_robust <- function(x, y, trim=.20){
  n_x <- length(x)
  n_y <- length(y)

  x_trim_mean <- winsor.mean(x, trim=trim, na.rm=TRUE)
  y_trim_mean <- winsor.mean(y, trim=trim, na.rm=TRUE)
  x_win_var <- winsor.var(x, trim=trim, na.rm=TRUE)
  y_win_var <- winsor.var(y, trim=trim, na.rm=TRUE)

  win_var_pooled <- ((n_x-1)*x_win_var+(n_y-1)*y_win_var)/(n_x+n_y-2)
  #The formula for the adjustment as a function of the level of trimming was determined through extensive simulation
  #for trim levels ranging from .01 to .48.
  adj <- .9993 -1.7402*trim -.0512*trim^2 -.9266*trim^3

  d_r <- adj*(x_trim_mean - y_trim_mean)/sqrt(win_var_pooled)
  d_r
}
