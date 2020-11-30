#' Wrapper function to estimate an lm() model in lavaan under full information maximum likelihood to account for missing data.
#'
#' @param data The dataset for the analysis.
#' @param model A regression model from lm()
#' @return Returns (1) a dataset with the ML regression estimates under FIML assuming either missing at random
#' or missing completely at random, standard errors, t-test statistic, p-values under t-distribution,
#' gamma (estimated fraction of missing data), N.effective (estimated equivalent complete data sample size),
#' and df = n*(1-gamma) where n is the number of rows in the dataset. Both N.effective and df are rounded down.
#' (2) sigma which estimates the residual standard error.
#' @export
#' @examples
#' \dontrun{
#' x <- c(1,2,3,4,5,NA,NA,7,7,7,7)
#' y <- c(2.1,NA,2.1,1.8,2,2.2,4,NA,7,7,7)
#' temp_data <- as.data.frame(cbind(y, x))
#' lm_model <- lm(y~x, data=temp_data)
#' fiml.regression(data=temp_data, model=lm_model)
#' }
fiml.regression <- function(data, model) {
  sem_model <- lavaan::sem(format(formula(model)), data, missing = "fiml", estimator = "ML", fixed.x=F)

  #Extract the name of the dependent variable.
  dv <- model$terms[[2]]

  #Trimming output to just the regression coefficients and residual SD
  sem.output <- subset(lavaan::parameterEstimates(sem_model),
                       lavaan::parameterEstimates(sem_model)$lhs == dv)
  #sem.output <- subset(sem.output, sem.output$rhs != dv)
  colnames(sem.output)[colnames(sem.output) == "rhs"] <- "Term"
  colnames(sem.output)[colnames(sem.output) == "est"] <- "Coefficient"
  colnames(sem.output)[colnames(sem.output) == "z"] <- "t.test"
  sem.output$Coefficient[sem.output$Term == dv] <- sqrt(sem.output$Coefficient[sem.output$Term == dv])
  sem.output$Term[sem.output$Term == dv] <- "(Residual.SD)"

  sem.output$lhs <- NULL
  sem.output$Term[sem.output$op == "~1"] <- "(Intercept)"
  sem.output$op <- NULL

  #removing the CI endpoints since they are based on normal theory not
  #the t-distribution
  sem.output$ci.lower <- NULL
  sem.output$ci.upper <- NULL

  #Obtaining gamma (fraction of missing information)
  #Following Savalei & Rhemtulla (2012, SEM)
  #Step 1. Obtain the model-implied covariance matrix and means
  #        Use the finite sample adjustment for the covariance matrix.
  n <- lavaan::nobs(sem_model)
  Sigma.hat <- lavaan::fitted.values(sem_model)$cov * n/(n - 1)
  mu.hat <- lavaan::fitted.values(sem_model)$mean

  model.implied <- lavaan::sem(format(formula(model)), sample.cov = Sigma.hat,
                       sample.mean = mu.hat, sample.nobs = n,
                       std.lv = TRUE, meanstructure = TRUE,
                       information = "observed")
  model.implied.output <- subset(lavaan::parameterEstimates(model.implied),
                                 lavaan::parameterEstimates(model.implied)$lhs == dv)
  #model.implied.output <- subset(model.implied.output, model.implied.output$rhs != dv)

  sem.output$gamma <- 1 - (model.implied.output$se^2/sem.output$se^2)
  sem.output$N.effective <- trunc(n * (1 - sem.output$gamma))

  #Obtaining df
  complete.data.df <- n - nrow(as.matrix(coef(lm(formula(model), data = data))))
  sem.output$df <- trunc(complete.data.df * (1 - sem.output$gamma))

  #recomputing p-values based on the t-test with df given above.
  sem.output$pvalue <- 2 * pt(abs(sem.output$t.test),
                              df = sem.output$df, lower.tail = FALSE)

  #return the output sorted by term so (Intercept) matches lm() output and norm.regression()
  reordered.output <- sem.output[order(sem.output$Term), ]
  rownames(reordered.output) <- c(1:nrow(reordered.output))
  reordered.output[, c(1:4, 8, 5:7)]

  sigma <- reordered.output[which(grepl("(Residual.SD)", reordered.output$Term)),2]
  #Remove row with (Residual.SD)
  fiml.summary <- reordered.output[-which(grepl("(Residual.SD)", reordered.output$Term)),]

  list(fiml.summary=fiml.summary, sigma=sigma)
}

#' Function to standardize all variables based on supplied mean and sd vectors
#'
#' @param data The dataset
#' @param mean Vector of means.
#' @param sd Vector of standard deviations
#' @return Returns a dataset where all of the variables have been scaled.
#' @export
#' @examples
#' \dontrun{
#' x <- c(1,2,3,4,5)
#' y <- c(2.1,2.5,4,5,6)
#' temp_data <- as.data.frame(cbind(y, x))
#' scale.all.variables(data=temp_data, mean=sapply(temp_data, mean), sd=sapply(temp_data, sd))
#' }
scale.all.variables <- function(data, mean, sd) {
  for (i in 1:length(data)) {
    data[, i] <- (data[, i] - means[i])/sds[i]
  }
  data
}

#' Function to return the EM means and standard deviations using the norm package.
#'
#' @param data The dataset for the analysis.
#' @return Returns (1) the EM means and (2) the EM standard deviations for the dataset.
#' @export
#' @import norm
#' @examples
#' \dontrun{
#' x <- c(1,2,3,4,5,NA,NA,7,7,7,7)
#' y <- c(2.1,NA,2.1,1.8,2,2.2,4,NA,7,7,7)
#' temp_data <- as.data.frame(cbind(y, x))
#' em.summary(temp_data)
#' }
em.summary <- function(data) {
  rngseed(sample(1:1e+08, 1))
  n <- nrow(data)
  s <- prelim.norm(as.matrix(data))
  thetahat <- em.norm(s, showits=FALSE)
  mean <- as.data.frame(t(getparam.norm(s, thetahat)$mu))
  sd <- as.data.frame(t(sqrt((diag(getparam.norm(s, thetahat)$sigma)))))
  colnames(mean) <- names(data)
  colnames(sd) <- names(data)
  output <- list(mean, sd)
  names(output) <- c("mean", "sd")
  return(output)
}

#' Bootstraps the full information maximum likelihood regression model under lavaan.
#'
#' @param data The dataset for the analysis.
#' @param model A regression model from lm()
#' @param R The number of resamples to be conducted.
#' @param conf The confidence level for intervals.
#' @param z.function A user supplied function to standardize variables and create nonlinear terms. Called internally.
#' @return Returns (1) maximum likelihood regression coefficient estimates
#' (s) percentile and BCa confidence intervals using case-wise resampling
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
boot.fiml <- function(data, model, z.function = FALSE, R=9999, conf=.95) {
  formula <- format(formula(model))
  #Extract estimates based on the original
  fiml.reg <- fiml_boot_casewise(data=data,formula=formula, indices=seq(1:nrow(data)),z.function=z.function)
  resampling <- boot(data=data, statistic= fiml_boot_casewise,
                     R=R, formula=formula, z.function =z.function,
                     parallel ="multicore", ncpus=(detectCores()-1))
  boot_out <- matrix(NA, nrow=length(fiml.reg), ncol=5)
  for (i in 1:length(fiml.reg)){
    perc <- boot.ci(resampling, index=i,conf=conf, type="perc")
    bca <-  boot.ci(resampling, index=i,conf=conf, type="bca")
    boot_out[i,1] <- fiml.reg[i]
    boot_out[i,2] <- perc$percent[4]
    boot_out[i,3] <- perc$percent[5]
    boot_out[i,4] <- bca$bca[4]
    boot_out[i,5] <- bca$bca[5]
  }
  colnames(boot_out)<-c("Estimate","perc.lower","perc.upper","bca.lower","bca.upper")
  rownames(boot_out) <- fiml.regression(data, model)$fiml.summary$Term
  boot_out
}

fiml_boot_casewise <- function(data, indices, formula, z.function = FALSE) {
  bootdata <- data[indices, ]

  #Using lm to extract the name of the dv.
  dv <- lm(formula, data = bootdata)$terms[[2]]

  if (is.function(z.function) == FALSE){
    b.model <- lavaan::sem(formula, bootdata, missing = "fiml", estimator = "ML")

    #Trimming output to just the regression coefficients
    sem.output <- subset(lavaan::parameterEstimates(b.model), lavaan::parameterEstimates(b.model)$lhs == dv)
    sem.output <- subset(sem.output, sem.output$rhs != dv)
  }
  else {
    b.model <- lavaan::sem(formula, bootdata, missing = "fiml", estimator = "ML")
    n <- lavaan::nobs(b.model)
    ml.sigma <- sqrt(diag(lavaan::fitted.values(b.model)$cov * n/(n - 1)))
    ml.sigma <- as.data.frame(t(ml.sigma), col.names= names(ml.sigma))
    ml.mu <- lavaan::fitted.values(b.model)$mean
    ml.mu <- as.data.frame(t(as.vector(ml.mu)))
    names(ml.mu) <- names(ml.sigma)
    zbootdata <- z.function(bootdata, mean = ml.mu, sd = ml.sigma)
    z.model <- lavaan::sem(formula, zbootdata,  missing = "fiml", estimator = "ML")
    sem.output <- subset(lavaan::parameterEstimates(z.model), lavaan::parameterEstimates(z.model)$lhs == dv)
    sem.output <- subset(sem.output, sem.output$rhs != dv)
  }
  colnames(sem.output)[colnames(sem.output) == "rhs"] <- "Term"
  colnames(sem.output)[colnames(sem.output) == "est"] <- "Coefficient"
  sem.output$Term[sem.output$op == "~1"] <- "(Intercept)"
  reordered.output <- sem.output[order(sem.output$Term), ]
  #Remove row with (Residual.SD)
  boot.fiml.summary <- reordered.output[-which(grepl("(Residual.SD)", reordered.output$Term)),]
  return(reordered.output$Coefficient)
}

#' Function to automate multiple imputation (MI) for missing data using the norm package for an lm() regression model.
#'
#' @param data The dataset for the analysis.
#' @param model A regression model from lm()
#' @param m Number of imputations to conduct
#' @param standardize Whether or not to standardize each imputed dataset before rerunning the lm() model on that imputed dataset.
#' @param digits Number of digits to print.
#' @return Returns a dataset with the lm() regression terms, average regression coefficients under MI, standard errors, t.test,
#' df calculated using Barnard and Rubin (1999), gamma   ML regression estimates under FIML assuming either missing at random
#' or missing completely at random, standard errors, t-test statistic, p-values under t-distribution,
#' gamma (estimated fraction of missing data), N.effective (estimated equivalent complete data sample size),
#' and efficiency = 1/(1 + gamma/m).
#' and df = n*(1-gamma) where n is the number of rows in the dataset. Both N.effective and df are rounded down.
#' @export
#' @import norm
#' @examples
#' \dontrun{
#' x <- c(1,2,3,4,5,NA,NA,7,7,7,7)
#' y <- c(2.1,NA,2.1,1.8,2,2.2,4,NA,7,7,7)
#' temp_data <- as.data.frame(cbind(y, x))
#' lm_model <- lm(y~x, data=temp_data)
#' norm.regression(data=temp_data, model=lm_model)
#' }
norm.regression <- function(data, model, m = 100, standardize = F, digits=6){
  formula <- formula(model)
  if (ncol(data) > 30) cat("Consider reducing your dataset down to a smaller number of variables (e.g., less than 30) \nas norm can experience difficulty with large datasets.")
  if (sum(sapply(data, is.numeric)) < ncol(data)) cat("Only numeric variables are appropriate for norm. Please reduce the dataset to numeric or integeter variables.\n")
  if (rankMatrix(cov(data, use="pairwise.complete.obs"))[1] == ncol(data)){
    rngseed(sample(1:1e+08, 1))
    model_vars <- as.data.frame(model.matrix(as.formula(formula), model.frame(formula, data, na.action = NULL)))
    if (var(model_vars[, 1]) == 0)
      model_vars$"(Intercept)" <- NULL
    model_vars$rows <- rownames(model_vars)
    data$rows <- rownames(data)
    commonvars <- intersect(colnames(model_vars), colnames(data))
    newdata <- merge(model_vars, data, by = commonvars)
    newdata$rows <- NULL

    p1 <- nrow(as.matrix(coef(lm(formula, data = newdata))))
    s <- prelim.norm(as.matrix(newdata)) # this line has no output
    thetahat <- em.norm(s, showits = FALSE)
    theta <- da.norm(s, thetahat, steps = 50000)

    b <- matrix(nrow = m, ncol = p1)
    vb <- matrix(nrow = m, ncol = p1)

    for (j in 1:m) {
      theta1 <- da.norm(s, theta, steps = 100)
      ximp <- imp.norm(s, theta1, newdata)
      ximp <- as.data.frame(ximp)

      if (standardize == TRUE)
        model <- lm(formula, data = as.data.frame(scale(ximp)))
      else model <- lm(formula, data = ximp)

      for (k in 1:p1) {
        b[j, k] <- coefficients(model)[k]
        vb[j, k] <- vcov(model)[k, k]
      }
      theta <- theta1
    }

    # Combining results and computing df from Barnard & Rubin (1999)
    MI_Output <- as.data.frame(matrix(NA, nrow = p1, ncol = 9))
    colnames(MI_Output) <- c("Term", "Coefficient", "se", "t.test", "df", "pvalue", "gamma", "N.effective", "efficiency")
    for (c in 1:p1) {
      n <- NROW(ximp)
      vcom <- n - p1
      Qbar <- mean(b[, c])
      B <- var(b[, c])
      Ubar <- mean(vb[, c])
      Tvar <- Ubar + (1 + 1/m) * B
      rel.incr <- (1 + 1/m) * B/Ubar
      vm <- (m - 1) * (1 + 1/rel.incr)^2
      gamma <- (rel.incr + 2/(vm + 3))/(rel.incr + 1)
      vobs <- ((vcom + 1)/(vcom + 3)) * vcom * (1 - gamma)
      dfmi <- 1/(1/vm + 1/vobs)
      MI_Output$Term[c] <- rownames(as.data.frame(coef(model)))[c]
      MI_Output$Coefficient[c] <- round(Qbar, digits)
      MI_Output$se[c] <- round(sqrt(Tvar), digits)
      MI_Output$t.test[c] <- round(Qbar/sqrt(Tvar), digits)
      MI_Output$df[c] <- trunc(dfmi)
      MI_Output$pvalue[c] <- round(2 * pt(abs(Qbar/sqrt(Tvar)), df = dfmi, lower.tail = FALSE), digits)
      MI_Output$gamma[c] <- round(gamma, digits)
      MI_Output$N.effective[c] <- trunc(n * (1 - gamma))
      MI_Output$efficiency[c] <- round(1/(1 + gamma/m), digits)
    }
    MI_Output[order(MI_Output$Term), ]
  } else {
    cat("Error: Norm requires a dataset with only continuous variables and a dataset of full rank (e.g., no sets of variables that are perfectly correlated).")
  }
}
