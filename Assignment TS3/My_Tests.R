
library(ggplot2)
library(DescTools)
library(lmtest)
library(ggplot2)
library(fBasics)
library(fpp3)
library(forecast)
library(urca)
library(fUnitRoots)
library(TSA)

# - the functions below are just simple wrappers for various time-series data 
#   tests to help me understand the results in a more layman's perspective =)
#
# - (most) tests below are tested against their respective null hypothesis (H0):
#     TRUE: fail to reject H0
#     FALSE: reject H0



myttest <- function(c, print_output=TRUE) {
  # T-test: test for mean 0, linear trend
  ttest_rm <- t.test(c)
  # T-Test for mean 0: check if 95% confidence interval ($conf.int[1] lower, $conf.int[2] upper) 
  # contains 0 and p-value ($p.value) >= 0.05
  # TRUE (H0): mean is zero, removed linear trend -> fail to reject H0
  # FALSE (HA): mean not zero, linear trend present -> reject H0
  result <- ttest_rm$conf.int[1] <=0 && ttest_rm$conf.int[2] >=0 && ttest_rm$p.value >= 0.05
  if (print_output) {
    print(ttest_rm)
    if(result) print("T-Test: mean is zero, removed linear trend -> fail to reject H0")
    else print("T-Test: mean not zero, linear trend present -> reject H0")
  }
  # TRUE: fail to reject H0
  # FALSE: reject H0
  return(result)
}

myskewtest <- function(c, print_output=TRUE) {
  # Skew: test for no skew
  skew_rm <- Skew(c, method = 3, conf.level = 0.05, ci.type = "norm", R = 1000)
  # Skewness check for 0: check if 95% confidence interval 
  # ([2] lower, [3] upper) contains zero
  # TRUE: no skewness, property conforms to normality and Gaussian PDF
  # FALSE: has skewness, property does not to conform normality and Gaussian PDF
  # NOTE: does not check left/right skewness
  result <- skew_rm[2] <=0 && skew_rm[3] >=0
  if (print_output) {
    print(skew_rm)
    if(result) sprintf("Skew: no skewness, property conforms to normality and Gaussian PDF")
    else sprintf("Skew: has skewness, property does not conform to normality and Gaussian PDF")
  }
  # TRUE: no skewness
  # FALSE: has skewness
  return(result)
}

mykurttest <- function(c, print_output=TRUE) {
  # Kurtosis: test for no (excess) Kurtosis
  kurt_rm <- Kurt(c, method = 3, conf.level = 0.05, ci.type = "norm", R = 1000)
  # (excess) Kurtosis check for 0: check if 95% confidence interval 
  # ([2] lower, [3] upper) contains zero
  # TRUE: no kurtosis, property conforms to normality and Gaussian PDF
  # FALSE: has skewness, property does not to conform normality and Gaussian PDF
  # NOTE: does not check flat (thin tails)/high (fat tails) Kurtosis
  result <- kurt_rm[2] <=0 && kurt_rm[3] >=0
  if (print_output) {
    print(kurt_rm)
    if(result) sprintf("Kurt: no kurtosis, property conforms to normality and Gaussian PDF")
    else sprintf("Kurt: has (excess) kurtosis, property does not conform to normality and Gaussian PDF")
  }
  # TRUE: no kurtosis
  # FALSE: has kurtosis
  return(result)
}


myboxljungtest <- function(c,print_output=TRUE,lags=30) {
  # Independence: Box-Ljung test
  box_rm <- Box.test(c,lag=lags,type='Ljung') 
  # Box-Ljung: p-value ($p.value) > 0.05
  # TRUE (H0): independent, no autocorrelation -> fail to reject H0
  # FALSE (H1): not independent, autocorrelation present -> reject H0
  result <- box_rm $p.value > 0.05
  if (print_output) {
    print(box_rm)
    if(result) sprintf("Box-Ljung: implies independence over %i lags, no autocorrelation -> fail to reject H0", lags)
    else sprintf("Box-Ljung: implies dependency present over %i lags, autocorrelation present -> reject H0", lags)
  }
  # TRUE: fail to reject H0
  # FALSE: reject H0
  return(result)
}


mykpsstest <- function(c, print_output=TRUE) {
  # Stationarity: KPSS
  kpss_m <- ur.kpss(c,type="tau",lags="short")
  # KPSS: t-stat (@teststat) <= critical value for 5% (@cval[2])
  # TRUE (H0): no unit roots, no linear trend, series is trend stationary -> fail to reject H0
  # FALSE (HA): constains unit roots, has linear trend, series is not trend stationary -> reject H0
  result <- kpss_m@teststat <= kpss_m@cval[2]
  if (print_output) {
    print(summary(kpss_m))
    if(result) sprintf("KPSS: no unit roots, no linear trend, series is trend stationary -> fail to reject H0")
    else sprintf("KPSS: contains unit roots, has linear trend, series not trend stationary -> reject H0")
  }
  # TRUE: fail to reject H0
  # FALSE: reject H0
  return(result)
}


myadftest <- function(c,print_output=TRUE,lags=30) {
  # temporarily suppress warning message: 
  # In adfTest(rm, lags = lags, type = "nc"): p-value smaller than printed p-value
  defaultW <- getOption("warn") 
  options(warn = -1) 
  # Stationarity: ADF
  adf_m <- adfTest(rm,lags=lags,type="nc")
  # ADF: p-value (@test$p.value) >= 0.05
  # TRUE (H0): presence of unit roots, series is not random walk stationary -> fail to reject H0
  # FALSE (HA): contains no unit roots, series is random walk stationary -> reject H0
  result <- adf_m@test$p.value >= 0.05
  if (print_output) {
    print(adf_m)
    if(result) sprintf("ADF: presence of unit roots over %i lags, series is not random walk stationary -> fail to reject H0", lags)
    else sprintf("ADF: contains no unit roots over %i lags, series is random walk stationary -> reject H0", lags)
  }
  # turn warnings back on
  options(warn = defaultW)
  # TRUE: fail to reject H0
  # FALSE: reject H0
  return(result)
}


mymcleodlitest <- function(c, print_output=TRUE) {
  # Constant Variance: McLeod-Li Test
  ml_m <- McLeod.Li.test(c) 
  # McLeod.Li.test: check to see if all lags' p-values ($p.values) are greater than 0.05
  # TRUE (H0): constant variance, homoscedastic -> fail to reject H0
  # FALSE (HA): non-constant variance, heteroscedastic -> reject H0
  result <- length(which(ml_m$p.values >= 5)) > 0
  if (print_output) {
    ml_m
    if(result) sprintf("McLeod-Li: constant variance, homoscedastic -> fail to reject H0")
    else sprintf("McLeod-Li: non-constant variance, heteroscedastic -> reject H0")
  }
  # TRUE: fail to reject H0
  # FALSE: reject H0
  return(result)
}

mybptest <- function(c, print_output=TRUE) {
  # Constant Variance: Breusch-Pagan test
  bp_m <- bptest(lm(c ~ seq(1,length(c))))
  # Breusch-Pagan test: check to see if p-value ($p.value) >= 0.05
  # TRUE (H0): constant variance, homoscedastic -> fail to reject H0
  # FALSE (HA): non-constant variance, possible clustering, heteroscedastic -> reject H0
  result <- bp_m$p.value >= 0.05
  if (print_output) {
    print(bp_m)
    if(result) sprintf("Breusch-Pagan: constant variance, homoscedastic -> fail to reject H0")
    else sprintf("Breusch-Pagan: non-constant variance, possible clustering, heteroscedastic -> reject H0")
  }
  # TRUE: fail to reject H0
  # FALSE: reject H0
  return(result)
}
