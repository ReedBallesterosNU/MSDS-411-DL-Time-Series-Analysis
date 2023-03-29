
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
#   tests to help understand the results in a more layman's perspective
# - (most) tests below are tested against their respective null hypothesis (H0)
#     TRUE: FAIL to reject H0
#     FALSE: reject H0
#   check each specific function for details



# ** DISCLAIMER **
# - posted as-is without warranty
# - use at your discretion
# - may contain bugs/errors/changes
# - far from perfect; to be used as own learning tool/reference



# version history:
# v0.01:
#   - made parameter names more explanatory
#   - document inputs and returns for each test
#   - mymcleodlitest - fixed:
#     result <- length(which(ml_m$p.values <= 0.05)) == 0
# v0.02:
#   - emphasized NO/NOT/FAIL in documentation for visibility
#   - rewrote documentation for each test
#   - removed TRUE/FALSE comments for each test:
#       refer to print statement for each test based on (result)
#   - myadftest - added drift explanation based on prof's updates
#       (link accessible only to class):
#     https://canvas.northwestern.edu/courses/177824/pages/gaussian-based-time-series-models-ts2?module_item_id=2405635
# v0.03:
#   - emphasized *NO*/*NOT*/*FAIL* in documentation for visibility
#   - updated KPSS/ADF information
#   - myskewtest: added left/right Skewness output 
#   - mykurttest: added flat/tall (excess) Kurtosis output


# ----------------------------------



# input:
#   data: raw data or model residuals
#   print_output: output the test and the explanatory result against H0 (i.e. layman terms) - default TRUE
# returns:
#  TRUE: FAIL to reject H0
#  FALSE: reject H0
myttest <- function(data, print_output=TRUE) {
  # Mean 0: T-test
  # t.test() test: check if 95% confidence interval ($conf.int[1] lower, $conf.int[2] upper) 
  #   contains 0 and p-value ($p.value) >= 0.05 (H0)
  ttest_data <- t.test(data)
  result <- ttest_data$conf.int[1] <=0 && ttest_data$conf.int[2] >=0 && ttest_data$p.value >= 0.05
  if (print_output) {
    print(ttest_data)
    if(result)
      cat("T-Test: mean is statistically zero, linear trend *REMOVED* -> \n*FAIL* to reject H0\n")
    else
      cat("T-Test: mean *NOT* statistically zero, linear trend present -> \nreject H0\n")
  }
  # TRUE: FAIL to reject H0
  # FALSE: reject H0
  return(result)
}


# input:
#   data: raw data or model residuals
#   print_output: output the test and the explanatory result (i.e. layman terms) - default TRUE
# returns:
#  TRUE: NO skewness
#  FALSE: has skewness
myskewtest <- function(data, print_output=TRUE) {
  # Skew: test for NO skew
  # Skew() test for 0: check if 95% confidence interval ([2] lower, [3] upper) contains zero
  skew_data <- Skew(data, method = 3, conf.level = 0.05, ci.type = "norm", R = 1000)
  result <- skew_data[2] <=0 && skew_data[3] >=0
  if (print_output) {
    print(skew_data)
    if(result)
      cat("Skew: *NO* skewness, \nproperty conforms to normality and Gaussian PDF\n")
    else
    {
      skewType <- "LEFT"
      if (skew_data[1] > 0) skewType <- "RIGHT"
      cat(sprintf("Skew: has *%s* skewness, \nproperty does *NOT* conform to normality and Gaussian PDF\n", skewType))
    }
  }
  # TRUE: NO skewness
  # FALSE: has skewness
  return(result)
}


# input:
#   data: raw data or model residuals
#   print_output: output the test and the explanatory result (i.e. layman terms) - default TRUE
# returns:
#  TRUE: NO (excess) Kurtosis
#  FALSE: has (excess) Kurtosis
mykurttest <- function(data, print_output=TRUE) {
  # Kurtosis: test for NO (excess) Kurtosis
  # Kurt() test for 0: check if 95% confidence interval ([2] lower, [3] upper) contains zero
  kurt_data <- Kurt(data, method = 3, conf.level = 0.05, ci.type = "norm", R = 1000)
  result <- kurt_data[2] <=0 && kurt_data[3] >=0
  if (print_output) {
    print(kurt_data)
    if(result)
      cat("Kurt: *NO* (excess) kurtosis, \nproperty conforms to normality and Gaussian PDF\n")
    else {
      kurtType <- "FLAT thin-tailed"
      if (kurt_data[1] > 0) kurtType <- "TALL thick-tailed"
      cat(sprintf("Kurt: has *%s* (excess) kurtosis, \nproperty does *NOT* conform to normality and Gaussian PDF\n", kurtType))
    }
  }
  # TRUE: NO kurtosis
  # FALSE: has kurtosis
  return(result)
}


# input:
#   data: raw data or model residuals
#   print_output: output the test and the explanatory result against H0 (i.e. layman terms) - default TRUE
#   lags: number of lags - default 30
# returns:
#  TRUE: FAIL to reject H0
#  FALSE: reject H0
myboxljungtest <- function(data,print_output=TRUE,lags=30) {
  # Independence: Box-Ljung test
  # Box.test(type='Ljung') test: check to see if p-value ($p.value) >= 0.05 (H0)
  bl_data <- Box.test(data,lag=lags,type='Ljung') 
  result <- bl_data$p.value >= 0.05
  if (print_output) {
    print(bl_data)
    if(result)
      cat(sprintf("Box-Ljung: implies independence over %i lags, \n*NO* autocorrelation -> *FAIL* to reject H0\n", lags))
    else
      cat(sprintf("Box-Ljung: implies dependency present over %i lags, \nautocorrelation present -> reject H0\n", lags))
  }
  # TRUE: FAIL to reject H0
  # FALSE: reject H0
  return(result)
}


# input:
#   data: raw data or model residuals
#   print_output: output the test and the explanatory result against H0 (i.e. layman terms) - default TRUE
# returns:
#  TRUE: FAIL to reject H0
#  FALSE: reject H0
mykpsstest <- function(data, print_output=TRUE) {
  # Stationarity: KPSS Test
  # ur.kpss() test: check if t-stat (@teststat) < critical value for 5% (@cval[2]) (H0)
  kpss_data <- ur.kpss(data,type="tau",lags="short")
  result <- kpss_data@teststat < kpss_data@cval[2]
  if (print_output) {
    print(summary(kpss_data))
    if(result)
      cat("KPSS: *NO* unit roots, *NO* linear trend, slope zero, \nseries is trend stationary -> *FAIL* to reject H0\n")
    else
      cat("KPSS: contains unit roots, linear trend present, slope *NOT* zero, \nseries *NOT* trend stationary -> reject H0\n")
  }
  # TRUE: FAIL to reject H0
  # FALSE: reject H0
  return(result)
}


# input:
#   data: raw data or model residuals
#   print_output: output the test and the explanatory result against H0 (i.e. layman terms) - default TRUE
#   lags: number of lags - default 30
# returns:
#  TRUE: FAIL to reject H0
#  FALSE: reject H0
myadftest <- function(data,print_output=TRUE,lags=30) {
  # temporarily suppress warning message: 
  # In adfTest(rm, lags = lags, type = "nc"): p-value smaller than printed p-value
  defaultW <- getOption("warn") 
  options(warn = -1) 
  # Stationarity: ADF Test
  # adfTest() test: check to see if p-value (@test$p.value) >= 0.05 (H0)
  adf_data <- adfTest(data,lags=lags,type="nc")
  result <- adf_data@test$p.value >= 0.05
  if (print_output) {
    print(adf_data)
    if(result)
      cat(sprintf("ADF: presence of unit roots over %i lags, indicates mean drift, \nbusiness cycles present, series is *NOT* stationary -> *FAIL* to reject H0\n", lags))
    else
      cat(sprintf("ADF: contains *NO* unit roots over %i lags, indicates *NO* mean drift, \nbusiness cycles *NOT* present, series is stationary -> reject H0\n", lags))
  }
  # turn warnings back on
  options(warn = defaultW)
  # TRUE: FAIL to reject H0
  # FALSE: reject H0
  return(result)
}


# input:
#   model: the Arima() model
#   print_output: output the test and the explanatory result against H0 (i.e. layman terms) - default TRUE
# returns:
#  TRUE: FAIL to reject H0
#  FALSE: reject H0
mymcleodlitest <- function(model, print_output=TRUE) {
  # Constant Variance: McLeod-Li Test
  # McLeod.Li.test() test: check to see if all lags' p-values ($p.values) are greater than 0.05 (H0)
  ml_m <- McLeod.Li.test(model)
  lagsUnder <- which(ml_m$p.values < 0.05)
  result <- length(lagsUnder) == 0
  if (print_output) {
    if(result) {
      cat("McLeod-Li: constant variance, homoscedastic -> *FAIL* to reject H0\n")
      cat(sprintf("McLeod-Li: Lags >= 0.05:\n  (none)\n"))
    }
    else{
      cat("McLeod-Li: *NON*-constant variance, heteroscedastic -> reject H0\n")
      lags <- paste(lagsUnder,collapse = ',')
      cat(sprintf("McLeod-Li: Lags < 0.05:\n%s\n",lags))
    }
  }
  # TRUE: FAIL to reject H0
  # FALSE: reject H0
  return(result)
}


# input:
#   data: raw data or model residuals
#   print_output: output the test and the explanatory result against H0 (i.e. layman terms) - default TRUE
# returns:
#  TRUE: FAIL to reject H0
#  FALSE: reject H0
mybptest <- function(data, print_output=TRUE) {
  # Constant Variance: Breusch-Pagan Test
  # bptest() test: check to see if p-value ($p.value) >= 0.05 (H0)
  bp_data <- bptest(lm(data ~ seq(1,length(data))))
  result <- bp_data$p.value >= 0.05
  if (print_output) {
    print(bp_data)
    if(result)
      cat("Breusch-Pagan: constant variance, homoscedastic ->\n*FAIL* to reject H0\n")
    else
      cat("Breusch-Pagan: *NON*-constant variance, possible clustering, \nheteroscedastic -> reject H0\n")
  }
  # TRUE: FAIL to reject H0
  # FALSE: reject H0
  return(result)
}
