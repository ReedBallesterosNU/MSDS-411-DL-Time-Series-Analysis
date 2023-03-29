
# ProphetOpt.R
# Find optimum values of seasons, changepoints, and holidays
# Revisions by Sean Mason, 2022.05.26
# fit.stats <- apply(A, 1, m.opt, X.train=X.train, X.validation=X.validation, holidays=h)

m.opt <- function(row, X.train, X.validation, holidays) {
    X.train <- data.frame(X.train)
    X.train$cap <- as.numeric(row["capacity"])
    
    m <- prophet(X.train,
                 growth=as.character(row["growth"]),
                 holidays=holidays,
                 seasonality.prior.scale=as.numeric(row["sps"]),
                 changepoint.prior.scale=as.numeric(row["cps"]),
                 holidays.prior.scale=as.numeric(row["hps"]),
                 yearly.seasonality=T)

    future <- make_future_dataframe(m, periods = 184)   # 6 months daily forecasts
    future$cap <- as.numeric(row["capacity"])
    fc <- predict(m, future)
    
    # use forecast package for RMSE and MAPE
    fit.stats <- forecast::accuracy(fc[ymd(fc$ds) %in% X.validation$ds, ]$yhat, X.validation$y)[, c("RMSE", "MAPE")]
    return(fit.stats)
}

