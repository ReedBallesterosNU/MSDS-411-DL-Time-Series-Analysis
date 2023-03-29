
# By Mark Heiple


parameterTest <- function(model)  # model <- m; diff.order <- 0
{
  se = sqrt(diag(model$var.coef))
  t = model$coef/se
  signif(t,3)

  len = length(resid(model))   # set length of time series
  # model$arma gives number of model coefficients plus difference order. If ts model
  # from a differenced series, then diff.order gives this to adust model$arma
  df = len - sum(model$arma) + ifelse(model$arma[5]<2,1,0)
  pval = (1-pt(abs(t),df=df))*2
  signif(pval,3)

  rc = list()
  rc$t = t
  rc$pval_t = pval

  pval = (1-pnorm(abs(t)))*2
  signif(pval,3)
  rc$pval_z = pval

  rcdf = data.frame(rc)
  pr = rep(' ',length(rc$pval_z))
  pr[which(rc$pval_z<=.1)] = '.'
  pr[which(rc$pval_z<=.05)] = '*'
  pr[which(rc$pval_z<=.01)] = '**'
  pr[which(round(rc$pval_z,3)==0)] = '***'
  rcdf$`Pr(>|t|)` = pr

  return(rcdf)
}
