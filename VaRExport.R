library(moments)
library(rugarch)
library(tidyverse)

n <- c(100, 250, 500, 1000, 2500, 5000, 10000)
reps <- 5000

standard_historic_quantile <- function(data) {
  return(quantile(data, probs = c(0.05))[1])
}

modified_quantile <- function(data, gauss_VaR) {
  skew = skewness(data)
  kurt = kurtosis(data)
  VaR = (gauss_VaR 
         + (gauss_VaR^2 - 1) * skew/6
         + (gauss_VaR^3 - 3 * gauss_VaR) * (kurt - 3)/24
         - (2 * gauss_VaR^3 - 5 * gauss_VaR) * skew^2/36)
  return(VaR)
}

garch <- function(mu, omega, alpha, beta, model = 'norm',
                  shape = 0, skew = 0, estimate = 'norm') {
  gauss_VaR = qnorm(0.05)
  spec = ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = TRUE,
                                      archm = FALSE, archpow = 1,
                                      arfima = FALSE,
                                      external.regressors = NULL,
                                      archex = FALSE),
                    fixed.pars = list(mu = mu,
                                      omega = omega, alpha1 = alpha,
                                      beta1 = beta, shape = shape, skew = skew),
                    distribution.model = model)
  gauss_cross <<- c()
  historic_cross <<- c()
  modified_cross <<- c()
  for (j in n) {
    gaussian = numeric(reps)
    historic = numeric(reps)
    modified = numeric(reps)
    for (k in 1:reps) {
      fit = NULL
      while (is.null(fit) || convergence(fit) != 0) {
        data = fitted(ugarchpath(spec, j + 1))
        data_fit = data[1:j]
        test = data[j + 1]
        specfit = ugarchspec(mean.model = list(armaOrder = c(0, 0),
                                               include.mean = TRUE,
                                               archm = FALSE, 
                                               archpow = 1, arfima = FALSE,
                                               external.regressors = NULL,
                                               archex = FALSE),
                             fixed.pars = list(shape = shape, skew = skew),
                             distribution.model = estimate)
        fit = ugarchfit(spec = specfit, data = data_fit, solver = 'hybrid')
      }
      res = residuals(fit, standardize = TRUE)
      forecast = ugarchforecast(fit, n.ahead = 1)
      sig = sigma(forecast)[1]
      value = fitted(forecast)[1]
      gaussian[k] = value + sig*gauss_VaR > test
      historic[k] = value + sig*standard_historic_quantile(res) > test
      modified[k] = value + sig*modified_quantile(res, gauss_VaR) > test
    }
    gauss_cross <<- c(gauss_cross, mean(gaussian))
    historic_cross <<- c(historic_cross, mean(historic))
    modified_cross <<- c(modified_cross, mean(modified))
  }
}

calculate_VaR <- function(mu, omega, alpha, beta, dist,
                         shape = 0, skew = 0, model) {
  garch(mu, omega, alpha, beta, dist, shape, skew)
  df = data.frame(model = rep(c(model), 7), dist = rep(c(dist), 7),
                  n = n, gaussian = gauss_cross, historic = historic_cross,
                  modified = modified_cross)
  return(df)
}

result <- rbind(calculate_VaR(5, 1, 0.4, 0.5, 'std', 5, model = 1),
                calculate_VaR(5, 1, 0.8, 0.1, 'std', 5, model = 2),
                calculate_VaR(5, 1, 0.1, 0.8, 'std', 5, model = 3),
                calculate_VaR(5, 1, 0.4, 0.5, 'snorm', skew = 4, model = 1),
                calculate_VaR(5, 1, 0.8, 0.1, 'snorm', skew = 4, model = 2),
                calculate_VaR(5, 1, 0.1, 0.8, 'snorm', skew = 4, model = 3))

write.csv(result, file = 'Results/VaR.csv', row.names = FALSE)
