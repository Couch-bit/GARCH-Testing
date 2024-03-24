library(moments)
library(rugarch)
library(tidyverse)

n <- c(100, 250, 500, 1000, 2500, 5000, 10000)
reps <- 5000

garch <- function(mu, omega, alpha, beta, model = 'norm',
                  shape = 0, skew = 0, estimate = 'norm') {
  spec = ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = TRUE,
                                      archm = FALSE, archpow = 1,
                                      arfima = FALSE,
                                      external.regressors = NULL,
                                      archex = FALSE),
                    fixed.pars = list(mu = mu,
                                      omega = omega, alpha1 = alpha,
                                      beta1 = beta, shape = shape, skew = skew),
                    distribution.model = model)
  true_VaR_cross <<- c()
  for (j in n) {
    true_VaR = numeric(reps)
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
      forecast = ugarchforecast(fit, n.ahead = 1)
      true_VaR[k] = quantile(forecast, 0.05)[1] > test
    }
    true_VaR_cross <<- c(true_VaR_cross, mean(true_VaR))
  }
}

calculate_VaR <- function(mu, omega, alpha, beta, dist,
                         shape = 0, skew = 0, model, estimate) {
  garch(mu, omega, alpha, beta, dist, shape, skew, estimate)
  df = data.frame(model = rep(c(model), 7), dist = rep(c(dist), 7),
                  n = n, true_VaR = true_VaR_cross)
  return(df)
}

result <- rbind(calculate_VaR(5, 1, 0.4, 0.5, 'std', 5, model = 1,
                              estimate = 'std'),
                calculate_VaR(5, 1, 0.8, 0.1, 'std', 5, model = 2,
                              estimate = 'std'),
                calculate_VaR(5, 1, 0.1, 0.8, 'std', 5, model = 3,
                              estimate = 'std'),
                calculate_VaR(5, 1, 0.4, 0.5, 'snorm', skew = 4, model = 1,
                              estimate = 'snorm'),
                calculate_VaR(5, 1, 0.8, 0.1, 'snorm', skew = 4, model = 2,
                              estimate = 'snorm'),
                calculate_VaR(5, 1, 0.1, 0.8, 'snorm', skew = 4, model = 3,
                              estimate = 'snorm'))

write.csv(result, file = 'Results/VaRCorrect.csv', row.names = FALSE)
