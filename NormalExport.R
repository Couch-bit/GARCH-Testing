library(rugarch)
library(tidyverse)

n <- c(100, 250, 500, 1000, 2500, 5000, 10000)
reps <- 3000

garch <- function(mu, omega, alpha, beta, model = 'norm',
                  shape = 0, skew = 0, estimate = 'norm') {
  omega_values <<- c()
  alpha_values <<- c()
  beta_values <<- c()
  spec = ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = TRUE,
                                      archm = FALSE, archpow = 1,
                                      arfima = FALSE,
                                      external.regressors = NULL,
                                      archex = FALSE),
                    fixed.pars = list(mu = mu,
                                      omega = omega, alpha1 = alpha,
                                      beta1 = beta, shape = shape, skew = skew),
                    distribution.model = model)
  for (j in n) {
    omega = numeric(reps)
    alpha = numeric(reps)
    beta = numeric(reps)
    for (k in 1:reps) {
      fit = NULL
      while (is.null(fit) || convergence(fit) != 0) {
        data = fitted(ugarchpath(spec, j))
        specfit = ugarchspec(mean.model = list(armaOrder = c(0, 0),
                                               include.mean = TRUE,
                                               archm = FALSE, 
                                               archpow = 1, arfima = FALSE,
                                               external.regressors = NULL,
                                               archex = FALSE),
                             fixed.pars = list(shape = shape, skew = skew),
                             distribution.model = estimate)
        fit = ugarchfit(spec = specfit, data = data, solver = 'hybrid')
      }
      omega[k] = coef(fit)[2]
      alpha[k] = coef(fit)[3]
      beta[k] = coef(fit)[4]
    }
    omega_values <<- c(omega_values, omega)
    alpha_values <<- c(alpha_values, alpha)
    beta_values <<- c(beta_values, beta)
  }
}

calculate_error_and_variance <- function(mu, omega, alpha, beta, model) {
  omega_error = numeric(7)
  alpha_error = numeric(7)
  beta_error = numeric(7)
  omega_std = numeric(7)
  alpha_std = numeric(7)
  beta_std = numeric(7)
  garch(mu, omega, alpha, beta)
  for (i in 1:7) {
    omega_error[i] = mean(omega_values[((i-1)*reps + 1):(i*reps)]) - omega
    alpha_error[i] = mean(alpha_values[((i-1)*reps + 1):(i*reps)]) - alpha
    beta_error[i] = mean(beta_values[((i-1)*reps + 1):(i*reps)]) - beta
    omega_std[i] = sd(omega_values[((i-1)*reps + 1):(i*reps)])
    alpha_std[i] = sd(alpha_values[((i-1)*reps + 1):(i*reps)])
    beta_std[i] = sd(beta_values[((i-1)*reps + 1):(i*reps)])
  }
  df = data.frame(model = rep(c(model), 21),
                  n = rep(n, 3), parameter = rep(c('omega', 'alpha', 'beta'),
                                                 each = 7),
                  Me = c(omega_error, alpha_error, beta_error),
                  Mpe = c(omega_error / omega, alpha_error / alpha,
                          beta_error / beta),
                  std = c(omega_std, alpha_std, beta_std))
  return(df)
}

result <- rbind(calculate_error_and_variance(5, 1, 0.4, 0.5, 1),
                calculate_error_and_variance(5, 1, 0.8, 0.1, 2),
                calculate_error_and_variance(5, 1, 0.1, 0.8, 3))

write.csv(result, file = 'Results/ErrorNorm.csv', row.names = FALSE)
