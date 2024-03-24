library(rugarch)
library(tidyverse)

n <- c(100, 250, 500, 1000, 2500, 5000, 10000)
reps <- 100

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
    omega_sd = numeric(reps)
    alpha_sd = numeric(reps)
    beta_sd = numeric(reps)
    for (k in 1:reps) {
      fit = NULL
      while (is.null(fit) || 
             convergence(fit) != 0 ||
             sum(is.na(vcov(fit, robust = TRUE))) != 0) {
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
      omega_sd[k] = sqrt(vcov(fit, robust = TRUE)[2, 2])
      alpha_sd[k] = sqrt(vcov(fit, robust = TRUE)[3, 3])
      beta_sd[k] = sqrt(vcov(fit, robust = TRUE)[4, 4])
    }
    omega_values <<- c(omega_values, omega_sd)
    alpha_values <<- c(alpha_values, alpha_sd)
    beta_values <<- c(beta_values, beta_sd)
  }
}

generate_var <- function(mu, omega, alpha, beta, dist,
                         shape = 0, skew = 0, model) {
  
  garch(mu, omega, alpha, beta, dist, shape, skew)
  df = data.frame(model = rep(c(model), reps * 21),
                  dist = rep(c(dist), reps * 21),
                  n = rep(n, 3, each = reps),
                  parameter = rep(c('omega', 'alpha', 'beta'),
                                                 each = reps * 7),
                  std = c(omega_values, alpha_values, beta_values))
  return(df)
}

result <- rbind(generate_var(5, 1, 0.4, 0.5, 'std', 5, model = 1),
                generate_var(5, 1, 0.8, 0.1, 'std', 5, model = 2),
                generate_var(5, 1, 0.1, 0.8, 'std', 5, model = 3),
                generate_var(5, 1, 0.4, 0.5, 'snorm', skew = 4, model = 1),
                generate_var(5, 1, 0.8, 0.1, 'snorm', skew = 4, model = 2),
                generate_var(5, 1, 0.1, 0.8, 'snorm', skew = 4, model = 3))

write.csv(result, file = 'Results/Est.csv', row.names = FALSE)
