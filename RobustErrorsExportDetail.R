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

generate_var <- function(stds, mu, omega, alpha, beta, dis,
                         shape = 0, skew = 0, mod) {
  stds = stds %>%
    filter(model == mod, dist == dis)
  omega_sd = stds %>%
    filter(parameter == 'omega')
  omega_sd = omega_sd$std_err
  alpha_sd = stds %>%
    filter(parameter == 'alpha')
  alpha_sd = alpha_sd$std_err
  beta_sd = stds %>%
    filter(parameter == 'beta')
  beta_sd = beta_sd$std_err
  omega_error = numeric(7)
  alpha_error = numeric(7)
  beta_error = numeric(7)
  omega_std = numeric(7)
  alpha_std = numeric(7)
  beta_std = numeric(7)
  garch(mu, omega, alpha, beta, dis, shape, skew)
  for (i in 1:7) {
    omega_error[i] = mean(omega_values[((i-1)*reps + 1):(i*reps)]) - omega_sd[i]
    alpha_error[i] = mean(alpha_values[((i-1)*reps + 1):(i*reps)]) - alpha_sd[i]
    beta_error[i] = mean(beta_values[((i-1)*reps + 1):(i*reps)]) - beta_sd[i]
    omega_std[i] = sd(omega_values[((i-1)*reps + 1):(i*reps)])
    alpha_std[i] = sd(alpha_values[((i-1)*reps + 1):(i*reps)])
    beta_std[i] = sd(beta_values[((i-1)*reps + 1):(i*reps)])
  }
  df = data.frame(model = rep(c(mod), 21), dist = rep(c(dis), 21),
                  n = rep(n, 3), parameter = rep(c('omega', 'alpha', 'beta'),
                                                 each = 7),
                  Me = c(omega_error, alpha_error, beta_error),
                  Mpe = c(omega_error / omega_sd, alpha_error / alpha_sd,
                          beta_error / beta_sd),
                  std = c(omega_std, alpha_std, beta_std))
  return(df)
}

stds <- read.csv('Results/Std.csv')

result <- rbind(generate_var(stds, 5, 1, 0.4, 0.5, 'std', 5, mod = 1),
                generate_var(stds, 5, 1, 0.8, 0.1, 'std', 5, mod = 2),
                generate_var(stds, 5, 1, 0.1, 0.8, 'std', 5, mod = 3),
                generate_var(stds, 5, 1, 0.4, 0.5, 'snorm', skew = 4, mod = 1),
                generate_var(stds, 5, 1, 0.8, 0.1, 'snorm', skew = 4, mod = 2),
                generate_var(stds, 5, 1, 0.1, 0.8, 'snorm', skew = 4, mod = 3))

write.csv(result, file = 'Results/EstDetail.csv', row.names = FALSE)
