##### Packages #####
  library(ggplot2)
  library(reshape2)
  library(locpol)
  library(interp)
  library(stats)
  library(future)
  library(future.apply)
  library(parallel)
  library(tidyverse)

##### source codes #####
source("mean/functions_mean_derivative.R")
# contains the illustration of the mean function and the first derivative in 
# Figure 2.1 

##### Figure 2.1 - Illustration #####

x = seq(0,1, 0.005)
mu_deriv_tibble = tibble(x = x,
                         val = mu_1(x), 
                         fun = "mu'")
mu_tibble = tibble(x = x, 
                   val = mu(x), 
                   fun = "mu")

# estimation of the first derivative 
set.seed(234)
n = 10
obs = seq(0, 1, by = 0.01)

Y = mu(obs) + sapply(1:n, function(x){bm(obs, 0,  sigma = 1)}) + 
  matrix(rnorm(length(obs)*n, 0, 0.3), length(obs), n)

# estimation of derivative and function 
# same bandwidth tho
Y_est3 = locPolSmootherC(x = obs, y = rowMeans(Y), xeval = x, bw = 0.15, deg = 3, EpaK)
Y_est2 = locPolSmootherC(x = obs, y = rowMeans(Y), xeval = x, bw = 0.08, deg = 2, EpaK)

mu_tibble = rbind(mu_tibble, 
                  tibble(x = Y_est2$x, 
                         val = Y_est2$beta0, 
                         fun = "deg = 2"), 
                  tibble(x = Y_est3$x, 
                         val = Y_est3$beta0, 
                         fun = "deg = 3"), 
                  tibble(x = obs, 
                         val = rowMeans(Y), 
                         fun = "obs"))

mu_deriv_tibble = rbind(mu_deriv_tibble, 
                        tibble(x = Y_est2$x, 
                               val = Y_est2$beta1, 
                               fun = "deg = 2"), 
                        tibble(x = Y_est3$x, 
                               val = Y_est3$beta1, 
                               fun = "deg = 3") )

###### Mean Function #####
ggplot() + 
  geom_line(aes(x = x, y = val, col = fun, lty = fun), 
            mu_tibble |> filter(fun != "obs")) + 
  scale_linetype_manual(values = c(1,2,4), 
                        breaks = c("mu", "deg = 2", "deg = 3"),
                        labels = c("deg = 3" = "Cubic", 
                                   'mu'   = expression(mu),
                                   "deg = 2" = "Quad")) + 
  scale_colour_manual(values =  c("black", "green", "red"),
                      breaks = c("mu", "deg = 2", "deg = 3"),
                      labels = c("deg = 3" = "Cubic", 
                                 'mu'   = expression(mu),
                                 "deg = 2" = "Quad")) + 
  geom_point(aes(x = x, y = val), 
             mu_tibble |> filter(fun == "obs"), alpha = .2) + 
  labs(subtitle = "Mean function", y = NULL, x = NULL, colour = NULL, lty = NULL) + 
  deriv_est_theme

ggsave("mean/grafics/mean_est_n10_p101.pdf", device = "pdf", width = 5, height = 3.7, units = "in")

###### First Derivative ######


mu_deriv_tibble |> 
  ggplot() + 
  geom_line(aes(x = x, y = val, col = fun, lty = fun)) + 
  scale_linetype_manual(
    values = c(1, 2, 4),
    breaks = c("mu'", "deg = 2", "deg = 3"),
    labels = list(
      "mu'" = expression(mu * "'"),
      "deg = 2" = "Quad",
      "deg = 3" = "Cubic"
    )
  ) + 
  scale_colour_manual(
    values = c("black", "green", "red"),
    breaks = c("mu'", "deg = 2", "deg = 3"),
    labels = list(
      "mu'" = expression(mu * "'"),
      "deg = 2" = "Quad",
      "deg = 3" = "Cubic"
    )
  )  + 
  labs(subtitle = "First derivative", y = NULL, x = NULL, colour = NULL, lty = NULL) + 
  deriv_est_theme

ggsave("mean/grafics/derivative_est_n10_p101.pdf", device = "pdf", width = 5, height = 3.7, units = "in")
