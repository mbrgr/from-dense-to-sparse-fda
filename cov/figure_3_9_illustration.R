library(biLocPol) # please install this package from Github first. See "README.md" file for instructions
library(future.apply)
library(tidyverse)
library(tictoc)
library(plotly)

#### results ####
load("cov/data/results_figure_3_9a.RData") # for Illustration in Figure 3.9a
load("cov/data/results_figure_3_9b.RData") # for Illustration in Figure 3.9b
# without the weights w and w3

source("cov/functions.R")

##### Illustration #####


# Parameter OU Process
theta = 2; sigma = 3
sd = 0.75
p.eval = 150
p = 50
n = 100#
h = 0.6
x.design = (1:p - 0.5)/p
set.seed(99)
Y = biLocPol::OU(n, x.design, sigma = sigma, alpha = theta) + matrix(rnorm(n * p, 0, sd), n, p) # n x p
Z = observation_transformation(Y)
w = local_polynomial_weights(p, h, p.eval, F, m = 2, grid.type = "less", del = 1)
estimate = eval_weights(w, Z)[,,3]
x = observation_grid(p.eval, comp = "full")[1:p.eval, 1]
G10 = matrix(apply(observation_grid(p.eval, comp = "full"), 1, del10_cov_OU), p.eval, p.eval)

UP = upper.tri(G10, diag = T)
G10_UP = G10
G10_UP[!UP] = NA
G10_DN = G10
G10_DN[UP] = NA
estimate_UP = estimate
estimate_DN = estimate
estimate_DN[UP] = NA
estimate_UP[!UP] = NA
trim = which(x.design > 0.2 & x.design < 0.8)

cs2 = lisG10_UPcs2 = list(c(0, 1), c("lightblue", "darkred"))

##### Figure 3.9 (a) #####
figure39a = plot_ly() |> 
  add_surface(x = x, y = x, z = G10_UP, 
              alpha = .8, showscale = F) |> 
  add_surface(x = x, y = x, z = G10_DN, 
              alpha = .8, showscale = F) |> 
  add_surface(x = x, y = x, z = estimate_UP, alpha = .8, colorscale = cs2,
              showscale = F) |> 
  add_surface(x = x, y = x, z = estimate_DN, 
              alpha = .8, colorscale = cs2, showscale = F) |>
  front_layout(x = -.9, y = -2.3, z = .7)
figure39a
save_image(figure39a, 
           file = "cov/grafics/OU_G10.pdf", 
           width = 600, height = 750)

rm(w)
save.image("cov/data/results_figure_3_9a.RData")

# parameter second process

del10_cov_z_2rv = function(x){
  4/9 * pi * cos(pi * x[1]) * sin(pi * x[2]) - 10/9 * pi * sin(1.25 * pi * x[1]) * cos(pi * 1.25 * x[2])
}

del01_cov_z_2rv = function(x){
  del10_cov_z_2rv(c(x[2], x[1]))
}

# data 
n = 100
h2 = 0.3
#h3 = 0.2
w2 = local_polynomial_weights(p, h2, p.eval, T, m = 2, grid.type = "less", del = 1)
#w3 = local_polynomial_weights(p, h3, p.eval, F, m = 2, grid.type = "less", del = 1)
set.seed(45)
Y2 = biLocPol::z_2rv(n, p) + matrix(rnorm(n * p, 0, sd), n, p) # n x p
Z2 = observation_transformation(Y2)
estimate_2 = eval_weights(w2, Z2)[,,3]
G10_2 = matrix(apply(observation_grid(p.eval, comp = "full"), 1, del01_cov_z_2rv), p.eval, p.eval)

UP = upper.tri(G10_2, diag = T)
DN = lower.tri(G10_2, diag = F)
estimate_2UP = estimate_2
estimate_2DN = estimate_2
estimate_2DN[UP] = NA
estimate_2UP[DN] = NA

##### Figure 3.9 (b) #####
figure39b = plot_ly() |> 
  add_surface(x = x, y = x, z = G10_2, alpha = .8, showscale = F) |> 
  add_surface(x = x, y = x, z = estimate_2UP, alpha = .8, colorscale = cs2, showscale = F) |> 
  add_surface(x = x, y = x, z = estimate_2DN, alpha = .8, colorscale = cs2, showscale = F) |>
  front_layout(x = -.9, y = -2.3, z = .7)
figure39b

save_image(figure39a, 
           file = "cov/grafics/2rv_G10.pdf", 
           width = 600, height = 750)
rm(w)
rm(w2)
save.image("cov/data/results_figure_3_9b.RData")

rm(list = setdiff(ls(), 
                  c("figure39a", "figure39b")))
save.image("cov/data/figure_3_9.RData")
