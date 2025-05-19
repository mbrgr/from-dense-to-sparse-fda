library(plotly)
library(tidyverse)
library(biLocPol) # please install this package from Github first. See "README.md" file for instructions

##### Description #####
# Code for images of the two covariance kernels
# Figure 3.1, 3.2, 3.7

#### results ####
# instead of evaluating the functions all again the results can be loaded with
load("cov/data/illustrations.RData")

##### Processes #####
# Figure of paths of the processes of the second process 
# not contained in the paper
set.seed(134)
p = 40
n =  5
obs = OU(n, t = c(0,(1:p - 0.5)/p), x0 = 0) 
obs_tibble = tibble(Y = as.vector(t(obs)), x = rep(c(0, (1:p - 0.5)/p), n), n = gl(n, p + 1)) 
obs_tibble |> 
  ggplot(aes(x, Y, col = n)) + 
  geom_line(aes(lty = n), size = .6) + 
  geom_point(size = .9) + 
  geom_abline(slope = 0, intercept = 0, lty = 4) + 
  labs(y = NULL, x = NULL) + 
  theme(legend.position = "none", 
        text = element_text(size = 16)) 


##### Figure of (empirical) covariances #####
set.seed(513)
n = 100
p = 40
x.design = (1:p - 1/2)/p
p.eval = 150
h = 0.5

# Target function: Orntein-Uhlenbeck process with following parameters: 
sigma = 2; theta = 3
x = observation_grid(p.eval, comp = "full")[1:p.eval, 1]
Y = FDA_observation(n, x.design, 
                    r.process = OU, process.arg = list(alpha = theta, sigma = sigma, x0 = 0), 
                    eps.arg = list(sd = 0.5))
Z = observation_transformation(Y)

temp = apply(Y, 1, tcrossprod) - as.vector(tcrossprod(colMeans(Y)))     
Z.all = rowSums(temp)/(n-1)
df.all = data.frame(observation_grid(p, comp = "full"), Z.all)

# empirical covariances only 
# Figure not containes in the paper
plot_ly(df.all, x = ~Var1, y = ~Var2, z = ~Z.all, size = .4) |> add_markers()

# evaluate target covariance kernel
cov.ou.eval = matrix(apply(observation_grid(p.eval, comp = "full"), 1, function(x){
  cov_ou(x, sigma, theta)
}), p.eval, p.eval)

# plot empirical covariances with covariance kernel
###### Figure 3.1 (a) & (b) ######
figure31 = plot_ly(df.all, x = ~Var1, y = ~Var2, z = ~Z.all, size = .4) |> 
  add_markers() |> 
  add_surface(x = x, y = x, z = cov.ou.eval, 
              colors = c('#BF382A', '#0C4B8E'), alpha = .3, showscale = FALSE) |> 
  layout() 
figure31b = figure31 |> back_layout()
figure31a = figure31 |> front_layout()

# the following needs the setup to python... otherwise use the export function from the viewer
# front shot
save_image(figure31a, 
           file = "cov/grafics/OU_observation_n100p40theta3sigma2sd05_front.pdf", 
           width = 600, height = 750)
# back shot
save_image(figure31b, 
           file = "cov/grafics/OU_observation_n100p40theta3sigma2sd05_back.pdf", 
           width = 600, height = 750)


##### Estimation #####
# calculate local polynomial weights
p.eval = 150
H = seq(1, 3/p, -0.05)
#h_cv = k_fold_cv(Y, H, h.parallel = T, h.parallel.environment = T)
W = local_polynomial_weights(p, 0.3, p.eval, T) # Calculation of weights (takes a while, data can be loaded instead)
est = eval_weights(W, Z)

# plot with estimation
cs2 = list(c(0, 1), c("lightblue", "darkred"))

##### Figure 3.2 (a) #####
figure32a = plot_ly(df.all, x = ~Var1, y = ~Var2, z = ~Z.all, size = .4) |> 
  add_surface(x = x, y = x, z = cov.ou.eval, alpha = .3, showscale = F) |> 
  add_surface(x = x, y = x, z = est, colorscale = cs2, alpha = .3, showscale = F)
figure32a |> back_layout(x = 2, y = .8, z = .6)
save_image(figure32a |> back_layout(x = 2, y = .8, z = .6), 
           file = "cov/grafics/ou_estimate_m1_h03_sd05.pdf", 
           width = 800, height = 750)

# calculate estimator that does not mirror the results on the diagonal (without diagonal)
#h0_cv =  k_fold_cv(Y, H, m = 0, h.parallel = T, 
#                   h.parallel.environment = T)
W0 = local_polynomial_weights(p, 0.2, p.eval, T, m = 1, grid.type = "without diagonal")
Z0 = observation_transformation(Y, grid.type = "without diagonal")
est_standard = eval_weights(W0, Z0)

##### Figure 3.2 (b) #####
# compare estimate without the diagonal and without mirroring with actual kernel
figure32b = plot_ly(df.all, x = ~Var1, y = ~Var2, z = ~Z.all, size = .4) |> 
  add_surface(x = x, y = x, z = cov.ou.eval, alpha = .3, showscale = F) |> 
  add_surface(x = x, y = x, z = est_standard, colorscale = cs2, alpha = .3,
              showscale = F) |> 
  back_layout(x = 2, y = .8, z = .6)

figure32b
save_image(figure32b, 
           file = "cov/grafics/ou_estimate_m1_h02_full_sd05.pdf", 
           width = 800, height = 750)

#### Comparison with smooth process ####
# Process consisting of two random variables #
set.seed(24)
n = 10
p = 40
z_two_rv(n, p) |> dim()
path_simulation = tibble(paths = z_two_rv(n, p) |> t() |> as.vector(),
                         x     = rep(((1:p) - 1/2)/p, n),
                         count = gl(n, p))
path_simulation |> 
  ggplot(aes(x = x, y = paths, col = count, lty = count)) + 
  geom_line(size = .6)  + 
  labs(y = NULL, x = NULL) + 
  theme(legend.position = "none", 
        text = element_text(size = 16)) 

# figure of paths of second process
# not included in paper
path_simulation |> 
  summarise(.by = x, mean = mean(paths)) |> 
  pull(mean)


# Visualization of covarinace kernel
p.eval = 150
point_grid = observation_grid(p.eval, comp = "full")
x = (1:p.eval - 0.5)/p.eval
cov_val = apply(point_grid, 1, cov_z_2rv)

# estimation with mirroring
n = 100
p = 40
set.seed(24)
Y2 = z_two_rv(n, p)
# h_2rv = k_fold_cv(Y, H, m = 1, h.parallel = T, h.parallel.environment = T)
Z2 = observation_transformation(Y2, grid.type = "less")
W01 = local_polynomial_weights(p, h = 0.1, p.eval = p.eval, m = 1, parallel = T)
est2 = eval_weights(W01, Z2)

cov_val2 = apply(point_grid, 1, cov_z_2rv)
c_val2 = matrix(cov_val2, sqrt(length(cov_val)))

##### Figure 3.7 (a) #####
df_37 = data.frame(observation_grid(p.eval, comp = "full"), 
                   as.vector(c_val2), 
                   as.vector(est2), 
                   as.vector(est_wd))
colnames(df_37) = c("x1", "x2", "c_val2", "est2", "est_wd")

figure37a = plot_ly(df_37) |>
  add_surface(x = x,
              y = x, 
              z = c_val2, 
              size = .4, alpha = .7, 
              showscale = F) |>
  add_surface(x = x,
              y = x, 
              z = est2, 
              size = .4, alpha = .2, colorscale = cs2, 
              showscale = F)  |> back_layout(x = 2.2, y = 1.1, z = 1.1)
figure37a

save_image(figure37a, 
           file = "cov/grafics/2rv_estimate_m1_h01.pdf", 
           width = 700, height = 800)

# estimation without mirroring
Z_wd = observation_transformation(Y2, grid.type = "without diagonal")
# h_2rv = k_fold_cv(Y, H, m = 1, h.parallel = T, 
#                   h.parallel.environment = T)
W_wd = local_polynomial_weights(p, h = 0.1, p.eval = p.eval, m = 1, parallel = T, 
                                grid.type = "without diagonal")
est_wd = eval_weights(W_wd, Z_wd)

##### Figure 3.7 (b) #####
figure37b = plot_ly(df_37) |>
  add_surface(x = x,
              y = x, 
              z = c_val2, 
              size = .4, alpha = .7, 
              showscale = F) |>
  add_surface(x = x,
              y = x, 
              z = est_wd, 
              size = .4, alpha = .7, colorscale = cs2, 
              showscale = F) |>
  back_layout(x = 2.2, y = 1.1, z = 1.1)
figure37b
save_image(figure37b, 
           file = "cov/grafics/2rv_estimate_m1_h01_full.pdf", 
           width = 700, height = 800)


save.image("cov/data/illustrations.RData")

rm(list = setdiff(ls(), 
                  c("figure31a", "figure31b", 
                    "figure32a", "figure32b", 
                    "figure37a", "figure37b")))
save.image("cov/data/plotly_illustration_figures.RData")
