library(biLocPol) # please install this package from Github first. See "README.md" file for instructions
library(future.apply)
library(tidyverse)
library(tictoc)
library(plotly)

#### results ####
load("cov/data/results_figure_3_9.RData") # for Illustration in Figure 3.9
# without the weights w and w3
load("cov/data/bandwidth_comparison_derivative_cov_OU.RData") # Figure 3.10


##### Illustration #####

del10_cov_OU = function(t, sigma = 3, theta = 2){
  if(t[1] > t[2]){
    return(sigma^2/2 * (exp(-theta*(t[1] - t[2])) + exp(-theta*(t[1]+t[2])) ))
  }
  #if(t[1] == t[2]){return(sigma^2 * exp(-theta*2*t[1]))}
  else{
    return(sigma^2/2 * (-exp(theta*(t[1] - t[2])) + exp(-theta*(t[1]+t[2])) ))
  }
}

del01_cov_OU = function(t, sigma = 3, theta = 2){
  return(del10_cov_OU(c(t[2], t[1]), sigma, theta))
}

# Parameter OU Process
theta = 2; sigma = 3
sd = 0.75
p.eval = 150
p = 50
n = 100#
h = 0.6
x.design = (1:p - 0.5)/p
#set.seed(99)
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
plot_ly() |> 
  add_surface(x = x, y = x, z = G10_UP, alpha = .8) |> 
  add_surface(x = x, y = x, z = G10_DN, alpha = .8) |> 
  add_surface(x = x, y = x, z = estimate_UP, alpha = .8, colorscale = cs2) |> 
  add_surface(x = x, y = x, z = estimate_DN, alpha = .8, colorscale = cs2) |> 
  layout(scene = list(xaxis = list(title = ""), 
                      yaxis = list(title = ""), 
                      zaxis = list(title = ""))) 
#plot_ly() |> 
#  add_surface(x = x[trim], y = x[trim], z = abs(G10_UP - estimate_UP)[trim, trim]) |> 
#  add_surface(x = x[trim], y = x[trim], z = abs(G10_DN - estimate_DN)[trim, trim]) |> 
#  layout(scene = list(xaxis = list(title = ""), 
#                     yaxis = list(title = ""), 
#                    zaxis = list(title = ""))) 

n 
# parameter second process

del10_cov_z_2rv = function(x){
  4/9 * pi * cos(pi * x[1]) * sin(pi * x[2]) - 10/9 * pi * sin(1.25 * pi * x[1]) * cos(pi * 1.25 * x[2])
}

del01_cov_z_2rv = function(x){
  del10_cov_z_2rv(c(x[2], x[1]))
}

# data 
n = 100
#h2 = 0.3
h3 = 0.2
#w2 = local_polynomial_weights(p, h2, p.eval, F, m = 2, grid.type = "less", del = 1)
w3 = local_polynomial_weights(p, h3, p.eval, F, m = 2, grid.type = "less", del = 1)
set.seed(45)
Y2 = biLocPol::z_2rv(n, p) + matrix(rnorm(n * p, 0, sd), n, p) # n x p
Z2 = observation_transformation(Y2)
estimate_2 = eval_weights(w3, Z2)[,,3]
G10_2 = matrix(apply(observation_grid(p.eval, comp = "full"), 1, del01_cov_z_2rv), p.eval, p.eval)

UP = upper.tri(G10_2, diag = T)
DN = lower.tri(G10_2, diag = F)
estimate_2UP = estimate_2
estimate_2DN = estimate_2
estimate_2DN[UP] = NA
estimate_2UP[DN] = NA

##### Figure 3.9 (b) #####
plot_ly() |> 
  add_surface(x = x, y = x, z = G10_2, alpha = .8) |> 
  add_surface(x = x, y = x, z = estimate_2UP, alpha = .8, colorscale = cs2) |> 
  add_surface(x = x, y = x, z = estimate_2DN, alpha = .8, colorscale = cs2) |> 
  layout(scene = list(xaxis = list(title = ""), 
                      yaxis = list(title = ""), 
                      zaxis = list(title = ""))) 

rm(w, w3)
save.image("cov/data/results_figure_3_9.RData")


##### Bandwidth Comparison for Covariance derivative estimation #####
N = 500 #N = 1000
#n.seq = c(50) #
n.seq = c(50, 100, 200, 400)
#p.seq = c(50) #
p.seq = c(15, 25, 50, 75)

p.eval= 50 #p.eval = 100
H = lapply(1:length(p.seq), function(l){seq(1, max(5/p.seq[l], 0.1), -0.05)})

# Parameter OU Process
theta = 2; sigma = 3
# Standard deviation for additional errors
sd = 0.25
set.seed(87)
bw_comparison = list()

cl = parallel::makeCluster(parallel::detectCores( ) - 1)
future::plan(future::multisession)

for(l in 1:length(p.seq)){
  tic()
  bw_comparison[[l]] = matrix(t(future_sapply(1:length(H[[l]]), FUN = function(k)
  {
    bandwidth_evaluation_derivative_OU(H[[l]][k], p.seq[l], p.eval, n.seq, 
                                       N, boundary_correction = T)
  },
  future.seed = T)),ncol = 4)
  cat("p =", p.seq[l], "done.")
  toc()
}

parallel::stopCluster(cl)

save.image("cov/data/bandwidth_comparison_derivative_cov_OU.RData")

bw_comparison_tibble = Reduce(rbind, bw_comparison)
colnames(bw_comparison_tibble) = c("n", "p", "h", "sup.err")
bw_comparison_tibble = bw_comparison_tibble %>% 
  as_tibble() %>% 
  mutate(p = as.factor(p))


##### Figure 3.10
ggplot(bw_comparison_tibble) + 
  geom_point(aes(x = h, y = sup.err, col = p)) + 
  lims(y = c(0, 10.4)) +
  facet_wrap(n ~., nrow = 1) 

ggsave("grafics/bw_comp_cov_deriv_OU.png", device = "png", unit = "in", width = 10, height = 5)

bw_comparison_tibble |> 
  group_by(n, p) |> 
  slice_min(sup.err)

##### Bandwidth Evaluation for the smooth process #####

bandwidth_evaluation_derivative_2rv = function(h, p, p.eval, n.seq, N, eps.sd = 0.75, 
                                               boundary_correction = F, correction = 0.1){
  if(!boundary_correction){
    correction = 0
  }
  x.design = (1:p - 0.5)/p
  x.eval.design = (1:p.eval - 0.5)/p.eval
  s = length(h)
  x.eval.grid = observation_grid(p.eval, comp = "full")
  G10 = matrix(apply(x.eval.grid, 1, del10_cov_z_2rv), p.eval, p.eval)
  
  sup.err10 = numeric(length(n.seq))
  w = local_polynomial_weights(p, h, p.eval, F, m = 2, grid.type = "less", del = 1)
  for(k in 1:length(n.seq)){
    n = n.seq[k]
    sup.err10[k] =  mean(replicate(N, {
      Y = biLocPol::z_2rv(n, p) + matrix(rnorm(n * p, 0, eps.sd), n, p) # n x p
      Z = observation_transformation(Y)
      trim = which(x.eval.design > correction & x.eval.design < (1-correction))
      estimate = eval_weights(w, Z)[,,3]
      max(abs(estimate - G10)[trim, trim])
    }))
  }
  rm(w)
  matrix(c(n.seq, rep(c(p, h), each = length(n.seq)), sup.err10), ncol =  4)
}


##### Bandwidth Comparison for Covariance derivative estimation #####
N = 500 #N = 1000
#n.seq = c(50) #
n.seq = c(50, 100, 200, 400)
#p.seq = c(50) #
p.seq = c(15, 25, 50, 75)

p.eval= 50 #p.eval = 100
H = lapply(1:length(p.seq), function(l){seq(1, max(5/p.seq[l], 0.1), -0.05)})

# Standard deviation for additional errors
sd = 0.75
set.seed(137)
bw_comparison = list()

cl = parallel::makeCluster(parallel::detectCores( ) - 1)
future::plan(future::multisession)

for(l in 1:length(p.seq)){
  tic()
  bw_comparison[[l]] = matrix(t(future_sapply(1:length(H[[l]]), FUN = function(k)
  {
    bandwidth_evaluation_derivative_2rv(H[[l]][k], p.seq[l], p.eval, n.seq, N, 
                                        boundary_correction = T)
  },
  future.seed = T)),ncol = 4)
  cat("p =", p.seq[l], "done.")
  toc()
}

parallel::stopCluster(cl)

save.image("data/bandwidth_comparison_derivative_cov_2rv.RData")
load("data/bandwidth_comparison_derivative_cov_2rv.RData")

bw_comparison_tibble = Reduce(rbind, bw_comparison)
colnames(bw_comparison_tibble) = c("n", "p", "h", "sup.err")
bw_comparison_tibble = bw_comparison_tibble %>% 
  as_tibble() %>% 
  mutate(p = as.factor(p))

ggplot(bw_comparison_tibble) + 
  geom_point(aes(x = h, y = sup.err, col = p, pch = p)) + 
  lims(y = c(0, 7.5)) +
  facet_wrap(n ~., nrow = 1)

ggsave("grafics/cov_deriv_2rv_bw_comp.png", device = "png", units = "in", width = 10, height = 5)

