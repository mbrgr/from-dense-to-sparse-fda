#### packages, code ####
library(ggplot2)
library(tidyverse)
library(biLocPol)
library(future.apply)
library(tictoc)

my_theme = theme_grey(base_size = 15) + 
  theme(plot.title = element_text(size = 14))

#### results ####
load("cov/data/estimator_comparison.RData")

#  Figures 3.6
source("cov/functions.r")

#### calculations ####
# Examine a final sample error for mirrored process

N = 1000
n.seq = c(100, 300)
p.seq = c( 50,  75)

p.eval = 75
H = lapply(1:length(p.seq), function(l){seq(1, 3/p.seq[l], -0.05)})


##### OU #####

# Parameter OU Process
theta = 2; sigma = 3
# Standard deviation for additional errors
sd = 0.75

##### m = 0 #####
m = 0

bw_comparison_OU_m0 = list()
bw_comparison_OU_m0_full = list()

cl = parallel::makeCluster(parallel::detectCores( ) - 1)
future::plan(future::cluster)

for(l in 1:length(p.seq)){
  tic()
  bw_comparison_OU_m0[[l]]      = matrix(t(future_sapply(1:length(H[[l]]), FUN = function(k)
  {
    bandwidth_evaluation(H[[l]][k], p.seq[l], p.eval, n.seq, N, 
                         cov_ou, list(theta = theta, sigma = sigma), 
                         OU, list(alpha = theta, sigma = sigma, x0 = 0), 
                         eps.arg = list(sd = sd), m = m)
  },
  future.seed = T)),ncol = 4)
  cat("p =", p.seq[l], "done.")
  
  bw_comparison_OU_m0_full[[l]] = matrix(t(future_sapply(1:length(H[[l]]), FUN = function(k)
  {
    bandwidth_evaluation(H[[l]][k], p.seq[l], p.eval, n.seq, N, 
                         cov_ou, list(theta = theta, sigma = sigma), 
                         OU, list(alpha = theta, sigma = sigma, x0 = 0), 
                         eps.arg = list(sd = sd), m = m, grid.type = "without diagonal")
  },
  future.seed = T)),ncol = 4)
  
  toc()
}

parallel::stopCluster(cl)


##### m = 1 #####
m = 1

bw_comparison_OU_m1 = list()
bw_comparison_OU_m1_full = list()

cl = parallel::makeCluster(parallel::detectCores( ) - 1)
future::plan(future::cluster)

for(l in 1:length(p.seq)){
  tic()
  bw_comparison_OU_m1[[l]]      = matrix(t(future_sapply(1:length(H[[l]]), FUN = function(k)
  {
    bandwidth_evaluation(H[[l]][k], p.seq[l], p.eval, n.seq, N, 
                         cov_ou, list(theta = theta, sigma = sigma), 
                         OU, list(alpha = theta, sigma = sigma, x0 = 0), 
                         eps.arg = list(sd = sd), m = m)
  },
  future.seed = T)),ncol = 4)
  cat("p =", p.seq[l], "done.")
  
  bw_comparison_OU_m1_full[[l]] = matrix(t(future_sapply(1:length(H[[l]]), FUN = function(k)
  {
    bandwidth_evaluation(H[[l]][k], p.seq[l], p.eval, n.seq, N, 
                         cov_ou, list(theta = theta, sigma = sigma), 
                         OU, list(alpha = theta, sigma = sigma, x0 = 0), 
                         eps.arg = list(sd = sd), m = m, grid.type = "without diagonal")
  },
  future.seed = T)),ncol = 4)
  
  toc()
}

parallel::stopCluster(cl)

##### m = 2 #####
m = 2

H2 = lapply(1:length(p.seq), function(l){seq(1, 4/p.seq[l], -0.05)})
bw_comparison_OU_m2 = list()
bw_comparison_OU_m2_full = list()

cl = parallel::makeCluster(parallel::detectCores( ) - 1)
future::plan(future::cluster)

for(l in 1:length(p.seq)){
  tic()
  bw_comparison_OU_m2[[l]]      = matrix(t(future_sapply(1:length(H2[[l]]), FUN = function(k)
  {
    bandwidth_evaluation(H2[[l]][k], p.seq[l], p.eval, n.seq, N, 
                         cov_ou, list(theta = theta, sigma = sigma), 
                         OU, list(alpha = theta, sigma = sigma, x0 = 0), 
                         eps.arg = list(sd = sd), m = m)
  },
  future.seed = T)),ncol = 4)
  cat("p =", p.seq[l], "done.")
  
  bw_comparison_OU_m2_full[[l]] = matrix(t(future_sapply(1:length(H2[[l]]), FUN = function(k)
  {
    bandwidth_evaluation(H2[[l]][k], p.seq[l], p.eval, n.seq, N, 
                         cov_ou, list(theta = theta, sigma = sigma), 
                         OU, list(alpha = theta, sigma = sigma, x0 = 0), 
                         eps.arg = list(sd = sd), m = m, grid.type = "without diagonal")
  },
  future.seed = T)),ncol = 4)
  
  toc()
}

parallel::stopCluster(cl)




#### 2RV process ####

n.seq = 100
p.seq =  50
##### m = 0 #####
m = 0

bw_comparison_OU_m0_2rv = list()
bw_comparison_OU_m0_full_2rv = list()

cl = parallel::makeCluster(parallel::detectCores( ) - 1)
future::plan(future::cluster)

for(l in 1:length(p.seq)){
  tic()
  bw_comparison_OU_m0_2rv[[l]]      = matrix(t(future_sapply(1:length(H[[l]]), FUN = function(k)
  {
    bandwidth_evaluation(H[[l]][k], p.seq[l], p.eval, n.seq, N, 
                         Gamma = cov_z_2rv, r.process = z_two_rv, 
                         eps.arg = list(sd = sd), m = m)
  },
  future.seed = T)),ncol = 4)
  cat("p =", p.seq[l], "done.")
  
  bw_comparison_OU_m0_full_2rv[[l]] = matrix(t(future_sapply(1:length(H[[l]]), FUN = function(k)
  {
    bandwidth_evaluation(H[[l]][k], p.seq[l], p.eval, n.seq, N,  
                         Gamma = cov_z_2rv, r.process = z_two_rv, 
                         eps.arg = list(sd = sd), m = m, grid.type = "without diagonal")
  },
  future.seed = T)),ncol = 4)
  
  toc()
}

parallel::stopCluster(cl)
##### m = 1 #####
m = 1

bw_comparison_OU_m1_2rv = list()
bw_comparison_OU_m1_full_2rv = list()

cl = parallel::makeCluster(parallel::detectCores( ) - 1)
future::plan(future::cluster)

for(l in 1:length(p.seq)){
  tic()
  bw_comparison_OU_m1_2rv[[l]]      = matrix(t(future_sapply(1:length(H[[l]]), FUN = function(k)
  {
    bandwidth_evaluation(H[[l]][k], p.seq[l], p.eval, n.seq, N, 
                         Gamma = cov_z_2rv, r.process = z_two_rv, 
                         eps.arg = list(sd = sd), m = m)
  },
  future.seed = T)),ncol = 4)
  cat("p =", p.seq[l], "done.")
  
  bw_comparison_OU_m1_full_2rv[[l]] = matrix(t(future_sapply(1:length(H[[l]]), FUN = function(k)
  {
    bandwidth_evaluation(H[[l]][k], p.seq[l], p.eval, n.seq, N,  
                         Gamma = cov_z_2rv, r.process = z_two_rv, 
                         eps.arg = list(sd = sd), m = m, grid.type = "without diagonal")
  },
  future.seed = T)),ncol = 4)
  
  toc()
}

parallel::stopCluster(cl)


##### m = 2 #####
m = 2

bw_comparison_OU_m2_2rv = list()
bw_comparison_OU_m2_full_2rv = list()

cl = parallel::makeCluster(parallel::detectCores( ) - 1)
future::plan(future::cluster)

for(l in 1:length(p.seq)){
  tic()
  bw_comparison_OU_m2_2rv[[l]]      = matrix(t(future_sapply(1:length(H[[l]]), FUN = function(k)
  {
    bandwidth_evaluation(H[[l]][k], p.seq[l], p.eval, n.seq, N, 
                         Gamma = cov_z_2rv, r.process = z_two_rv, 
                         eps.arg = list(sd = sd), m = m)
  },
  future.seed = T)),ncol = 4)
  cat("p =", p.seq[l], "done.")
  
  bw_comparison_OU_m2_full_2rv[[l]] = matrix(t(future_sapply(1:length(H[[l]]), FUN = function(k)
  {
    bandwidth_evaluation(H[[l]][k], p.seq[l], p.eval, n.seq, N,  
                         Gamma = cov_z_2rv, r.process = z_two_rv, 
                         eps.arg = list(sd = sd), m = m, grid.type = "without diagonal")
  },
  future.seed = T)),ncol = 4)
  
  toc()
}

parallel::stopCluster(cl)



#### Estimator Comparison ####
to_tibble = function(bw_comp, est = "mir", m = 1, Z = "OU") {
  Reduce(rbind, bw_comp) |>
    as_tibble() |>
    rename(n = V1, p = V2, h = V3, sup.err = V4) |>
    mutate(p = as.factor(p), 
           n = as.factor(n), 
           estimator = est,
           m = m,
           Z = Z)
}

# alrady contained in dataset
OU_m0 = to_tibble(bw_comparison_OU_m0, m = 0)
OU_m1 = to_tibble(bw_comparison_OU_m1)
OU_m2 = to_tibble(bw_comparison_OU_m2, m = 2)
OU_m0_full = to_tibble(bw_comparison_OU_m0_full, m = 0, est = "wd")
OU_m1_full = to_tibble(bw_comparison_OU_m1_full, est = "wd")
OU_m2_full = to_tibble(bw_comparison_OU_m2_full, m = 2, est = "wd")
rv2_m0 = to_tibble(bw_comparison_OU_m0_2rv, m = 0, Z  = "2rv")
rv2_m0_full = to_tibble(bw_comparison_OU_m0_full_2rv, m = 0, est = "wd", Z  = "2rv")
rv2_m1 = to_tibble(bw_comparison_OU_m1_2rv, Z  = "2rv")
rv2_m1_full = to_tibble(bw_comparison_OU_m1_full_2rv, est = "wd", Z  = "2rv")
rv2_m2 = to_tibble(bw_comparison_OU_m2_2rv, Z  = "2rv", m = 2)
rv2_m2_full = to_tibble(bw_comparison_OU_m2_full_2rv, est = "wd", Z  = "2rv", m = 2)

est_comp = rbind(OU_m0, OU_m1,
                 OU_m0_full, OU_m1_full,
                 OU_m2, OU_m2_full, 
                 rv2_m0, rv2_m1,
                 rv2_m0_full, rv2_m1_full,
                 rv2_m2, rv2_m2_full) |> 
  mutate(estimator = as.factor(estimator), 
         m = as.factor(m), Z = as.factor(Z))

##### Figure 8a ##### 
est_comp |> 
  filter(n == 100, p == 50, Z == "OU") |> 
  ggplot(aes(x = h, y = sup.err, pch = estimator, col = m)) + 
  geom_point() + 
  lims(y = c(0.03, 1.2)) +
  labs(title = "Ornstein-Uhlenbeck", subtitle = "n = 100, p = 50") + 
  my_theme
ggsave("cov/grafics/est_comp_OU_points.pdf", device = "pdf", width = 5, height = 3.8, units = "in")


##### Figure 8b #####
est_comp |> 
  filter(n == 100, p == 50, Z == "2rv") |> 
  ggplot(aes(x = h, y = sup.err, pch = estimator, col = m)) + 
  geom_point() + 
  lims(y = c(0.03, 1.1))+
  labs(title = "Process 2", subtitle = "n = 100, p = 50") +
  my_theme
ggsave("cov/grafics/est_comp_2rv_points.pdf", device = "pdf", width = 5, height = 3.8, units = "in")


