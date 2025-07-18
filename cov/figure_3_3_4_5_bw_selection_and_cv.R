#### packages, code ####
library(ggplot2)
library(tidyverse)
library(biLocPol)

my_theme = theme_grey(base_size = 15) + 
  theme(plot.title = element_text(size = 14))

#### results ####
load("cov/data/bw_comp_OU.RData")

#  Figures 3.4, 3.5 and 3.6
# further Figure B.2
source("cov/functions.r")

#### Bandwidth comparison ####
##### Calculation #####
N = 1000
n.seq = c(50, 100, 200, 400)
p.seq = c(15, 25, 50, 75, 100)
# p.seq = 25

p.eval = 100
H = lapply(1:length(p.seq), function(l){seq(1, 3/p.seq[l], -0.05)})
# H = lapply(1:length(p.seq), function(l){seq(.6, 3/p.seq[l], -0.2)})

# Parameter OU Process
theta = 2; sigma = 3
# Standard deviation for additional errors
sd = 0.75

bw_comparison = list()

cl = parallel::makeCluster(parallel::detectCores( ) - 1)
future::plan(future::cluster)

for(l in 1:length(p.seq)){
  tic()
  bw_comparison[[l]] = matrix(t(future_sapply(1:length(H[[l]]), FUN = function(k)
  {
    bandwidth_evaluation(H[[l]][k], p.seq[l], p.eval, n.seq, N, 
                         cov_ou, list(theta = theta, sigma = sigma), 
                         OU, list(alpha = theta, sigma = sigma, x0 = 0), 
                         eps.arg = list(sd = sd))
  },
  future.seed = T)),ncol = 4)
  cat("p =", p.seq[l], "done.")
  toc()
}

parallel::stopCluster(cl)


##### Figure 3.3 #####
bw_comparison_tbl = Reduce(rbind, bw_comparison) |>
  as_tibble() |>
  rename(n = V1, p = V2, h = V3, sup.err = V4) |>
  mutate(p = as.factor(p), n = as.factor(n)) 

bw_comparison_tbl |> 
  filter(n == 400) |> 
  ggplot() + 
  geom_point(aes(x = h, y = sup.err, col = p, pch = p)) + 
  ylim(c(0.02, 0.52)) +
  labs(title = "n = 400") + 
  my_theme

ggsave("cov/grafics/cov_optimal_bw_n400.pdf", device = "pdf", width = 5, height = 3.8, units = "in")

##### Figure B.2  ######
bw_comparison_tbl |> 
  ggplot() + 
  geom_point(aes(x = h, y = sup.err, col = p, pch = p)) + 
  lims(y = c(0, .95)) +
  facet_wrap(n~., nrow = 1)  + 
  my_theme
ggsave("cov/grafics/cov_optimal_bw_various_n.pdf", device = "pdf", width = 9, height = 4, units = "in")

min_h_tibble = bw_comparison_tbl |>
  group_by(n, p) |> 
  slice_min(sup.err) 

#### Cross Validation ####
##### Calculation #####
five_fold_cv_n400 = list()

cl = parallel::makeCluster(parallel::detectCores( ) - 1)
future::plan(future::cluster)

for(l in 1:length(p.seq)){
  tic()
  five_fold_cv_n400[[l]] = k_fold_cv_simulation(N, n = 400, p.seq[l], H[[l]], 
                                                K = 5, m = 1, w.parallel = T, 
                                                theta = theta, sigma = sigma, sd = sd)
  cat("p =", p.seq[l], "done. ")
  toc()
}

parallel::stopCluster(cl)


##### Figure 3.4 #####
five_fold_tbl_n400 = tibble(h = five_fold_cv_n400 %>% unlist(), p = gl(5, 1000, labels = p.seq))

five_fold_tbl_n400 %>% 
  ggplot(aes(y = h, x = p, col = p)) + 
  geom_boxplot(size = .6) + 
  labs(title = "n = 400") +
  my_theme + 
  lims(y = c(0,1))

ggsave("cov/grafics/cov_5fcv_bw_n400.pdf", device = "pdf", width = 5, height = 3.8, units = "in")

five_fold_table_n400 = five_fold_tbl_n400 %>% 
  group_by_all() %>% 
  summarise(n = n()/1000)

##### Figure 3.5 #####
five_fold_table_n400 %>% 
  ggplot(aes(h, n)) + 
  geom_point(size = 3) + 
  facet_wrap(.~p, nrow = 1) +
  ylab(NULL) + 
  labs(title = "n = 400") + 
  my_theme

ggsave("cov/grafics/cov_5fcv_bw_n400_table.pdf", device = "pdf", width = 10, height = 5.5, units = "in")

five_fold_tbl_n400 %>% 
  summarise(.by = p, mean(h))

save_image("cov/data/bw_comp_OU.RData")
