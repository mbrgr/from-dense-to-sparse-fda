#### packages, code ####
library(ggplot2)
library(tidyverse)
library(biLocPol)
library(future.apply)
library(tictoc)

my_theme = theme_grey(base_size = 15) + 
  theme(plot.title = element_text(size = 14))

#  Figures 3.4, 3.6
source("cov/functions.r")

#### calculations ####

##### setup #####
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

##### One-Fold Cross Validation #####

one_fold_cv = list()

cl = parallel::makeCluster(parallel::detectCores( ) - 1)
future::plan(future::cluster)

for(l in 1:length(p.seq)){
  tic()
  one_fold_cv[[l]] =  lopocv_sim(N, n = 100, p.seq[l], H[[l]], 
                                 m = 1, w.parallel = T, 
                                 theta = theta, sigma = sigma, sd = sd)
  cat("p =", p.seq[l], "done.")
  toc()
}

parallel::stopCluster(cl)

##### 5-Fold Cross Validation #####

five_fold_cv = list()

cl = parallel::makeCluster(parallel::detectCores( ) - 1)
future::plan(future::cluster)

for(l in 1:length(p.seq)){
  tic()
  five_fold_cv[[l]] = k_fold_cv_simulation(N, n = 100, p.seq[l], H[[l]], 
                                           K = 5, m = 1, w.parallel = T, 
                                           theta = theta, sigma = sigma, sd = sd)
  cat("p =", p.seq[l], "done.")
  toc()
}

parallel::stopCluster(cl)



##### One-Fold Cross Validation #####
# different n

one_fold_cv_n400 = list()

cl = parallel::makeCluster(parallel::detectCores( ) - 1)
future::plan(future::cluster)

for(l in 1:3){
  tic()
  one_fold_cv_n400[[l]] =  lopocv_sim(N, n = 400, p.seq[l], H[[l]], 
                                      m = 1, w.parallel = T, 
                                      theta = theta, sigma = sigma, sd = sd)
  cat("p =", p.seq[l], "done.")
  toc()
}

parallel::stopCluster(cl)

##### 5-Fold Cross Validation #####

five_fold_cv = list()

cl = parallel::makeCluster(parallel::detectCores( ) - 1)
future::plan(future::cluster)

for(l in 1:length(p.seq)){
  tic()
  five_fold_cv[[l]] = k_fold_cv_simulation(N, n = 100, p.seq[l], H[[l]], 
                                           K = 5, m = 1, w.parallel = T, 
                                           theta = theta, sigma = sigma, sd = sd)
  cat("p =", p.seq[l], "done.")
  toc()
}

parallel::stopCluster(cl)

##### 5-Fold Cross Validation #####

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

##### 2-Fold Cross Validation #####

two_fold_cv = list()

cl = parallel::makeCluster(parallel::detectCores( ) - 1)
future::plan(future::cluster)

for(l in 1:length(p.seq)){
  tic()
  two_fold_cv[[l]] = k_fold_cv_simulation(N, n, p.seq[l], H[[l]], 
                                          K = 5, m = 1, w.parallel = T, 
                                          theta = theta, sigma = sigma, sd = sd)
  cat("p =", p.seq[l], "done.")
  toc()
}



#### evaluation ####
##### One -Fold CV #####
###### results ######
load("cov/data/one_fold_cv_n400_OU.RData")
one_fold_tbl_n400 = tibble(h = one_fold_cv_n400 %>% unlist(), p = gl(3, 1000, labels = p.seq))
one_fold_tbl_n400 %>% 
  ggplot(aes(y = h, x = p, col = p)) + 
  geom_boxplot()

n400_table = one_fold_tbl_n400 %>% 
  group_by_all() %>% 
  summarise(n = n()/1000)

n400_table %>% 
  ggplot(aes(h, n)) + 
  geom_point(size = 3) + 
  ylim(c(0, 0.1))+
  facet_wrap(.~p)

one_fold_tbl_n400 %>% 
  summarise(.by = p, mean(h))



##### five fold cv #####
#load("cov/data/five_fold_cv_OU.RData")
# five_fold_tbl = tibble(h = five_fold_cv %>% unlist(), p = gl(5, 1000, labels = p.seq))
# five_fold_tbl %>% 
#   ggplot(aes(y = h, x = p, col = p)) + 
#   geom_boxplot()
# 
# 
# five_fold_table = five_fold_tbl %>% 
#   group_by_all() %>% 
#   summarise(n = n()/1000)
# 
# # Figure in Appendix
# five_fold_table %>% 
#   ggplot(aes(h, n)) + 
#   geom_point(size = 3) + 
#   ylim(c(0, 0.1))+
#   facet_wrap(.~p) +
#   ylab(NULL) + 
#   labs(title = "n = 100")
# 
# five_fold_tbl %>% 
#   summarise(.by = p, mean(h))


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

ggsave("cov/grafics/cov_5fcv_bw_n400_table.pdf", device = "pdf", width = 8, height = 5, units = "in")

five_fold_tbl_n400 %>% 
  summarise(.by = p, mean(h))
