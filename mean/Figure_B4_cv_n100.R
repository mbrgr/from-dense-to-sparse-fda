#### results ####
load("cov/data/results_cv_n100.RData")


#### Packages ####

library(biLocPol)
library(future.apply)
library(tidyverse)
library(tictoc)

my_theme = theme_grey(base_size = 15) +
 theme(plot.title = element_text(size = 14))


#### Bandwidth Comparison ####

N = 500
n.seq = c(50, 100, 200, 400)
p.seq = c(15, 25, 50, 75, 100)
# p.seq = 25

p.eval = 100
H = lapply(1:length(p.seq), function(l){seq(1, 3/p.seq[l], -0.05)})
# H = lapply(1:length(p.seq), function(l){seq(.6, 3/p.seq[l], -0.2)})

# Parameter OU Process
theta = 2; sigma = 3
# Standard deviation for additional errors
sd = 0.5

##### function #####

k_fold_cv_simulation = function(N, n, p, h.seq, K = 5, m = 1, w.parallel = T, theta = 2, sigma = 3, sd = 0.75){
  
  grp = sample(rep(1:K, ceiling(n/K)), n)
  
  help = function(h){
    w_h = local_polynomial_weights(p, h, p.eval = p, m = m, parallel = w.parallel, parallel.environment = F)
    max_diff = numeric(K)
    
    future_replicate(N, {
      Y = FDA_observation(n, x.design = (1:p - 0.5)/p, f = biLocPol::mu,
                          r.process = OU, process.arg = list(alpha = theta, sigma = sigma, x0 = 0), eps.arg = list(sd = sd))
      for(kk in 1:K){
        test_grp     = matrix(observation_transformation(Y[grp == kk,], grid.type = "full"), p, p)
        train_grp    = observation_transformation(Y[grp != kk,], grid.type = "less")
        pred         = eval_weights(w_h, train_grp)
        max_diff[kk] = max(abs((test_grp - pred)[!as.logical(diag(p))]))
      }
      mean(max_diff)
    }, future.seed = T)
  }
  
  
  mean_sup = sapply(h.seq, help)
  
  h.seq[apply(mean_sup, 1, which.min)]
}

##### calculation #####

set.seed(43)
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


##### evaluation #####
five_fold_tbl = tibble(h = five_fold_cv %>% unlist(), p = gl(5, 1000, labels = p.seq))
five_fold_tbl %>% 
  ggplot(aes(y = h, x = p, col = p)) + 
  geom_boxplot()+
  theme(text = element_text(size = 18)) 

five_fold_table = five_fold_tbl %>% 
  group_by_all() %>% 
  summarise(n = n()/1000)

five_fold_table %>% 
  ggplot(aes(h, n)) + 
  geom_point(size = 3) + 
  ylim(c(0, 0.1))+
  facet_wrap(.~p, 1) +
  ylab(NULL) + 
  labs(title = "n = 100")  + 
  my_theme

ggsave("cov/grafics/5fold_cv_n100.pdf", 
       device = "pdf", 
       unit = "in", 
       width = 10, 
       height = 5.5)

five_fold_tbl %>% 
  summarise(.by = p, mean(h)) 

save.image("cov/data/results_cv_n100.RData")
