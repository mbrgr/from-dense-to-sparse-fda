library(ggplot2)
library(reshape2)
library(locpol)
library(interp)
library(stats)
library(future)
library(future.apply)
library(parallel)
library(tidyverse)

source("mean/functions_mean.R")
#my_theme = theme_grey(base_size = 15) + 
#  theme(plot.title = element_text(size = 14))



##### Comparison of different degrees for the local polynomial estimator #####
N = 1000
n = 600
p = c(50, 200, 500)
set.seed(1244) # 244

m = length(p)
H = sapply(p, function(x)(rev(seq(0.3, 3/x, -0.005))))
H

cl = makeCluster(detectCores( ) - 1);
plan(cluster);
erg_loc_lin = lapply(1:m, function(j){
  print(j)
  return(h_est_BB(n, p[j], N, H[[j]], deg = 1 )[-1]) 
})
erg_loc_quad = lapply(1:m, function(j){
  print(j)
  return(h_est_BB(n, p[j], N, H[[j]], deg = 2 )[-1]) 
})
erg_loc_cubic = lapply(1:m, function(j){
  print(j)
  return(h_est_BB(n, p[j], N, H[[j]][-1], deg = 3 )[-1]) 
})
erg_loc_4 = lapply(1:m, function(j){
  print(j)
  return(h_est_BB(n, p[j], N, H[[j]][-(1:5)], deg = 4 )[-1]) 
})
stopCluster(cl);


comp_tibble = tibble(sup.err  = c(unlist(erg_loc_lin),
                                  unlist(erg_loc_quad), 
                                  unlist(erg_loc_cubic),
                                  unlist(erg_loc_4)),
                     h = c(rep(unlist(H), 2), H[[1]][-1], H[[2]][-1], H[[3]][-1], H[[1]][-(1:5)], H[[2]][-(1:5)], H[[3]][-(1:5)]),
                     p = factor(c(rep(p, erg_loc_lin |> sapply(length)), 
                                  rep(p, erg_loc_quad |> sapply(length)),
                                  rep(p, erg_loc_cubic |> sapply(length)),
                                  rep(p, erg_loc_4 |> sapply(length)))),
                     est = as.factor(c(rep("lin", length(unlist(erg_loc_lin))),
                                       rep("quad", length(unlist(erg_loc_quad))),
                                       rep("cubic", length(unlist(erg_loc_cubic))),
                                       rep("quartic", length(unlist(erg_loc_4)))  ) ))

comp_tibble |> 
  ggplot(aes(x = h, y= sup.err, col = p, pch = est)) + 
#  geom_line() + 
  geom_point() + 
  lims( y = c(0.01, 0.31)) + 
  scale_linetype_manual(values = 2:5, 
                        breaks = c("lin", "quad", "cubic", "quartic")) + 
  labs(subtitle = "n = 600") + 
  deriv_est_theme


ggsave("degree_comparison.pdf", device = "pdf", width = 5, height = 3.8, units = "in")

save.image("mean/data/results_degree_comparison.RData")
