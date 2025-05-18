#### packages, code ####
library(ggplot2)
library(tidyverse)
library(biLocPol)
library(tictoc)
library(plotly)

my_theme = theme_grey(base_size = 15) + 
  theme(plot.title = element_text(size = 14))

#### results ####
load("cov/data/bandwidth_comparison_derivative_cov_OU.RData")

# 
source("cov/functions.r")


##### Bandwidth Comparison for Covariance derivative estimation #####
N = 500 
n.seq = c(50, 100, 200, 400)
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

bw_comparison_tibble = Reduce(rbind, bw_comparison)
colnames(bw_comparison_tibble) = c("n", "p", "h", "sup.err")
bw_comparison_tibble = bw_comparison_tibble %>% 
  as_tibble() %>% 
  mutate(p = as.factor(p))

##### Figure 3.10 #####
ggplot(bw_comparison_tibble) + 
  geom_point(aes(x = h, y = sup.err, col = p)) + 
  lims(y = c(0, 10.4)) +
  facet_wrap(n ~., nrow = 1) + 
  my_theme

ggsave("cov/grafics/bw_comp_cov_deriv_OU.pdf", device = "pdf", width = 10, height = 5.5, units = "in")

bw_comparison_tibble |> 
  group_by(n, p) |> 
  slice_min(sup.err)


save.image("cov/data/bandwidth_comparison_derivative_cov_OU.RData")




