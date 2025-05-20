library(biLocPol) # please install this package from Github first. See "README.md" file for instructions
library(future.apply)
library(tidyverse)
library(tictoc)
library(plotly)

my_theme = theme_grey(base_size = 15) + 
  theme(plot.title = element_text(size = 14))


#### results ####
load("cov/data/results_3_11.RData")

#### functions ####
source("cov/functions.R")


###### Error Decomposition #####
N = 500 
n.seq = c(50, 100, 200, 400)
p.seq = c(15, 25, 50, 75)

p.eval = 50 #p.eval = 100
H = lapply(1:length(p.seq), function(l){seq(1, max(5/p.seq[l], 0.1), -0.05)})

# Parameter OU Process
theta = 2; sigma = 3
# Standard deviation for additional errors
sd = 0.25

cl = parallel::makeCluster(parallel::detectCores( ) - 1)
future::plan(future::multisession)

error_decomp_arr = array(0, c(length(p.seq), 5, length(n.seq)))

set.seed(87)
for(i in 1:length(n.seq)){
  h_min = bw_comparison_tibble |> 
    group_by(n, p) |> 
    filter(n == n.seq[i]) |> 
    slice_min(sup.err) |> 
    pull(h)
  for (j in 1:length(p.seq)) {
    tic()
    weight = local_polynomial_weights(p.seq[j], h_min[j], p.eval, parallel = T, m = 2,  del = 1, 
                                      parallel.environment = F)
    error_decomp_arr[j,,i] = error_decomposition_deriv(weight, n.seq[i], N, parallel = T, 
                                                       parallel.environment = F, correction  = T)
    cat("n =", n.seq[i], "and p =", p.seq[j], "done. ")
    toc()
  }
  
}

parallel::stopCluster(cl)

rm(weight)

dimnames(error_decomp_arr) = list(p.seq, c("eps", "dsc", "prc", "mix", "sup"), n.seq)
err_dec_tib = error_decomp_arr %>%  
  apply(2, rbind) %>% 
  as_tibble() %>% 
  mutate(n = rep(n.seq, each = 4), 
         p = rep(p.seq, 4),
         h = bw_comparison_tibble |> 
           group_by(n, p) |> 
           slice_min(sup.err) |> 
           pull(h)) %>% 
  pivot_longer(1:5, names_to = "term", values_to = "sup.err")

err_dec_tib %>% 
  ggplot() + 
  geom_line(aes(x = p, y = sup.err, col = term, lty = term)) + 
  facet_wrap(n~., nrow = 1)

##### Figure 3.11 ####
err_dec_tib %>% 
  ggplot(aes(x = n, y = sup.err, col = term, lty = term, pch = term)) + 
  geom_line() + 
  geom_point() +
  facet_wrap(p ~., nrow = 1) + 
  my_theme
ggsave("cov/grafics/error_decomp_cov_deriv_OU.pdf", device = "pdf", unit = "in", width = 10, height = 4.5)

save.image("cov/data/results_3_11.RData")
