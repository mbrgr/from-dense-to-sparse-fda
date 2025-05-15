#### packages, code ####
library(ggplot2)
library(tidyverse)
library(biLocPol)

my_theme = theme_grey(base_size = 15) + 
  theme(plot.title = element_text(size = 14))

#### results ####
load("cov/data/estimator_comparison.RData")

#  Figures 3.6
source("cov/functions.r")


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


