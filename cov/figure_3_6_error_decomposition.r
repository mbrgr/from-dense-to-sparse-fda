#### packages, code ####
library(ggplot2)
library(tidyverse)
library(biLocPol)

my_theme = theme_grey(base_size = 15) + 
  theme(plot.title = element_text(size = 14))

#### results ####
load("cov/data/error_decomp.RData")

#  Figures 3.6
source("cov/functions.r")

##### Error Decomposition #####
error_decomp_arr |> dimnames() = list(p.seq, c("eps", "dsc", "prc", "mix", "sup"), n.seq)
error_decomp_tbl = error_decomp_arr |> as_tibble() |> 
  pivot_longer(cols = everything(),
               names_to = c("term", "n"), 
               names_pattern = "(...).(.*)",
               names_transform = list(term = as.factor, n = as.double),
               values_to = "error") |> 
  mutate(p = rep(p.seq, each = 20))

error_decomp_tbl |> 
  filter(n == 400, term == "sup") |> 
  print(n = 100) 

##### Figure 3.6 #####
error_decomp_tbl |> 
  filter(p != 75) |> 
  ggplot(aes(x = n, y = error, col = term, lty = term)) + 
  geom_line(size = .6) + 
  facet_wrap(.~p, nrow = 1) +
  my_theme + 
  labs(y = "sup.error")

ggsave("cov/grafics/error_decomp_different_p.png", device = "png", width = 8, height = 5, units = "in")

# Not in Paper
error_decomp_tbl |> 
  ggplot(aes(x = p, y = error, col = term, lty = term)) + 
  geom_line() + 
  facet_wrap(.~n)