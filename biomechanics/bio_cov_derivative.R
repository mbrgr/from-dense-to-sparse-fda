library(ffscb)
library(tidyverse)
library(reshape2)
library(locpol)
library(biLocPol)
library(plotly)
library(MASS)

my_theme = theme_grey(base_size = 15) + 
  theme(plot.title = element_text(size = 14))

#load("R-Codes/data/biomechanics_covariance_estimation_h2.RData")

#### load data ####
data("Biomechanics")
Biomechanics
head(Biomechanics)
?Biomechanics

#### preprocess data as for mean estimation  ####
stance_phase = Biomechanics[,1]
extra.mat = t(Biomechanics[, 2:19])
normal.mat = t(Biomechanics[,20:37])
diff.mat = extra.mat - normal.mat

##### estimation of the variance ####
p = 201
p.eval = 201
var_weights = biLocPol::local_polynomial_weights(p, 0.2, p.eval, F, m = 2, del = 1, eval.type = "diagonal")

Z_normal = observation_transformation(normal.mat)
Z_extra  = observation_transformation(extra.mat)
Z_diff  = observation_transformation(diff.mat)

normal_var = eval_weights(var_weights, Z_normal)
extra_var  = eval_weights(var_weights, Z_extra)
diff_var = eval_weights(var_weights, Z_diff)


var_tib = tibble( x   = rep(stance_phase, 3), 
                  var = c(normal_var[,1], 
                          extra_var[,1], 
                          diff_var[,1]), 
                  cushion = gl(3, p.eval, labels = c("normal", "extra", "diff")))

##### Figure 3.17(a) #####
var_tib |> 
  filter(cushion != "diff") |> 
  ggplot(aes(x = x)) + 
  geom_line(aes(y = var, color = cushion, linetype = cushion)) + 
  labs(y = "Nm/Kg", 
       x = "% of stance phase") + 
  scale_linetype_manual(values = c(2,3)) + 
  scale_colour_manual(values = c("red", "blue")) + 
  my_theme
ggsave("biomechanics/grafics/biomechanics_var_normal_extra.pdf", device = "pdf", unit = "in", width = 5, height = 3.8)

##### Figure 3.17(b) #####
var_tib |> 
  filter(cushion == "diff") |> 
  ggplot(aes(x = x)) + 
  geom_line(aes(y = var, linetype = cushion)) + 
  labs(y = "Nm/Kg", 
       x = "% of stance phase") + 
  my_theme
ggsave("biomechanics/grafics/biomechanics_var_diff.pdf", device = "pdf", unit = "in", width = 5, height = 3.8)

del_var_tib  = tibble(x    = rep(stance_phase, 6), 
                      est = c(normal_var[,2], normal_var[,3], 
                              extra_var[,2], extra_var[,3], 
                              diff_var[,2], diff_var[,3]), 
                      del = rep(gl(2, p.eval, labels = c("G10", "G01")), 3), 
                      cushion = gl(3, 2*p.eval, labels = c("normal", "extra", "diff")))

##### Figure 3.18(a) #####
del_var_tib |>
  filter(cushion != "diff") |> 
  ggplot(aes(x = x)) + 
  geom_line(aes(y = est, color = cushion, lty = del))+ 
  scale_linetype_manual(name = "Deriv.", values = c(2,3), labels = c("G10" = expression(italic(d)^{"(1,0)"}*Gamma), "G01" = expression(italic(d)^{"(0,1)"}*Gamma))) + 
  scale_colour_manual(values = c("red", "blue")) + 
  labs(y = "Nm/Kg", 
       x = "% of stance phase") + 
  my_theme
ggsave("biomechanics/grafics/biomechanics_del_var_normal_extra.pdf", device = "pdf", unit = "in", width = 5, height = 3.8)

##### Figure 3.18(b) #####
del_var_tib |>
  filter(cushion == "diff") |> 
  ggplot(aes(x = x)) + 
  geom_line(aes(y = est, lty = del)) +
  scale_linetype_manual(name = "Deriv.", values = c(2,3), labels = c("G10" = expression(italic(d)^{"(1,0)"}*Gamma), "G01" = expression(italic(d)^{"(0,1)"}*Gamma))) + 
  labs(y = "Nm/Kg", 
       x = "% of stance phase") + 
  my_theme

ggsave("biomechanics/grafics/biomechanics_del_var_diff.pdf", device = "pdf", unit = "in", width = 5, height = 3.8)


##### estimation of the covariance kernel of the difference #####
p.eval = 25
cov_weights = biLocPol::local_polynomial_weights(p, 0.2, p.eval, F, m = 2, del = 1)

Z = observation_transformation(diff.mat)
cov_est = eval_weights(cov_weights, Z)[,,1]
cov_01_est = eval_weights(cov_weights, Z)[,,2]

x.eval = (1:p.eval - 0.5)/p.eval

plot_ly() |> 
  add_surface(x = ~x.eval, y = ~x.eval, z = cov_est)


UP = upper.tri(cov_01_est, diag = T)
cov_01_est_up = cov_01_est
cov_01_est_dn = cov_01_est
cov_01_est_dn[UP] = NA
cov_01_est_up[!UP] = NA
plot_ly() %>% 
  add_surface(x = ~ time, y = ~ time, z = cov_01_est_up) %>% 
  add_surface(x = ~ time, y = ~ time, z = cov_01_est_dn)




