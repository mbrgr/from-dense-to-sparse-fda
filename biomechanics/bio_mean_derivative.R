# Real Data application: Biomechanics data set

#### source and packages ####
source("mean/functions_mean.R")
library(ffscb)
library(reshape2)
library(locpol)
library(interp)
library(stats)
library(future)
library(future.apply)
library(parallel)
library(tidyverse)
library(plotly)
library(biLocPol) # has to be installed from Github


#### results ####

#### load data ####
data("Biomechanics")
Biomechanics
head(Biomechanics)
?Biomechanics

#### preprocessing ####
# generating matrices with the uncushioned data, the cushioned data and the difference
extra.mat = cbind(Biomechanics[, 1:19], rowMeans(Biomechanics[, 2:19]))
colnames(extra.mat) = c("stance_phase", 1:18, "mean")
normal.mat = cbind(Biomechanics[, c(1, 20:37)], rowMeans(Biomechanics[, 20:37]))
colnames(normal.mat) = c("stance_phase", 1:18, "mean")
diff.mat = cbind(Biomechanics[, 1], 
                 Biomechanics[, 2:19] - Biomechanics[, 20:37], 
                 rowMeans(Biomechanics[, 2:19] - Biomechanics[, 20:37]))
colnames(diff.mat) = c("stance_phase", 1:18, "mean")

#### melt to form of tidy data ####
extra = melt(extra.mat, id.vars = "stance_phase", 
             variable.name = "person", 
             value.name = "pressure")
normal = melt(normal.mat, id.vars = "stance_phase", 
              variable.name = "person", 
              value.name = "pressure")
diff = melt(diff.mat, id.vars = "stance_phase", 
            variable.name = "person", 
            value.name = "pressure")

acast(extra, stance_phase ~ person, value.var = "pressure")

normal$shoe = "normal"
extra$shoe = "extra"
diff$shoe = "extra - normal"

df.both = rbind(normal, extra)

normal |> 
  filter(shoe == "normal", person == "mean") |> 
  mutate(stance_phase = stance_phase/100) |> 
  mutate(difference_quot = c(0,diff(pressure)/diff(stance_phase))) |> 
  pull(difference_quot)

#### Plots of the data ####
##### Figure 5 #####
ggplot() + geom_point(aes(x = stance_phase, y = pressure, group = person), 
                      df.both[df.both$person != "mean", ], size = .2) + facet_wrap(shoe ~.)+ 
  labs(x = "% of stance phase", y = "Nm/Kg", subtitle = "Pressure curves of normal and extra cushioned shoe") 
ggsave("biomechanics/grafics/pressure_curves.pdf", device = "pdf",
       width = 5, height = 3.8, units = "in")

normal |> 
  filter(shoe == "normal", person == "mean") |> 
  mutate(stance_phase = stance_phase/100) |> 
  mutate(difference_quot = c(0,diff(pressure)/diff(stance_phase))) |> 
  pull(difference_quot)

##### Figure 7 #####
ggplot()  +
  geom_point(aes(x = stance_phase, y = pressure, group = person), 
             diff[diff$person != "mean", ], size = .3, alpha= .3, color = "black") +
  geom_point(aes(x = stance_phase, y = pressure), diff[diff$person == "mean", ], 
             size = .9, col = "black") +
  labs(x = "% of stance phase", y = "Nm/Kg", subtitle = "Difference of pressure") + 
  deriv_est_theme
ggsave("biomechanics/grafics/diff_pressure_curves.pdf", device = "pdf",
       width = 5, height = 3.8, units = "in")


#### estimation of the derivative ####
diff
full_data = diff |> tibble() |>  mutate(shoe = shoe |> as_factor()) |> 
  rbind(df.both)
full_data |>  summary()
mean_extra = full_data |> filter(person == "mean", shoe == "extra")
mean_normal = full_data |> filter(person == "mean", shoe == "normal")
mean_diff = full_data |> filter(person == "mean", shoe == "extra - normal")
bw = 6.5
est_extra_derivative = locPolSmootherC(mean_extra$stance_phase, 
                                       mean_extra$pressure,
                                       mean_extra$stance_phase,
                                       bw = bw, deg = 3, kernel = EpaK) # bw = 5.5
est_normal_derivative = locPolSmootherC(mean_normal$stance_phase, 
                                        mean_normal$pressure,
                                        mean_normal$stance_phase,
                                        bw = bw, deg = 3, kernel = EpaK) # bw = 5.5
est_diff_derivative = locPolSmootherC(mean_diff$stance_phase, 
                                      mean_diff$pressure,
                                      mean_diff$stance_phase,
                                      bw = bw, deg = 3, kernel = EpaK) # bw = 5.5

# Illustrations: not contained in Paper
est_extra_derivative$beta1 |>  plot(x = mean_normal$stance_phase, type = "l")
lines(est_normal_derivative$beta1, col = "red", x= mean_normal$stance_phase)

est_diff_derivative$beta1 |> plot(x = mean_diff$stance_phase, type = "l")
abline(h = 0, lty = 3)

deriv_tibble = tibble(estimation = c(est_extra_derivative$beta1, 
                                     est_normal_derivative$beta1, 
                                     est_diff_derivative$beta1 ), 
                      x = rep(est_extra_derivative$x, 3), 
                      cushion = c(rep("extra", length(est_extra_derivative$beta1 )), 
                                  rep("normal",length(est_extra_derivative$beta1 ) ), 
                                  rep("diff", length(est_extra_derivative$beta1 ))))

##### Figure 6 #####
deriv_tibble |>
  filter(cushion != "diff") |> 
  ggplot(aes(x = x, y = estimation)) + 
  geom_line(aes(col = cushion, lty = cushion), size = 1) +
  labs(x = "% of stance phase", y = NULL, 
       title = "Derivative estimation for pressure curves")  + 
  scale_linetype_manual(values = c(1,4))+ 
  deriv_est_theme
ggsave("biomechanics/grafics/deriv_est.pdf", device = "pdf",
       width = 5, height = 3.8, units = "in")


##### Figure 9 #####
rbind(est_extra_derivative, est_normal_derivative) |> 
  mutate(shoe = rep(c("extra", "normal"), each = 201)) |> 
  rename("phase" = x) |> 
  ggplot(mapping = aes(x = beta1, y = beta0, color = phase, shape = shoe)) + 
  geom_point(size = .7) + 
  geom_point(size = .7) + 
  geom_hline(yintercept = 0, lty = 2) + 
  geom_vline(xintercept = 0, lty = 2) + 
  labs(y = "Nm/kg", x = "d/dx Nm/kg", title = "Phase-plane plot") + 
  scale_shape_manual(values = c(3,6))+ 
  deriv_est_theme
ggsave("biomechanics/grafics/phase_plane.pdf", device = "pdf",
       width = 5, height = 4.3, units = "in")

##### Figure 8 #####
deriv_tibble |>
  filter(cushion == "diff") |> 
  ggplot(aes(x = x, y = estimation)) + 
  geom_line() +
  labs(x = "% of stance phase", y = NULL, 
       title = "Derivative estimation for difference") + 
  deriv_est_theme
ggsave("biomechanics/grafics/deriv_est_diff.pdf", device = "pdf",
       width = 5, height = 3.8, units = "in")

#### estimation of derivative of the covariance kernel ####
p = seq(0, 1, 0.005) 
weights_cov_estimation  = biLocPol::local_polynomial_weights(p = length(p),
                                                             h = .2, 
                                                             p.eval = length(p)/2, 
                                                             parallel = T, 
                                                             m = 2, 
                                                             del = 2, 
                                                             eval.type = "diagonal")

# evaluation of the weights for the difference of the processes 
est = eval_weights(weights_cov_estimation, observation_transformation(t(as.matrix(diff.mat[,2:19]))))

diag_est_tib = tibble(x = rep(seq(0, 1, length.out = 100), 2) * 100, 
                      Gamma = c(est[,2],est[,3]), 
                      est = gl(2, 100, labels = c("Gamma10", "Gamma01")))

##### Figure 10 #####
diag_est_tib |> 
  ggplot(aes( x = x, y = Gamma, color = est, linetype = est)) + 
  geom_line(size = .9) +
  deriv_est_theme + 
  labs(y = NULL, x = "% of stance phase", title = expression(Partial~derivatives~of~Gamma~on~the~diagonal) ) + 
  scale_discrete_manual(
    aesthetics = c("color", "linetype"),
    values = c(2,4), 
    name = "Deriv.",
    labels = c(expression(italic(d)^{"(1,0)"}*Gamma), expression(italic(d)^{"(0,1)"}*Gamma))
  ) + 
  deriv_est_theme

ggsave("biomechanics/grafics/part_deriv_gamma_diagonal.pdf", device = "pdf", width = 5, height = 3.8, units = "in")
save.image("biomechanics/data/bio_derivative.RData")


##### Evaluation with the full estimation on [0,1]^2: not contained in the paper #####
# calculation of weights is cost expensive
load("mean/data/biomechanics_covariance_estimation_h35.RData")
est[,,1] |> diag() > 0
gamma_11 = est[,,5]
z = est[,,3]
x = weights_cov_estimation$x.eval
x = seq(0,1, 0.005)
plot_ly() |> 
  add_surface(x = x, y = x, z = ~est[,,1])
plot_ly() |> 
  add_surface(x = x, y = x, z = ~est[,,3]) # Processes seem to be not differentiable
plot_ly() |> 
  add_surface(x = x, y = x, z = ~est[,,5])

save.image("mean/data/biomechanics_derivative_estimation_weights_cov.RData")


#### Data Reduction ####
bw = 5.5
est.normal = locPolSmootherC(diff$stance_phase[diff$person == "mean"], 
                             diff$pressure[diff$person == "mean"], 
                             diff$stance_phase[diff$person == "mean"], 
                             bw = bw, deg = 2, kernel = EpaK)


ggplot()  +
  geom_point(aes(x = stance_phase, y = pressure, group = person), 
             diff[diff$person != "mean", ], size = .3, alpha= .7) +
  geom_point(aes(x = stance_phase, y = pressure), diff[diff$person == "mean", ], 
             size = .9, col = "red") +
  geom_line(aes(x = stance_phase, y = pressure), diff[diff$person == "mean", ], 
            linewidth = .4, alpha = .8,col = "red") +
  geom_line(aes(x = x, y = beta0), est.normal, col = "blue") + 
  labs(x = "% of stance phase", y = "Nm/Kg", subtitle = "Difference curves")


subs = seq(1, 201, 4)
length(subs)
bw.red = cv.locpol((diff$stance_phase[diff$person == "mean"])[subs]/100,
                   diff.mat[subs,2:19], h = seq(0.05, 0.2, 0.005))*100
est.normal.red = locPolSmootherC((diff$stance_phase[diff$person == "mean"])[subs], 
                                 diff$pressure[diff$person == "mean"][subs], 
                                 xeval = seq(0, 100, 0.5), 
                                 bw = bw.red+2, deg = 2, kernel = EpaK)

subs2 = seq(1, 201, 5)
length(subs2)
bw.red2 = cv.locpol((diff$stance_phase[diff$person == "mean"])[subs2]/100,
                    diff.mat[subs2,2:19], h = seq(0.1, 0.2, 0.005))*100
est.normal.red2 = locPolSmootherC((diff$stance_phase[diff$person == "mean"])[subs2], 
                                  diff$pressure[diff$person == "mean"][subs2], 
                                  #xeval = seq(1, 100, 1), 
                                  xeval = seq(0, 100, 0.5), 
                                  bw = bw.red2+2, deg = 2, kernel = EpaK)

subs3 = seq(1, 201, 10)
length(subs3)
bw.red3 = cv.locpol((diff$stance_phase[diff$person == "mean"])[subs3]/100,
                    diff.mat[subs3,2:19], h = seq(0.1, 0.2, 0.005))*100
est.normal.red3 = locPolSmootherC((diff$stance_phase[diff$person == "mean"])[subs3], 
                                  diff$pressure[diff$person == "mean"][subs3], 
                                  #xeval = seq(1, 100, 1), 
                                  xeval = seq(0, 100, 0.5), 
                                  bw = bw.red3+7, deg = 2, kernel = EpaK)


ggplot()  +
  geom_point(aes(x = stance_phase, y = pressure), diff[diff$person == "mean", ], 
             size = .3) +
  geom_line(aes(x = stance_phase, y = pressure), diff[diff$person == "mean", ], 
            linewidth = .3, alpha = .6) +
 # geom_line(aes(x = x, y = beta0), est.normal, col = "blue") +
  geom_line(aes(x = x, y = beta0), est.normal.red, col = "orange") +
  geom_line(aes(x = x, y = beta0), est.normal.red2, col = "red") +
  geom_line(aes(x = x, y = beta0), est.normal.red3, col = "yellow") 


##### FNF confidence Bands #####
n = 18
cov.est = crossprod(t(diff.mat[,2:19] - diff.mat[,20]))/((n-1)*n)
#cov.est = est[,,1] passt von den Dimensionen her nicht. 
hat.tau = tau_fun(diff.mat[,2:19]) 

b = confidence_band(x = est.normal$beta0, cov = cov.est, tau = hat.tau,
                    df = n-1, type = "Bs", conf.level = 0.90, 
                    n_int = 8)
est.normal$lo = b[,2]
est.normal$up = b[,3]

temp.df = rbind(
  data.frame(x = est.normal$x, pressure = est.normal$beta0, p = factor("200")),
  data.frame(x = est.normal.red$x, pressure = est.normal.red$beta0, p = factor("50")),
  #data.frame(x = est.normal.red2$x, pressure = est.normal.red2$beta0, p = factor("40")),
  data.frame(x = est.normal.red3$x, pressure =est.normal.red3$beta0, p = factor("20"))
)

temp.df

ggplot() + 
  geom_point(aes( x = x, y = pressure, group = type), size = .8,
             data = data.frame(x = diff.mat[,1], pressure = diff.mat[,20], type = factor("mean"))) + 
  geom_ribbon(aes(x = x, ymin = lo, ymax = up), est.normal, col = "grey", alpha = .4) +
  geom_line(aes(x = x, y = pressure, col = p, lty = p), data = temp.df, lwd = .7) + 
  labs(x = "% of stance phase", y = "Nm/Kg", subtitle = "Mean difference estimate") +
  scale_color_manual(values=c("red", "green", "blue")) + 
  deriv_est_theme


ggsave("biomechanics/grafics/discretization_effect_biomechanics.pdf", device = "pdf", width = 5, height = 3.8, units = "in")

