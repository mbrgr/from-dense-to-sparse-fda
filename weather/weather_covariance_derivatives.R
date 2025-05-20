#### Packages and Theme ####
library(biLocPol)
library(tidyverse)
library(lubridate)
library(hms)
library(interp)
library(reshape2)
library(plotly)
library(locpol)
library(ffscb)

#### results ####
load("weather/data/weather_estimation.RData")


#### Load and process the data ####
load("weather/data/weather_data_nuremberg.RData")
#source("mean/functions_mean.R")

head(N, n = 7)
dim(N)

N0 = N[7:length(N$MESS_DATUM), ] 
which(N$MESS_DATUM == "2012-01-01 00:00:00")
which(N$MESS_DATUM == "2017-01-01 00:00:00")

N1 = N0 |> 
  mutate(UHRZEIT = as.character(UHRZEIT)) |> 
  dplyr::select(JAHR, MONAT, TAG, UHRZEIT, TT_10) |> 
  pivot_wider(names_from = UHRZEIT,
              values_from = TT_10)
# extend data to improve estimation at the boundaries
N2 = cbind(N1[, 1:3], 
           rbind(NA, N1[-length(N1$JAHR), -(1:3)]), 
           N1[, -(1:3)], 
           rbind(N1[-1, -(1:3)], NA))
# remove days that have wrong values by the previous procedure
N3 = N2[!(((N2[,3] %in% 28:29) & (N2[,2] == 2)) |
            (N2[,3] == 1) |
            ((N2[,3] == 31) & (N2[,2] %in% c(1,3,5,7,8,10,12))) |
            ((N2[,3] == 30) & (N2[,2] %in% c(4,6,9,11)))), ]
sum(is.na(N3))

# remove NA entries
to_rm = apply(N3, 1, function(v){any(is.na(v))}) %>% which()
N3 = N3[-to_rm, ]
sum(is.na(N3))

#### (Co)-Variance Estimation ####

p.eval = 144
p = 144

#evaluation set with extended interval to next an previous day
x_test = seq(-1, 2, length.out = 3*p)

W = local_polynomial_weights(3*p, 0.3, p, T, m = 2, del = 1, eval.type = "diagonal", x.design.grid = x_test) # watch out for evaluation

# evaluation per month
g_hat = array(0, dim = c(p.eval, 3, 12))
for ( m in 1:12) {
  Y = N3[N3[,2] == m, -(1:3)] |>  
    observation_transformation(na.rm = T)
  g_hat[,,m] = eval_weights(W, Y)
}
est_per_month_G = g_hat[,1,]   # variance
est_per_month_G10 = g_hat[,2,] # del01 G
est_per_month_G01 = g_hat[,3,] # del10 G

colnames(est_per_month_G) = factor(1:12, labels = 1:12)
colnames(est_per_month_G10) = factor(1:12, labels = 1:12)
colnames(est_per_month_G01) = factor(1:12, labels = 1:12)
est_G = est_per_month_G %>% 
  as_tibble() %>% 
  pivot_longer(cols = 1:12) %>% 
  mutate(x = rep(W$x.eval, each = 12),
         g = rep("G", each =length(x))) 
est_G10 = est_per_month_G10 %>% 
  as_tibble() %>% 
  pivot_longer(cols = 1:12) %>% 
  mutate(x = rep(W$x.eval, each = 12), 
         g = rep("G10", each =length(x))) 
est_G01 = est_per_month_G01 %>% 
  as_tibble() %>% 
  pivot_longer(cols = 1:12) %>% 
  mutate(x = rep(W$x.eval, each = 12), 
         g = rep("G01", each = length(x))) 
est_per_month = rbind(est_G01, est_G10) %>% mutate(g = factor(g, levels = c("G10", "G01")))



est_per_month %>% 
  mutate(name = factor(name, levels = 1:12)) %>% 
  ggplot(aes( x = x, y = value, color = g)) + 
  geom_line() + 
  facet_wrap(name~.)


est_G %>% 
  mutate(name = factor(name, levels = 1:12)) %>% 
  ggplot(aes( x = x, y = value, color = g)) + 
  geom_line() + 
  facet_wrap(name~.)


time = N$time[7:150]
time 

final_est = est_per_month %>% 
  mutate(time = rep(rep(time, each = 12), 2)) %>% 
  mutate(name = factor(name, levels = 1:12)) 

##### Figure 3.16 #####
final_est %>% 
  ggplot(aes( x = time, y = value, color = g, linetype = g)) + 
  facet_wrap(name~.)  + 
  geom_line(linewidth = .9) +
  labs(y = NULL, x = "time", title = expression(Partial~derivatives~of~Gamma~on~the~diagonal) ) + 
  scale_discrete_manual(
    aesthetics = c("color", "linetype"),
    values = c("G10" = 2,"G01" = 4), 
    name = "Deriv.",
    labels = c("G10" = expression(italic(d)^{"(1,0)"}*Gamma), "G01" = expression(italic(d)^{"(0,1)"}*Gamma))
  ) +
  scale_x_continuous(breaks = c(time[25], time[121]))

ggsave("mean/grafics/weather_part_deriv_gamma_diagonal_all_months.pdf", device = "pdf", width = 9, height = 6, units = "in")

#### 3 dimensional plots ####
add = 44
p.eval = 100
p = 144 + 2*add
x.design = (0:(144-1))/144
x.design.extended = c( x.design[ (144-add+1):144 ] - 1, x.design, x.design[1:add] + 1)

W = local_polynomial_weights(p, 0.3, p.eval, T, m = 2, del = 1, 
                             eval.type = "full", x.design.grid = x.design.extended) # watch out for evaluation

g_hat = array(0, dim = c(p.eval, p.eval, 12))
for ( m in 1:12) {
  Y = N3[N3[,2] == m, -(1:3)] |>  
    observation_transformation()
  g_hat[,,m] = eval_weights(W, Y)[,,2]
}


time = (1:p.eval - 0.5)/p.eval * 24


cs2 = lisG10_UPcs2 = list(c(0, 1), c("lightblue", "darkred"))
UP = upper.tri(g_hat[,,3], diag = T)
estimate_03up = g_hat[,,3]
estimate_03dn = g_hat[,,3]
estimate_03dn[UP] = NA
estimate_03up[!UP] = NA
plot_ly() %>% 
  add_surface(x = ~ time, y = ~ time, z = estimate_03up, alpha = 0.9, colorscale = cs2, lighting = list(
    ambient = 0.7, diffuse = 0.8, specular = 0.1, roughness = 0.9
  )) %>% 
  add_surface(x = ~ time, y = ~ time, z = estimate_03dn, alpha = 0.9, colorscale = cs2, lighting = list(
    ambient = 0.7, diffuse = 0.8, specular = 0.1, roughness = 0.9
  )) %>% 
  layout(
    scene = list(
      xaxis = list(title = "time"),
      yaxis = list(title = "time"),
      zaxis = list(title = "")
    )
  )%>%
  layout(
    scene = list(
      aspectratio = list(x = 1, y = 1, z = 1)  # Controls axis scaling
    )
  )
# save: cov01_march.png

UP = upper.tri(g_hat[,,11], diag = T)
estimate_11up = g_hat[,,11]
estimate_11dn = g_hat[,,11]
estimate_11dn[UP] = NA
estimate_11up[!UP] = NA
plot_ly() %>% 
  add_surface(x = ~ time, y = ~ time, z = estimate_11up, alpha = 0.9, lighting = list(
    ambient = 0.7, diffuse = 0.8, specular = 0.1, roughness = 0.9
  )) %>% 
  add_surface(x = ~ time, y = ~ time, z = estimate_11dn, alpha = 0.9, lighting = list(
    ambient = 0.7, diffuse = 0.8, specular = 0.1, roughness = 0.9
  )) %>% 
  layout(
    scene = list(
      xaxis = list(title = "time"),
      yaxis = list(title = "time"),
      zaxis = list(title = "")
    )
  )  %>%
  layout(
    scene = list(
      aspectratio = list(x = 1, y = 1, z = 1)  # Controls axis scaling
    )
  )


#save: cov01_nov.png

deriv_cov_weather_diag = tibble(x = rep(time, 2), 
                                est = c(diag(g_hat[,,2]), diag(g_hat[,,3])), 
                                deriv = gl(2, 72, labels = c("G10", "G01")))
deriv_cov_weather_diag |> 
  ggplot(aes( x = x, y = est, color = deriv, linetype = deriv)) + 
  geom_line(linewidth = .9) +
  deriv_est_theme + 
  labs(y = NULL, x = "time", title = expression(Partial~derivatives~of~Gamma~on~the~diagonal) ) + 
  scale_discrete_manual(
    aesthetics = c("color", "linetype"),
    values = c(2,4), 
    name = "Deriv.",
    labels = c(expression(italic(d)^{"(1,0)"}*Gamma), expression(italic(d)^{"(0,1)"}*Gamma))
  ) 
