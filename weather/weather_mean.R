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
load("weather/data/results_weather_mean.RData")
# includes Figure 2.9 and 2.10

#### Load and process the data ####
load("weather/data/weather_data_nuremberg.RData")
source("mean/functions_mean.R")

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

#### Mean Estimation ####
days = c(1, 4, 8, 12, 15, 18, 22, 25, 29) 
colnames(N) = c("id", "time_stamp", "qn", "pp_10", "TT_10", "tm5_10", "rf_10", "td_10", "eor", "year", "month", "day", "time", "date")

N |> 
  mutate(month = factor(month.name[month], levels = month.name)) |> 
  filter(day %in% days) |> 
  ggplot() +
  geom_line(aes(x = time, y = TT_10, group = year*day, colour = year), alpha = .4) +
  facet_wrap(month ~.)


# number of observations per month
N |> 
  filter(day %in% days, 
         !is.na(TT_10)) |>
  summarise(.by = c(year, month, day), 
            n = n()) |> 
  filter(n == 144) |>  
  summarise(.by = month, 
            m = n()) |>
  arrange(month)

# data reduction, average temperature per time
avg_tt = N |> 
  tibble() |> 
  filter(day %in% days) |> 
  group_by(month, time) |> 
  summarise(mean_tt = mean(TT_10, na.rm = T))
avg_tt 

p = length(avg_tt$time |> unique())
degree = 1 

L = list()
bandwidths = numeric(12)

# save data in n times p matrix
N_mat = N |>
  filter(day %in% days) |> 
  mutate(time = as.character(time)) |> 
  dplyr::select(year, month, day, time, TT_10) |> 
  pivot_wider(names_from = time,
              values_from = TT_10) |> 
  relocate(`00:00:00`, .after = day)
head(N_mat)

##### estimation and fnf conf bands ####
set.seed(43)
for(m in 1:12){
  M = as.matrix(N_mat[N_mat$month == m, -(1:3)]) # n x p 
  M = M[rowSums(is.na(M)) == 0, ]
  print(dim(M))
  
  bw = cv.locpol((0:143)/144, t(M),  h = seq(0.03, 0.15, 0.005), deg = degree)
  #  bw = k_fold_cv_locpol((0:143)/144, t(M),  h = seq(0.02, 0.2, 0.005), deg = degree)
  bandwidths[m] = bw
  print(bw)
  
  L[[m]] = N %>%
    filter(month == m, day %in% days) |> 
    group_by(time) %>%
    summarise_at(vars(TT_10), list(avg_temp = function(x){mean(x, na.rm = T)}))
  L[[m]]$est = locPolSmootherC(x = (0:143)/143, y = colMeans(M, na.rm = T), 
                               xeval = (0:143)/143, bw = bw, 
                               deg = degree, EpaK)$beta0
  L[[m]]$month = m
  
  # CONFIDENCE BANDS
  n = length(M[,1])
  
  cov.est = crossprod( t(t(M) - colMeans(M)) ) / (n * (n-1) )
  hat.tau = tau_fun(t(M)) 
  
  b = confidence_band(x = L[[m]]$est, 
                      cov = cov.est, tau = hat.tau,
                      df = n-1, type = "FFSCB.t", conf.level = 0.90, 
                      n_int = 2)
  L[[m]]$fnf_up = b[,2]
  L[[m]]$fnf_lo = b[,3]
  
}
bandwidths * 60 * 24
temp_estimation = Reduce(rbind, L)
temp_estimation |> 
  ggplot(aes(x = time)) + 
  geom_ribbon(aes(ymin = fnf_lo, ymax = fnf_up), color = "grey", alpha = .3) + 
  geom_line(aes(y = est)) + 
  facet_wrap(month ~.)


class(as.POSIXct(N$time[1:10]))
str(as.POSIXct(N$time[1:10]))

as.POSIXct(N$time[1:10])
as.POSIXct(temp_estimation$time[1:10])


##### Figure 2.9 #####
N |> 
  mutate(month = factor(month.name[month], levels = month.name)) |> 
  mutate(time  = as.POSIXct(time)) |> 
  filter(day %in% days) |> 
  ggplot() +
  geom_line(aes(x = time, y = TT_10, group = interaction(year,day), colour = year), alpha = .4) +
  geom_line(data = temp_estimation |> 
              mutate(month = factor(month.name[month], levels = month.name))|> 
              mutate(time  = as.POSIXct(time)), mapping = aes(y = est, x = time), color = "orange") + 
  facet_wrap(month ~.) +
  labs(y = "°C") +
  scale_x_datetime(date_breaks = "8 hours", date_labels = "%H:%M") + 
  deriv_est_theme


ggsave("weather/grafics/weather_curves.pdf", device = "pdf", unit = "in", height = 6, width = 9)


##### Data reduction #####
###### one-hour interval ######
L = list()
bandwidths = numeric(12)
set.seed(160)

for(m in 1:12){
  M = as.matrix(N_mat[N_mat$month == m, -(1:3)])[, seq(4,144,6)]
  print(dim(M)) # n x p 
  M = M[rowSums(is.na(M)) == 0, ]
  bw = cv.locpol((0.5:23.5)/24, t(M),  h = seq(0.1, 0.2, 0.005), deg = degree)
  bandwidths[m] = bw
  print(bw)
  
  L[[m]] = N |> 
    filter(month == m, day %in% days) |> 
    group_by(time) %>%
    summarise_at(vars(TT_10), list(avg_temp = function(x){mean(x, na.rm = T)}))
  L[[m]]$est = locPolSmootherC(x = (0.5:23.5)/24, 
                               y = colMeans(M, na.rm = T), 
                               xeval = (0:143)/143, bw = bw, 
                               deg = degree, EpaK)$beta0
  L[[m]]$month = m
}

temp_estimation$one_hour_est = Reduce(rbind, L)$est

###### two-hour interval ######
L = list()
bandwidths = numeric(12)
set.seed(199)

for(m in 1:12){
  M = as.matrix(N_mat[N_mat$month == m, -(1:3)])[, seq(7,144,12)]
  print(dim(M)) # n x p 
  M = M[rowSums(is.na(M)) == 0, ]
  bw = 0.15 # cv.locpol((0:11)/12, t(M),  h = seq(0.15, 0.25, 0.005), deg = degree)
  bandwidths[m] = bw
  print(bw)
  
  L[[m]] = N |> 
    filter(month == m, day %in% days) |> 
    group_by(time) %>%
    summarise_at(vars(TT_10), list(avg_temp = function(x){mean(x, na.rm = T)}))
  L[[m]]$est = locPolSmootherC(x = (0.5:11.5)/12, y = colMeans(M, na.rm = T), 
                               xeval = (0:143)/143, bw = bw, 
                               deg = degree, EpaK)$beta0
  L[[m]]$month = m
}

temp_estimation$two_hour_est = Reduce(rbind, L)$est

temp_estimation = temp_estimation |> 
  mutate(month = factor(month.name[month], levels = month.name))
selected_month = c("January", "July")


##### Figure 2.10 ######

temp_estimation = temp_estimation |> 
  mutate(time  = as.POSIXct(time))

ggplot() + 
  geom_point(data = temp_estimation |> 
               filter(month %in% selected_month) , aes(x = time, y = avg_temp), color = "orange", pch = 3, size = .5) +
  geom_point(data = temp_estimation[seq(4, 1728, 6), ] |> 
               filter(month %in% selected_month) , aes(x = time, y = avg_temp), color = "blue", pch = 3, size = .5) +
  geom_point(data = temp_estimation[seq(7, 1728, 12), ] |> 
               filter(month %in% selected_month) , aes(x = time, y = avg_temp), color = "green", pch = 3, size = .5) +
  geom_ribbon(data = temp_estimation |> 
                filter(month %in% selected_month) , aes(x = time, ymin = fnf_lo, ymax = fnf_up), alpha = .4, color = "grey") + 
  geom_line(data = temp_estimation |> 
              filter(month %in% selected_month) , aes(x = time, y = est), color = "orange") + 
  geom_line(data = temp_estimation |> 
              filter(month %in% selected_month) , aes(x = time, y = one_hour_est), color = "blue", linetype = 2) + 
  geom_line(data = temp_estimation |> 
              filter(month %in% selected_month) , aes(x = time, y = two_hour_est), color = "green", linetype = 3) + 
  labs(y = "°C") + 
  facet_wrap(month ~ .) +
  scale_x_datetime(date_breaks = "8 hours", date_labels = "%H:%M") + 
  deriv_est_theme

ggsave("weather/grafics/weather_fnf_jan_july.pdf", device = "pdf", unit = "in", height = 7, width = 9)

temp_estimation |> 
  filter(one_hour_est > fnf_up | one_hour_est < fnf_lo) |> 
  group_by(month) |> 
  summarise(n = n())

temp_estimation |> 
  filter(two_hour_est > fnf_up | two_hour_est < fnf_lo) |> 
  group_by(month) |> 
  summarise(n = n())

save.image("weather/data/results_weather_mean.RData")
