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
load("weather/data/results_weather_mean_derivative.RData")
# includes Figure 2.11 and 2.12


#### Load and process the data ####
load("weather/data/weather_data_nuremberg.RData")
my_theme = theme_grey(base_size = 15) + 
  theme(plot.title = element_text(size = 14))

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



##### Derivative Figures #####
p = 144
N_bar = matrix(0, 12, 3 * p)

for(m in 1:12){
  N_bar[m,] = N3[N3[,2] == m,-(1:3)] |> 
    apply(2, mean, na.rm = T)
}


x_temp = seq(0,1,length.out = 145)[-145]
x_test = seq(-1, 2, length.out = 3*p)

farben = c( "#a6cee3", "#03396c",  # Blau-Töne
            "#33a02c", "#66c21f", "#006400",  # Grün-Töne
            "#e31a1c", "#fb9a99", "#990000",  # Rot-Töne
            "#ff7f00", "#ffb300", "#b15928",  # Gelb-Orange-Töne
            "#1f78b4")

weights_mean_derivative = locPolWeights(x_test, x_temp, 3, 0.2, EpaK)$allWeig[,2,]
monthly_weather_deriv = matrix(0, p, 12)

for(m in 1:12){
  Y = N_bar[m,]
  monthly_weather_deriv[, m] = weights_mean_derivative %*% Y
}

dim(monthly_weather_deriv)
colnames(monthly_weather_deriv) = 1:12
deriv_tibble = monthly_weather_deriv |> 
  as_tibble() |> 
  cbind(x_temp) |> 
  pivot_longer(1:12, values_to = "deriv", names_to = "month") |> 
  mutate(month = as.numeric(month) ) 


time = N$UHRZEIT[7:150]
time 



###### Figure 11 ######
cbind(deriv_tibble, rep(time, each = 1728/144)) |> 
  mutate(month = as.factor(month)) |> 
  rename(time = "rep(time, each = 1728/144)") |> 
  mutate(time = as.POSIXct(time)) |>  
  ggplot(aes(x = time, y = deriv, col = month, linetype = month)) +
  geom_line() +
  labs(x = "time", y = NULL, subtitle = "Derivatives of mean temperatur") + 
  scale_colour_manual(values = farben) + 
  scale_linetype_manual(values = c(2,3,4,4,4,1:3,5,5,5,1)) + 
  scale_x_datetime(date_breaks = "8 hours", date_labels = "%H:%M") + 
  my_theme

ggsave("weather/grafics/derivatives_mean_temperature.pdf", device = "pdf",
       width = 5, height = 3.8, units = "in")

# test if the estimation produces reliable results 
rbind(deriv_tibble, deriv_tibble |> mutate(x_temp = x_temp +1)) |> 
  mutate(month = as.factor(month)) |> 
  ggplot(aes(x = x_temp, y = deriv, col = month, linetype = month)) +
  geom_line() +
  labs(x = "time", y = NULL, subtitle = "Derivatives of mean temperatur") + 
  scale_colour_manual(values = farben) + 
  scale_linetype_manual(values = c(2,3,4,4,4,1:3,5,5,5,1)) 

##### Mean Figure #####
weights_mean = locPolWeights(x_test, x_temp, 3, 0.2, EpaK)$allWeig[,1,]
p = length(x_temp)
monthly_weather = matrix(0, p, 12)

for(m in 1:12){
  Y = N_bar[m,]
  monthly_weather[, m] = weights_mean %*% Y
}

dim(monthly_weather)
colnames(monthly_weather) = 1:12
mean_tibble = monthly_weather |> 
  as_tibble() |> 
  cbind(x_temp) |> 
  pivot_longer(1:12, values_to = "temp", names_to = "month") |> 
  mutate(month = as.numeric(month) ) 

###### Figure 12 ######
ggplot() +
  geom_line(data = cbind(mean_tibble, rep(time, each = 1728/144)) |> 
              mutate(month = as.factor(month)) |> 
              rename(time = "rep(time, each = 1728/144)") |> 
              mutate(time = as.POSIXct(time))  , 
            aes(x = time, y = temp, col = month, linetype = month)) +
  labs(x = "time", y = NULL, subtitle = "Estimated daily mean temperature") + 
  scale_colour_manual(values = farben) + 
  scale_linetype_manual(values = c(2,3,4,4,4,1:3,5,5,5,1)) + 
  geom_text(data = data.frame(time = c(hms(0, -15, 0), hms(0, -15, 0), hms(0, -15, 0), hms(0, -15, 0), hms(0, -15, 0), hms(0, -15, 0),
                                    hms(0, -15, 0), hms(0, -15, 0), hms(0, -15, 0), hms(0, -22, 0), hms(0, -22, 0), hms(0, -22, 0)), 
                              y = c(-.1, 0.3, 2.5, 6.6, 10.2, 14.25, 16, 15.5, 11.8, 7.6, 4.2, 1.5), 
                              month = gl(12,1)) |> mutate(time = as.POSIXct(time)), aes(label = month, x = time, y = y), hjust = 0.4) + 
  scale_x_datetime(date_breaks = "8 hours", date_labels = "%H:%M") + 
  my_theme

ggsave("weather/grafics/mean_temperature.pdf", device = "pdf",
       width = 5, height = 3.8, units = "in")


save.image("weather/data/results_weather_mean_derivative.RData")
