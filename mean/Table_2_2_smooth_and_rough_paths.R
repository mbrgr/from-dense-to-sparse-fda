##### Packages #####
library(ggplot2)
library(reshape2)
library(locpol)
library(interp)
library(stats)
library(future)
library(future.apply)
library(parallel)
library(tidyverse)
library(biLocPol)

#### results ####
load("mean/data/rate_smooth_and_rough_processes.RData")


#### sim ####
set.seed(22)
n = c(10, 20, 40, 80, 160, 240, 480, 800, 1600)
p = 800
p.eval = 50
x = (1:p - 0.5)/p
x.eval = (1:p.eval - 0.5)/p.eval
N = 1000

h = seq(0.34, .1, by = -0.03)

cl = makeCluster(detectCores(logical = F) - 1)
plan(multisession)
res = future_sapply(1:length(n), function(i){
  w  = locPolWeights(x, x.eval, 2,  h[i], EpaK)$allWeig[,2,] |> as.matrix()
  max.errors = sapply(1:N, function(k){
    Y_rough  = biLocPol::BM(n = 1, t = x,  sigma = 1/sqrt(n[i])) |> as.vector() 
    Y_smooth = biLocPol::z_2rv(n[i], p) |> colMeans() 
    e_rough  = max( abs( as.vector(w  %*% Y_rough) ))
    e_smooth = max( abs( as.vector(w %*% Y_smooth) ))
    c(e_rough, e_smooth)
  })
  print(n[i])     
  rowMeans(max.errors)     
}, future.seed = T)
stopCluster(cl);

round(res[1,] * sqrt(n)*h^0.5, 3)
round(res[2,] * sqrt(n), 3)

plot(res[1,] ~ n, type = "l", ylim = c(0, 3))
lines(res[2,] ~n, lty = 2)
lines(3.5/sqrt(n) ~n, col = "red")
lines(2.5*h_rgh^(-0.5)/sqrt(n) ~ n, col = "red")
plot(res[1,]*h_rgh^0.5/res[2,] ~ n, type = "l", ylim = c(0, 1.5))

save.image(("mean/data/rate_smooth_and_rough_processes.RData"))
