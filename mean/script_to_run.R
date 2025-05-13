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

# load functions
source("Mean Simulation/functions.R")


##### Bandwidth Selection #####
# read functions.r first


##### Comparison of different degrees for the local polynomial estimator #####
N = 1000
n = 600
p = c(50, 200, 500)
set.seed(1244) # 244

m = length(p)
H = sapply(p, function(x)(rev(seq(0.3, 3/x, -0.005))))
H

cl = makeCluster(detectCores( ) - 1);
plan(cluster);
erg_loc_lin = lapply(1:m, function(j){
  print(j)
  return(h_est_BB(n, p[j], N, H[[j]], deg = 1 )[-1]) 
})
erg_loc_quad = lapply(1:m, function(j){
  print(j)
  return(h_est_BB(n, p[j], N, H[[j]], deg = 2 )[-1]) 
})
erg_loc_cubic = lapply(1:m, function(j){
  print(j)
  return(h_est_BB(n, p[j], N, H[[j]][-1], deg = 3 )[-1]) 
})
erg_loc_4 = lapply(1:m, function(j){
  print(j)
  return(h_est_BB(n, p[j], N, H[[j]][-(1:5)], deg = 4 )[-1]) 
})
stopCluster(cl);


comp_tibble = tibble(sup.err  = c(unlist(erg_loc_lin),
                                  unlist(erg_loc_quad), 
                                  unlist(erg_loc_cubic),
                                  unlist(erg_loc_4)),
       h = c(rep(unlist(H), 2), H[[1]][-1], H[[2]][-1], H[[3]][-1], H[[1]][-(1:5)], H[[2]][-(1:5)], H[[3]][-(1:5)]),
       p = factor(c(rep(p, erg_loc_lin |> sapply(length)), 
                    rep(p, erg_loc_quad |> sapply(length)),
                    rep(p, erg_loc_cubic |> sapply(length)),
                    rep(p, erg_loc_4 |> sapply(length)))),
       est = as.factor(c(rep("lin", length(unlist(erg_loc_lin))),
               rep("quad", length(unlist(erg_loc_quad))),
               rep("cubic", length(unlist(erg_loc_cubic))),
               rep("quartic", length(unlist(erg_loc_4)))  ) ))

comp_tibble |> 
  ggplot(aes(x = h, y= sup.err, col = p,lty = est)) + 
  geom_line() + 
  lims( y = c(0.01, 0.31)) + 
  scale_linetype_manual(values = 2:5, 
                        breaks = c("lin", "quad", "cubic", "quartic")) + 
  labs(subtitle = "n = 600")


ggsave("degree_comparison.png", device = "png", width = 6, height = 3.8, units = "in")

##### Error Decomp #####

which.h = sapply(erg, which.min)
h = numeric(m)
for(j in 1:m){
  h[j] = H[[j]][which.h[j]]
}
h

set.seed(543)
res = sapply(1:m, function(j){
  sampleAndDecompositionBB(n, p[j], h[j], N, deg = 1)
})
df.erg = data.frame(sup.error = as.vector(res), 
                    p = rep(p, each = 4), 
                    term = rep(c("sup.err", "bias", "eps", "Z"), length(p)))
df.erg

ggplot(df.erg, aes(p, sup.error, lty = term, pch = term, col = term)) + 
  geom_point() + geom_line() + 
  labs(subtitle = "n = 600")



#ggave("bandwidth_comparison.pdf", device = "pdf", width = 5, height = 3.8, units = "in")


##### Estimator Comparison #####
N = 1000
n = c(50, 75, 100, 150, 200, 300, 400, 500, 660, 830, 1000)
alpha = 2.1
p = ceiling((n*log(n))^(3/(2*alpha))) # dense setting
# p = c(20, 35, 50, 60, 80, 110, 150, 220, 400)
deg = floor(alpha)

set.seed(937)

h = (log(p)/(n*p))^(1/(2*alpha+1))/4
#h = rep(0.3, length(p))

sup.err = matrix(0, length(n), 3)
colnames(sup.err) = c("LocPol", "Interp", "Spline")

for(i in 1:length(p)){
  sup.err[i,] = h_est_BB.all(n[i], p[i], N, h[i], deg = 2)
}



df = data.frame(LocPol = sup.err[,1] * sqrt(n) , 
                Interp = sup.err[,2] * sqrt(n) , 
                Spline = sup.err[,3] * sqrt(n) ,
                n.p = factor(paste0("n = ", n, ", p = ",p)),n)
df = melt(df, id.vars = c("n.p", "n"), measure.vars = c("LocPol", "Interp", "Spline"), 
          value.name = "sup.err", variable.name = "method") 
df

ggplot(df, aes(n, sup.err, col = method, pch = method)) + geom_point() + lims(y = c(0,4)) +
  labs(y = expression("sup.err"*{}%.%{}*sqrt(n)), 
       subtitle = expression("p =" ~ (n~log(n))^{0.71}) )


err.decomp = matrix(0, length(p), 4)
err.decomp.interp = matrix(0, length(p), 4)
err.decomp.spline = matrix(0, length(p), 4)



set.seed(423)
spar = sapply(1:length(p), function(j){spl.smth.para(n, p[j], N)})

h = c(0.095, 0.095, 0.085, 0.090, 0.085 ,0.080 ,0.075, 0.080, 0.075, 0.085, 0.080, 
      0.075)
h[c(4,7,10,11)] = c(0.085, 0.08, 0.075, 0.075)

p
####
set.seed(675)
for(i in 1:length(p)){
  err.decomp.interp[i, ] = sampleAndDecompositionBB.interp(n, p[i], N)
  err.decomp[i, ] = sampleAndDecompositionBB(n, p[i], h[i], N, deg = 2)
  err.decomp.spline[i, ] = sampleAndDecompositionBB.spline(n, p[i], N, spar[i])
}


df.erg.locpol = data.frame(sup.error = as.vector(t(err.decomp)), 
                           p = rep(p, each = 4), 
                           term = rep(c("sup.err", "bias", "eps", "Z"), length(p)), method = "LocPol")
df.erg.interp = data.frame(sup.error = as.vector(t(err.decomp.interp)), 
                           p = rep(p, each = 4), 
                           term = rep(c("sup.err", "bias", "eps", "Z"), length(p)), method = "Interp")
df.erg.spline = data.frame(sup.error = as.vector(t(err.decomp.spline)), 
                           p = rep(p, each = 4), 
                           term = rep(c("sup.err", "bias", "eps", "Z"), length(p)), method = "Spline")


df.erg = rbind(df.erg.locpol, df.erg.interp, df.erg.spline)

ggplot(df.erg[df.erg$term != "sup.err",], aes(p, sup.error, lty = term, pch = method, col = term)) + 
  geom_point() + geom_line() + labs(subtitle = "n = 600")  
  lims(y = c(0, 0.25)) 

ggsave("supErr_comp_interp_lp_spline-pdf", device = "pdf",
       width = 5, height = 3.8, units = "in")


ggplot(df.erg[df.erg$term == "sup.err",],
       aes(p, sup.error, pch = method, lty = method, col = method)) + 
  geom_point() + geom_line() + 
  labs(subtitle = "n = 600") + 
  lims(y = c(0, 0.4)) 

ggsave("supErrDecomp_comp_interp_lp_spline.pdf", device = "pdf",
       width = 5, height = 3.8, units = "in")

##### Display curve with some paths #####
    
set.seed(2374)
n = 10
x = seq(0, 1, 0.01)
obs = seq(0, 1, by = 0.05)
Y = mu(obs) + sapply(1:n, function(x){bm(obs, 0,  sigma = 1)}) + 
  matrix(rnorm(length(obs)*n, 0, 0.3), length(obs), n)
df = data.frame(x = x, y = mu(x), fac = "mu")
Y.df = melt(data.frame(obs, Y), id.vars = "obs", 
            variable.name = "i", value.name = "Y")
Y.mean = data.frame(obs, Y = rowMeans(Y))

bw = cv.locpol(obs, Y, seq(0.1, 0.4, 0.01))
bw
Y.est = locPolSmootherC(obs, Y.mean$Y, x, bw, 2, EpaK)
Y.est = data.frame(x = x, Y.hat = Y.est$beta0)

