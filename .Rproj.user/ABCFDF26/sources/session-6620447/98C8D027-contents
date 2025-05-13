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


##### data functions #####
mu = function(x){
  x = 2*x - 1
  return(sin(3*pi*x)*exp(-2*(x)^2))
}

mu_1 = function(x){
  x = 2*x - 1
  y= 6*pi * cos(3*pi*x) * exp(-2*x^2) - 
    8*x * sin(3*pi*x) * exp(-2*x^2)
  return(y)
}

bm = function(grid, x0 = 0, sigma = 1){
  p = sum(grid>0)
  var = diff(c(0,grid[grid>0]))
  cumsum(c(x0, rnorm(p, 0, sigma*var^0.5)))
} 

##### Cross Validation #####

cv.locpol = function(x, Y.mat, h = seq(0.03, 0.15, 0.005), deg = 1, ignore.na = T){
  n = length(Y.mat[1,])
  p = length(x)
  
  W = lapply(h, FUN = function(bw){locPolWeights(x, x[-c(1,p)], deg, bw, EpaK)$locWeig})
  sup.err = sapply(1:length(h), function(k){
    tmp = numeric(n)
    for(i in 1:n){
      Y.i = rowMeans(Y.mat[, -i], na.rm = ignore.na)
      Y.i.hat = locWeightsEval(W[[k]], Y.i)
      tmp[i] = max(abs(Y.i.hat-Y.mat[-c(1,p),i]), na.rm = ignore.na)
    }
    mean(tmp)
  })
  
  return(h[which.min(sup.err)])
}


k_fold_cv_locpol = function(x, Y.mat, k = 10, h = seq(0.03, 0.15, 0.005), deg = 1, ignore.na = T){
  n      = length(Y.mat[1,])
  groups = sample(n) |> split(rep(1:k, length.out = n))
  p      = length(x)
  
  W = lapply(h, FUN = function(bw){locPolWeights(x, x[-c(1,p)], deg, bw, EpaK)$locWeig})
  sup.err = sapply(1:length(h), function(r){
    tmp = numeric(n)
    for(i in 1:k){
      Y_minus_k     = rowMeans(Y.mat[, -groups[[i]]], na.rm = ignore.na)
      Y_minus_k_hat = locWeightsEval(W[[r]], Y_minus_k)
      tmp[i]        = max(abs(Y_minus_k_hat - rowMeans(Y.mat[-c(1,p), groups[[i]]])), na.rm = ignore.na)
    }
    mean(tmp)
  })
  
  return(h[which.min(sup.err)])
}
##### Cross validation simulation #####
sim_k_fold_cv = function(n, p, N, h.seq, k = 10, 
                         sigma = 1, sigma.bm = 1, deg = 2, grid = seq(0,1,0.001)){
  
  x=(0:p)/p
  
  max.errors.bw.cv =  future_replicate(N, {
    Z = replicate(n, bm(x, sigma = sigma.bm))
    e = replicate(n, rnorm(p+1, mean = 0, sd = sigma))
    Y.mat = mu(x) + e + Z
    cv.bw = k_fold_cv_locpol(x, Y.mat, h.seq, deg = deg, k = k)
    cv.bw
  }, future.seed = T)
  
  print(paste("p =", p, "done"))
  timestamp()
  
  return(max.errors.bw.cv)
}
sim.lococv = function(n, p, N, h.seq,
                      sigma = 1, sigma.bm = 1, deg = 2, grid = seq(0,1,0.001)){
  
  x=(0:p)/p
  
  max.errors.bw.cv =  future_replicate(N, {
    Z = replicate(n, bm(x, sigma = sigma.bm))
    e = replicate(n, rnorm(p+1, mean = 0, sd = sigma))
    Y.mat = mu(x) + e + Z
    cv.bw = cv.locpol(x, Y.mat, h.seq, deg = deg)
    cv.bw
  }, future.seed = T)
  
  print(paste("p =", p, "done"))
  timestamp()
  
  return(max.errors.bw.cv)
}




##### Bandwidth Selection #####
h_est_BB = function(n, p, N, seq,
                    sigma = 1, sigma.bm = 1, deg = 1, grid=seq(0,1,0.001)){
  
  x=(0:p)/p
  
  ### Locpol
  sup.error = future_sapply(1:length(seq), function(i){
    weights = locPolWeights(x=x, bw=seq[i], deg=deg, 
                            xeval=grid, kernel=EpaK)$locWeig
    max.errors = sapply(1:N, function(k){
      Z_mean = bm(x, sigma = sigma.bm/sqrt(n))
      e_mean = rnorm(p+1, mean=0, sd = sigma/sqrt(n))
      Y = mu(x) + e_mean + Z_mean
      max( abs( as.vector( weights%*%(Y))
                -mu(grid)))
    })
    mean(max.errors)
  }, future.seed = T)
  
  
  ### Interpolation 
  max.errors.interp = sapply(1:N, function(k){
    Z_mean = bm(x, sigma = sigma.bm/sqrt(n))
    e_mean = rnorm(p+1, mean=0, sd = sigma/sqrt(n))
    Y = mu(x) + e_mean + Z_mean
    max(abs( splinefun(x, Y)(grid) - mu(grid) ))
  })
  
  sup.error.interp = mean(max.errors.interp)
  print("done")
  return(c(sup.error.interp, sup.error))
}

##### Bandwidth comparison for derivative estimation #####

h_est_BB_derivative = function(n, p, N, seq,
                    sigma = 1, sigma.bm = 1, deg = 2, grid=seq(0,1,0.001), 
                    rm_boundary_effect = F){
  
  x=(0:p)/p
  
  ### Locpol
  sup.error = future_sapply(1:length(seq), function(i){
    if(rm_boundary_effect){
      grid = grid[(grid > seq[i]) & (grid < (1-seq[i]))]
    }
    weights = locPolWeights(x=x, bw=seq[i], deg=deg, 
                            xeval=grid, kernel=EpaK)$allWeig[,2,]
    max.errors = sapply(1:N, function(k){
      Z_mean = bm(x, sigma = sigma.bm/sqrt(n))
      e_mean = rnorm(p+1, mean=0, sd = sigma/sqrt(n))
      Y = mu(x) + e_mean + Z_mean
      max( abs( as.vector( weights%*%(Y))
                -mu_1(grid)))
    })
    mean(max.errors)
  }, future.seed = T)
  
  
  print("done")
  return(sup.error)
}

h_est_BB.all = function(n, p, N, h,
                        sigma = 1, sigma.bm = 1, deg = 1, grid=seq(0,1,0.001)){
  
  max.errors = rep(0,N)
  max.errors.interp = rep(0,N)
  max.errors.spl = rep(0,N)
  
  x=(0:p)/p
  weights=locPolWeights(x=x, bw=h, deg=deg, 
                        xeval=grid, kernel=EpaK)$locWeig
  
  for(k in 1:N){
    Z_mean = bm(x, sigma = sigma.bm/sqrt(n))
    e_mean = rnorm(p+1, mean=0, sd = sigma/sqrt(n))
    Y = mu(x) + e_mean + Z_mean
    max.errors[k] = max( abs( as.vector( weights%*%(Y))
                              -mu(grid)))
    max.errors.interp[k] = max(abs( splinefun(x, Y)(grid) - mu(grid) ))
    sm.spl = smooth.spline(x, Y, all.knots = T, keep.data = F)
    max.errors.spl[k] = max(abs( predict(sm.spl, grid)$y-mu(grid)))
  }
  
  sup.error = mean(max.errors)
  sup.error.interp = mean(max.errors.interp)
  sup.error.spl = mean(max.errors.spl)
  
  return(c(sup.error, sup.error.interp, sup.error.spl))
}


##### Error Decomposition #####

sampleAndDecompositionBB = function(n, p, h, N,
                                    sigma = 1, sigma.bm = 1, 
                                    grid = seq(0,1,0.001), deg = 1){
  x=(0:p)/p
  max.error = matrix(0,nrow=N,ncol=4)
  weights = locPolWeights(x=x,bw=h,deg=deg,xeval=grid,kernel=EpaK)$locWeig
  I_1 = as.vector((weights%*%mu(x))-mu(grid))
  for(k in 1:N){
    Z_mean = bm(x, sigma = sigma.bm/sqrt(n))
    e_mean = rnorm(p+1,mean=0,sd=sigma/sqrt(n))
    I_2 = as.vector(weights%*%e_mean)
    I_3 = locWeightsEval(weights, Z_mean)
    max.error[k,] = c(max(abs(I_1+I_2+I_3)), 
                      max(abs(I_1)),
                      max(abs(I_2)),
                      max(abs(I_3)))
  }
  return(colMeans(max.error))
}




sampleAndDecompositionBB.spline = function(n, p, N, spar,
                                           sigma = 1, sigma.bm = 1, 
                                           grid = seq(0,1,0.001)){
  x=(0:p)/p
  max.error = matrix(0, nrow=N, ncol=4)
  sm.spl1 = smooth.spline(x, mu(x), spar = spar)
  I_1 = predict(sm.spl1, grid)$y - mu(grid) 
  for(k in 1:N){
    Z_mean = bm(x, 0,  sigma = sigma.bm/sqrt(n))
    e_mean = rnorm(p+1, mean=0, sd=sigma/sqrt(n))
    I_2 = predict(smooth.spline(x, e_mean, spar = spar), grid)$y
    I_3 = predict(smooth.spline(x, Z_mean, spar = spar), grid)$y
    max.error[k,] = c(max(abs(I_1+I_2+I_3)), 
                      max(abs(I_1)),
                      max(abs(I_2)),
                      max(abs(I_3)))
  }
  return(colMeans(max.error))
}

sampleAndDecompositionBB.interp = function(n, p, N,
                                           sigma = 1, sigma.bm = 1, 
                                           grid = seq(0,1,0.001)){
  x=(0:p)/p
  max.error = matrix(0, nrow=N, ncol=4)
  I_1 = splinefun(x, mu(x))(grid) - mu(grid) 
  for(k in 1:N){
    Z_mean = bm(x, 0,  sigma = sigma.bm/sqrt(n))
    e_mean = rnorm(p+1, mean=0, sd=sigma/sqrt(n))
    I_2 = splinefun(x, e_mean)(grid)
    I_3 = splinefun(x, Z_mean)(grid)
    max.error[k,] = c(max(abs(I_1+I_2+I_3)), 
                      max(abs(I_1)),
                      max(abs(I_2)),
                      max(abs(I_3)))
  }
  return(colMeans(max.error))
}

spl.smth.para = function(n, p, N,
                         sigma = 1, sigma.bm = 1, 
                         grid = seq(0,1,0.001)){
  x=(0:p)/p
  spar = numeric(N)
  for(k in 1:N){
    Y = mu(x) + bm(x, 0,  sigma = sigma.bm/sqrt(n)) + rnorm(p+1, mean=0, sd=sigma/sqrt(n))
    spar[k] = smooth.spline(x, Y)$spar
  }
  return(mean(spar))
  
}


