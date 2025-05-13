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


#### ggplot2 theme ####
library(ggplot2)
deriv_est_theme = theme_grey(base_size = 15) + 
  theme(plot.title = element_text(size = 14))


#### target functions #####
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

##### Bandwidth Comparison #####
# Function to perfom a bandwidth selection with known first derivative for multiple sizes 
# of n but with the same p
h_derivative_multiple_n = function(n, p, p.eval, N, h_seq,
                                   sigma = 0.5, sigma.bm = 1, m = 2){
  x = (1:p - 0.5)/p
  x.eval = (1:p.eval - 0.5)/p.eval
  
  
  sup.error = future_sapply(1:length(h_seq), function(i){
    weights = locPolWeights(x, x.eval, m, h_seq[i], EpaK)$allWeig[,2,]
    result_max_error = matrix(0, 2, length(n))
    rownames(result_max_error) = c("rgh", "smth")
    colnames(result_max_error) = n
    for(l in 1:length(n)){
      max.errors = sapply(1:N, function(k){
        mu_x     = f(x)
        Y_rough  = BM(n = 1, t = x,  sigma = 1/sqrt(n[l])) |> as.vector() +
          rnorm(p, mean = 0, sd = sigma/sqrt(n[l])) +
          mu_x
        Y_smooth = biLocPol::z_2rv(n[l], p) |> colMeans() +
          rnorm(p, mean = 0, sd = sigma/sqrt(n[l])) +
          mu_x
        e_rough  = max( abs( as.vector(weights %*% Y_rough)  - f_deriv(x.eval)))
        e_smooth = max( abs( as.vector(weights %*% Y_smooth) - f_deriv(x.eval)))
        c(e_rough, e_smooth)
      })
      result_max_error[,l] = rowMeans(max.errors)
    }
    rm(weights)
    result_max_error
  }, future.seed = T)
  print(p)
  colnames(sup.error) = h_seq
  data.frame(sup.error, n = rep(n, each = 2), process = rep(c("rgh", "smth"), length(n)))
}




# Bandwidth comparison for one single p and n
h_derivative = function(n, p, p.eval, N, h_seq,
                        sigma = .25, sigma.bm = 1, m = 3){
  x = (1:p - 0.5)/p
  x.eval = (1:p.eval - 0.5)/p.eval
  
  sup.error = future_sapply(1:length(h_seq), function(i){
    weights = locPolWeights(x, x.eval, m, h_seq[i], EpaK)$allWeig[,2,]
    max.errors = sapply(1:N, function(k){
      mu_x     = f(x)
      Y_rough  = BM(n = 1, t = x,  sigma = 1/sqrt(n[i])) |> as.vector() +
        rnorm(p, mean = 0, sd = sigma/sqrt(n)) +
        mu_x
      Y_smooth = biLocPol::z_2rv(n, p) |> colMeans() +
        rnorm(p, mean = 0, sd = sigma/sqrt(n)) +
        mu_x
      e_rough  = max( abs( as.vector(weights %*% Y_rough)  - f_deriv(x.eval)))
      e_smooth = max( abs( as.vector(weights %*% Y_smooth) - f_deriv(x.eval)))
      c(e_rough, e_smooth)
    })
    rowMeans(max.errors)
  }, future.seed = T)
  print(p)
  colnames(sup.error) = h_seq
  rownames(sup.error) = c("rgh","smth")
  sup.error
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


sampleAndDecompositionBB_derivative = function(n, p, h, N,
                                               sigma = 1, sigma.bm = 1, 
                                               grid = seq(0,1,0.001), deg = 2, 
                                               rm_boundary_effect = F, 
                                               bound_h = 0.125){
  x=(0:p)/p
  max.error = matrix(0, nrow = N, ncol = 4)
  if(rm_boundary_effect){
    grid = grid[(grid > bound_h) & (grid < (1-bound_h))]
  }
  weights = locPolWeights(x = x, bw = h, deg = deg, xeval = grid, kernel = EpaK)$allWeig[,2,]
  I_1 = as.vector((weights%*%mu(x)) - mu_1(grid))
  for(k in 1:N){
    Z_mean = bm(x, sigma = sigma.bm/sqrt(n))
    e_mean = rnorm(p + 1 , mean = 0, sd = sigma/sqrt(n))
    I_2 = as.vector(weights %*% e_mean)
    I_3 = locWeightsEval(weights, Z_mean)
    max.error[k,] = c(max(abs(I_1 + I_2 + I_3)), 
                      max(abs(I_1)),
                      max(abs(I_2)),
                      max(abs(I_3)))
  }
  return(colMeans(max.error))
}


