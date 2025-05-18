library(biLocPol) # please install this package from Github first. See "README.md" file for instructions

##### Error Decomposition #####

#' @param W weights
#' @param n amount of curves
#' @param N amount of simulations
#' @param Gamma covariance kernel function form [0,1]^2 to R
#' @param Gamma.arg list with further arguments that shall be passed to Gamma
#' @param lower.order.errors shall the lower order errors be computed?
#' @param f mean function -> not needed for FDA covariance estimation with synchronous design
#' @param r.process distribution of the processes Z
#' @param process.arg further arguments for the function of the process
#' @param r.error distribution of the additional errors
#' @param eps.arg further arguments for the additional errors
#' @param parallel shall the code be executed parallely?
#'
#' @return TODO
#' @importFrom stats rnorm
#' @export
#'
#' @examples
#' 0 # TODO why is this function so large??
error_decomposition = function(W, n, N, Gamma = cov.ou, Gamma.arg = list(theta = 2, sigma = 3), 
                               lower.order.errors = F, f = mu, 
                               r.process = OU, process.arg = list(alpha = 2, sigma = 3), 
                               r.error = rnorm, eps.arg = list(sd = 0.75), 
                               parallel = T, parallel.environment = T){
  if (W$del != 0) {stop("error decomposition for derivatives not implemented (yet)")}
  
  w = W$weights
  p = W$p
  
  up = as.vector(upper.tri(diag(numeric(p))))
  dw = as.vector(lower.tri(diag(numeric(p))))
  
  code = function(useless){
    # simulate data
    # f.eval = f(x.design) # R^p --> not neccessary, since it cancels out anyway
    process.arg = c(list(n = n, t = W$x.design), process.arg)
    Process = do.call(r.process, process.arg) # R^{n x p}
    eps.arg = c(list(n = n*p), eps.arg)
    eps = matrix(do.call(r.error, eps.arg), n, p) # R^{n x p}
    Z = observation_transformation(Process + eps)  # + f.eval if mean function shall be considered
    
    # calculate parts of error decomposition 
    # (10): 1/n * sum_{i = 1}^n e_{i,j} e_{i,k}
    E = rowMeans(apply(eps, 1, tcrossprod))[up] # length p*(p-1)/2
    # T10 = as.vector(crossprod(w, E))
    T10 = eval_weights(W, E) |> abs() |> max()
      
    # (12): 1/n*sum_{i=1}^n (Z_ij Z_ik - Gamma_jk)
    Gamma.w.args = function(x){do.call(Gamma, append(list(t = x), Gamma.arg))}
    G = apply(W$design, 1, Gamma.w.args)
    ZZ = rowMeans(apply(Process, 1, tcrossprod))[up]
    T12 = W |> eval_weights(ZZ - G) |> abs() |> max()
    # T12 = as.vector(crossprod(w, ZZ - G))
    
    # (13): Mixture term: 1/n sum_{i=1}^n Z_ij e_ik + Z_ik e_ij
    temp = rowMeans(sapply(1:n, function(i){ tcrossprod(Process[i,], eps[i,]) }))
    ZE = temp[up] + temp[dw]
    T13 = W |> eval_weights(ZE) |>abs() |> max()
    # T13 = as.vector(crossprod(w, ZE))
    
    # (11): Discretization Error: Gamma_{j,k} - Gamma(x,y)  
    temp = apply(W$eval, 1, Gamma.w.args)
    T11 = max(abs((crossprod(w, G) - temp)))
    # T11 as.vector.data.frame()# T11 = as.vector(crossprod(w, G) - temp)
    
    # Overall Error 
    SUM = max(abs((crossprod(w, Z) - temp)))
    rm(temp)
    
    T14 = 0
    if(lower.order.errors){
      E2 = rowSums(sapply(1:n, function(i){
        crossprod( matrix(eps[i,], n-1, p, T), eps[-i, ] )
      }))[up]
      Z2 = rowSums(sapply(1:n, function(i){
        crossprod( matrix(Process[i,], n-1, p, T), Process[-i, ] )
      }))[up]
      ZE2 = rowSums(sapply(1:n, function(i){
        crossprod( matrix(Process[i,], n-1, p, T), eps[-i, ] )
      }))[up]
      EZ2 = rowSums(sapply(1:n, function(i){
        crossprod( matrix(eps[i,], n-1, p, T), Process[-i, ] )
      }))[up]
      T14 = as.vector(crossprod(w, (E2 + Z2 + EZ2 + ZE2)/(n*(n-1))))
    }
    
    cbind(err.2 = T10, Discr = T11, Process = T12, Mix = T13, Rest = T14, sup.err = SUM)
  }
  
  if(parallel){
    if(parallel.environment)
    {
      cl = makeCluster(detectCores( ) - 1)
      plan(future::cluster)
    }
    erg = future_sapply(1:N, code, future.seed = T)
    if (parallel.environment) { stopCluster(cl) }
  }else{
    erg = sapply(1:N, code)
  }
  rowMeans(erg)
}

##### Bandwidth Evaluation #####

bandwidth_evaluation = function(h, p, p.eval, n.seq, N, 
                                Gamma = min, Gamma.arg = list(), 
                                r.process = BM, process.arg = list(), 
                                r.error = rnorm, eps.arg = list(), m = 1, 
                                grid.type = "less"){
  
  x.design = (1:p - 0.5)/p
  s = length(h)
  x.eval.grid = observation_grid(p.eval, comp = "lesseq")
  up.eval = as.vector(upper.tri(matrix(0, p.eval, p.eval), diag = T))
  temp = as.vector(matrix( (1:p.eval-1/2)/p.eval, p.eval, p.eval, byrow = F))[up.eval]
  
  
  Gamma.w.args = function(x){do.call(Gamma, append(list(t = x), Gamma.arg))}
  G = apply(x.eval.grid, 1, Gamma.w.args)
  
  
  sup.err = numeric(length(n.seq))
  w = local_polynomial_weights(p, h, p.eval, F, m = m, grid.type = grid.type)
  for(k in 1:length(n.seq)){
    n = n.seq[k]
    sup.err[k] =  mean(replicate(N, {
      # f.eval = f(x.design) # R^p --> not neccessary, since it cancels out anyway
      process.arg = c(list(n = n, t = x.design), process.arg)
      Process = do.call(r.process, process.arg) # R^{n x p}
      eps.arg = c(list(n = n*p), eps.arg)
      eps = matrix(do.call(r.error, eps.arg), n, p) # R^{n x p}
      if(grid.type == "less"){
        Z = observation_transformation(Process + eps)
      } else {
        Z = observation_transformation(Process + eps, grid.type = "without diagonal")
      }
      rm(Process); rm(eps)
      apply(crossprod(w$weights, Z) - G, 2, function(m){max(abs(m))})
    }))
  }
  rm(w)
  matrix(c(n.seq, rep(c(p, h), each = length(n.seq)), sup.err), ncol =  4)
}





##### Cross Validation #####
#K-Fold Cross Validation for bivariate local polynomial estimator 
lopocv_sim = function(N, n, p, h.seq, m = 1, w.parallel = F,
                      theta = 2, sigma = 3, sd = 0.75, ...){
  
  
  help = function(h){
    
    w_h = local_polynomial_weights(p = p, h = h, p.eval = p, parallel = w.parallel, m = m, parallel.environment = F,...)
    max_diff = numeric(n)
    
    future_replicate(N, {
      Y = FDA_observation(n, x.design = (1:p - 0.5)/p, f = biLocPol::mu,
                          r.process = OU, process.arg = list(alpha = theta, sigma = sigma, x0 = 0), eps.arg = list(sd = sd))
      Y = apply(Y, 2, function(x){x - mean(x)})
      for(l in 1:n){
        Z_l         = tcrossprod(Y[l,], Y[l,])
        Z_minus_l   = observation_transformation(Y[-l,])
        estimate    = eval_weights(w_h, Z_minus_l)
        max_diff[l] = max( abs(Z_l - estimate)[!as.logical(diag(p))])
      }
      rm(Z_l);rm(Z_minus_l)
      mean(max_diff)
    }, future.seed = T)
  }
  mean_sup = sapply(h.seq, help)
  h.seq[apply(mean_sup, 1, which.min)]
}

k_fold_cv_simulation = function(N, n, p, h.seq, K = 5, m = 1, w.parallel = T, theta = 2, sigma = 3, sd = 0.75){
  
  grp = sample(rep(1:K, ceiling(n/K)), n)
  
  help = function(h){
    w_h = local_polynomial_weights(p, h, p.eval = p, m = m, parallel = w.parallel, parallel.environment = F)
    max_diff = numeric(K)
    
    future_replicate(N, {
      Y = FDA_observation(n, x.design = (1:p - 0.5)/p, f = biLocPol::mu,
                          r.process = OU, process.arg = list(alpha = theta, sigma = sigma, x0 = 0), eps.arg = list(sd = sd))
      for(kk in 1:K){
        test_grp     = matrix(observation_transformation(Y[grp == kk,], grid.type = "full"), p, p)
        train_grp    = observation_transformation(Y[grp != kk,], grid.type = "less")
        pred         = eval_weights(w_h, train_grp)
        max_diff[kk] = max(abs((test_grp - pred)[!as.logical(diag(p))]))
      }
      mean(max_diff)
    }, future.seed = T)
  }
  
  
  mean_sup = sapply(h.seq, help)
  
  h.seq[apply(mean_sup, 1, which.min)]
}

k_fold_cv_simulation_2rv = function(N, n, p, h.seq, K = 5, m = 1, 
                                    w.parallel = T){
  
  grp = sample(rep(1:K, ceiling(n/K)), n)
  
  help = function(h){
    w_h = local_polynomial_weights(p, h, p.eval = p, m = m, 
                                   parallel = w.parallel, 
                                   parallel.environment = F)
    max_diff = numeric(K)
    
    future_replicate(N, {
      Y = FDA_observation(n, x.design = (1:p - 0.5)/p,  
                          f = biLocPol::mu, rprocess = z_two_rv, 
                          eps.arg = list(sd = 0.75))
      for(kk in 1:K){
        test_grp     = matrix(observation_transformation(Y[grp == kk,], grid.type = "full"), p, p)
        train_grp    = observation_transformation(Y[grp != kk,], grid.type = "less")
        pred         = eval_weights(w_h, train_grp)
        max_diff[kk] = max(abs((test_grp - pred)[!as.logical(diag(p))]))
      }
      mean(max_diff)
    }, future.seed = T)
  }
  
  
  mean_sup = sapply(h.seq, help)
  
  h.seq[apply(mean_sup, 1, which.min)]
}


##### Process with smooth diagonal in cov kernel ######
# Z(x) = sqrt(2)/6 * (N^2 - 1) * sin(pi * x)  + 2/3 * (M - 1) * (x - 0.5)

z_two_rv = function(n = 1, p = 100, t = NULL){
  if( is.null(t) ) {
    t = (1:p - 1/2)/p
  }
  N = matrix(rnorm(2*n, 0, 1), n, 2)
  apply(N, 1, function(m){
    2/3 * m[1] * sin(pi * t) + 
      sqrt(2)*2/3 * m[2] * cos(1.25*pi*t)
  }) |> t()
}


cov_z_2rv = function(t){
  4/9*sin(pi*t[1])*sin(pi*t[2]) + 
    8/9*cos(1.25*pi*t[1])*cos(1.25*pi*t[2])
}

partial01_cov_z_2rv = function(t){
  4/9*pi*sin(pi*t[1])*cos(pi*t[2]) -
    10/9*pi*cos(1.25*pi*t[1])*sin(1.25*pi*t[2])
}

partial11_cov_z_2rv = function(t){
  4/9*pi^2*cos(pi*t[1])*cos(pi*t[2]) +
    25/18*pi^2*sin(1.25*pi*t[1])*sin(1.25*pi*t[2])
}



##### Bandwidth Evaluation Function #####

bandwidth_evaluation_derivative_OU = function(h, p, p.eval, n.seq, N, eps.sd = 0.25, boundary_correction = F, correction = 0.1){
  if(!boundary_correction){
    correction = 0
  }
  x.design = (1:p - 0.5)/p
  x.eval.design = (1:p.eval - 0.5)/p.eval
  s = length(h)
  x.eval.grid = observation_grid(p.eval, comp = "full")
  G10 = matrix(apply(x.eval.grid, 1, del10_cov_OU), p.eval, p.eval)
  
  sup.err10 = numeric(length(n.seq))
  w = local_polynomial_weights(p, h, p.eval, F, m = 2, grid.type = "less", del = 1)
  for(k in 1:length(n.seq)){
    n = n.seq[k]
    sup.err10[k] =  mean(replicate(N, {
      Y = biLocPol::OU(n, x.design, sigma = sigma, alpha = theta) + matrix(rnorm(n * p, 0, eps.sd), n, p) # n x p
      Z = observation_transformation(Y)
      trim = which(x.design > correction & x.design < (1-correction))
      trim = which(x.eval.design > correction & x.eval.design < (1-correction))
      estimate = eval_weights(w, Z)[,,3]
      max(abs(estimate - G10)[trim, trim])
    }))
  }
  rm(w)
  matrix(c(n.seq, rep(c(p, h), each = length(n.seq)), sup.err10), ncol =  4)
}


error_decomposition_deriv = function(W, n, N, eps.sd = .75, parallel = T,
                                     parallel.environment = T, correction = F, h_c = 0.1){
  
  w = W$weights
  p = W$p
  
  up = as.vector(upper.tri(diag(numeric(p))))
  dw = as.vector(lower.tri(diag(numeric(p))))
  
  if(correction){
    x.eval = (1:p.eval - 0.5)/p.eval
    trim = which(x.eval > h_c & x.eval < (1-h_c))
  }else{
    trim = T
  }
  
  code = function(useless){
    # simulate data
    # f.eval = f(x.design) # R^p --> not neccessary, since it cancels out anyway
    process = biLocPol::OU(n, (1:p - 0.5)/p, sigma = sigma, alpha = theta)
    eps     = matrix(rnorm(n*p, 0, eps.sd), n, p)
    Z = observation_transformation(process + eps)  # + f.eval if mean function shall be considered
    
    # calculate parts of error decomposition 
    # (10): 1/n * sum_{i = 1}^n e_{i,j} e_{i,k}
    E = rowMeans(apply(eps, 1, tcrossprod))[up] # length p*(p-1)/2
    # T10 = as.vector(crossprod(w, E))
    T10 = eval_weights(W, E)[trim,trim,3] |> abs() |> max()
    
    G   = apply(observation_grid(p, comp = "full"), 1, biLocPol::cov_ou, sigma = 3, theta = 2)
    ZZ  = rowMeans(apply(process, 1, tcrossprod))
    T12 = eval_weights(W, (ZZ - G)[up])[trim,trim,3] |> abs() |> max()
    
    # (13): Mixture term: 1/n sum_{i=1}^n Z_ij e_ik + Z_ik e_ij
    temp = rowMeans(sapply(1:n, function(i){ tcrossprod(process[i,], eps[i,]) }))
    ZE = temp[up] + temp[dw]
    T13 = eval_weights(W, ZE)[trim,trim,3] |>abs() |> max()
    # T13 = as.vector(crossprod(w, ZE))
    
    # (11): Discretization Error: Gamma_{j,k} - Gamma(x,y)  
    temp = apply(biLocPol::observation_grid(p.eval, comp = "full"), 1, del10_cov_OU)
    T11 = max(abs((eval_weights(W, G[up])[trim,trim,3] - matrix(temp, p.eval, p.eval, byrow = F)[trim, trim] )))
    # T11 as.vector.data.frame()# T11 = as.vector(crossprod(w, G) - temp)
    
    # Overall Error 
    SUM = max(abs((eval_weights(W, Z)[trim,trim,3] - matrix(temp, p.eval, p.eval, byrow = F)[trim, trim] )))
    rm(temp)
    
    cbind(err.2 = T10, Discr = T11, Process = T12, Mix = T13, sup.err = SUM)
  }
  
  if(parallel){
    if(parallel.environment)
    {
      cl = makeCluster(detectCores( ) - 1)
      plan(future::cluster)
    }
    erg = future_sapply(1:N, code, future.seed = T)
    if (parallel.environment) { stopCluster(cl) }
  }else{
    erg = sapply(1:N, code)
  }
  rowMeans(erg)
}



#### Derivatives of the OU Cov Kernel ####
del10_cov_OU = function(t, sigma = 3, theta = 2){
  if(t[1] > t[2]){
    return(sigma^2/2 * (exp(-theta*(t[1] - t[2])) + exp(-theta*(t[1]+t[2])) ))
  }
  #if(t[1] == t[2]){return(sigma^2 * exp(-theta*2*t[1]))}
  else{
    return(sigma^2/2 * (-exp(theta*(t[1] - t[2])) + exp(-theta*(t[1]+t[2])) ))
  }
}

del01_cov_OU = function(t, sigma = 3, theta = 2){
  return(del10_cov_OU(c(t[2], t[1]), sigma, theta))
}

