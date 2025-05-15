##### source codes #####
source("mean/functions_mean_derivative.R")

#### simulation results #####
load("mean/data/bandwidth_comparison_per_degree.RData")

##### Bandwidth Comparison: Different degrees  #####
set.seed(264) 
N = 1000
n = 600

p = c(65, 115, 175, 275, 400, 550, 750, 1000)  # setting in dis
sigma = 1                                                          
m = length(p)
H = sapply(p, function(x)(rev(seq(0.3, 4/x, -0.005))))
H[[7]] = H[[7]][1:27] # not included in Figure 2.4
H[[8]] = H[[8]][1:28] # not included in Figure 2.4
H

###### Degree = 2 ######
alpha = 2.1
deg = floor(alpha)

cl = makeCluster(detectCores( ) - 1);
plan(multisession);
erg_deg2 = lapply(1:m, function(j){
  print(j)
  return(h_est_BB_derivative(n, p[j], N, H[[j]], deg = deg, rm_boundary_effect = T, 
                             sigma = sigma )) 
})
stopCluster(cl);

erg_deg_2_df = data.frame(sup.err = unlist(erg_deg2), p = 
                            factor(unlist(sapply(1:m, function(j){rep(p[j], length(H[[j]]))}))),
                          h = unlist(H)
) 

# Local quadratic estimator bandwidth comparison 

Figure24a = erg_deg_2_df |> filter(p %in% c(65, 115, 175, 275, 400, 550)) |> 
  ggplot(aes(x = h, y = sup.err, color = p, pch = p)) + 
  geom_point() + labs(subtitle = "n = 600") + lims(y = c(0.5, 20)) +
  deriv_est_theme
Figure24a

ggsave("mean/grafics/derivative_bandwidth_comparison_quad.pdf", device = "pdf", width = 5, height = 3.8, units = "in")

###### Degree = 3 ######
set.seed(264) # same seed as for degree = 2
alpha = 3.1
deg = floor(alpha)

cl = makeCluster(detectCores( ) - 1);
plan(multisession);
erg_deg3 = lapply(1:m, function(j){
  print(j)
  return(h_est_BB_derivative(n, p[j], N, H[[j]], deg = deg, rm_boundary_effect = T, 
                             sigma = sigma )) 
})
stopCluster(cl);

erg_deg_3_df = data.frame(sup.err = unlist(erg_deg3), p = 
                            factor(unlist(sapply(1:m, function(j){rep(p[j], length(H[[j]]))}))),
                          h = unlist(H)
) 

# Local cubic estimator bandwidth comparison 
Figure24b = erg_deg_3_df |> filter(p %in% c(65, 115, 175, 275, 400, 550)) |> 
  ggplot(aes(x = h, y = sup.err, color = p, pch = p)) + 
  geom_point() + lims(y = c(0.5, 20)) + labs(subtitle = "n = 600") +
  deriv_est_theme
Figure24b 
ggsave("mean/grafics/derivative_bandwidth_comparison_cubic.pdf", device = "pdf", width = 5, height = 3.8, units = "in")

erg_deg_2_df |> group_by(p) |> slice_min(sup.err)# results for table 2.1
erg_deg_3_df |> group_by(p) |> slice_min(sup.err)

save.image("mean/data/bandwidth_comparison_per_degree.RData")

