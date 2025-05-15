# Includes Figure 2.7 and 2.8

##### source code #####
source("mean/functions_mean_derivative.R")

##### results #####
load("mean/data/error_decomposition.RData")

#### data from Figure 2.4 ####
load("mean/data/bandwidth_comparison_per_degree.RData")

##### Error Decomp #####
erg_deg_3_df |> 
  group_by(p) |>
  slice_min(sup.err) |>
  pull(h)
which.h = sapply(erg_deg3, which.min)
h = numeric(m)
for(j in 1:m){
  h[j] = H[[j]][which.h[j]]
}
h # optimal h in Table 2.1

set.seed(549)

error_decomp = sapply(1:m, function(j){
  sampleAndDecompositionBB_derivative(n, p[j], h[j], N, deg = deg, 
                                      rm_boundary_effect = T, bound_h = h[j],
                                      sigma = sigma, sigma.bm = 1)
})


error_decomp_deg3 = data.frame(sup.error = as.vector(error_decomp), 
                               p = rep(p, each = 4), 
                               term = rep(c("sup.err", "bias", "eps", "Z"), length(p)))

error_decomp_deg3


##### Figure 2.7 #####
ggplot(error_decomp_deg3, aes(p, sup.error, lty = term, pch = term, col = term)) + 
  geom_point() +
  geom_line() +
  lims(x = c(0, 1000), y= c(0, 1.35)) + 
  labs(subtitle = "n = 600") +
  deriv_est_theme
ggsave("mean/grafics/derivative_error_decomp_n600_cubic.pdf", device = "pdf", width = 5, height = 3.8, units = "in")



##### fixed p; growing n #####
set.seed(65)
p = 300
n_seq = c(10, 50, 100, 200, 400, 800)
h_seq = seq(0.3, 5/200, -0.005)
N = 1000
sigma = .5

# best bw 
x = (0:p)/p
grid = seq(0, 1, 0.01)

multiple_weights = lapply(h_seq, function(h){
  grid = grid[(grid > h) & (grid < (1-h))]
  locPolWeights(x = x , bw = h, deg = 3, 
                xeval = grid, kernel = EpaK)$allWeig[,2,]
})

max_error_multiple_weights = lapply(n_seq, function(n){
  max.errors = sapply(1:N, function(k){
    Z_mean = bm(x, sigma = 1/sqrt(n))
    e_mean = rnorm(p+1, mean=0, sd = sigma/sqrt(n))
    Y = mu(x) + e_mean + Z_mean
    sapply(1:length(h_seq), function(k){
      grid = grid[(grid > h_seq[k]) & (grid < (1-h_seq[k]))]
      max( abs( as.vector( multiple_weights[[k]] %*% Y) - mu_1(grid)))
    })
  }) |> rowMeans()
})

optimal_h = h_seq[max_error_multiple_weights |> sapply(which.min)]

res = sapply(1:length(n_seq), function(i){
  sampleAndDecompositionBB_derivative(n_seq[i], p, optimal_h[i], N, deg = 3, sigma = sigma, rm_boundary_effect = T, bound_h = h_seq[i])
})
df.erg = data.frame(sup.error = as.vector(res), 
                    n = rep(n_seq, each = 4), 
                    term = rep(c("sup.err", "bias", "eps", "Z"), length(p)))
df.erg

ggplot(df.erg, aes(n, sup.error, lty = term, pch = term, col = term)) + 
  geom_point() +
  geom_line() +
  labs(subtitle = "p = 300") + 
  lims(y = c(0.1, 2.9))  + 
  deriv_est_theme

ggsave("mean/grafics/derivative_error_decomp_p300_cubic.pdf", device = "pdf",
       width = 5, height = 3.8, units = "in")
save.image("mean/data/error_decomposition.RData")





