##### source codes #####
source("Mean Simulation/functions_mean.R")

###### n = 600 ######
N = 1000
n = 600
alpha = 2.1
p = c(65, 115, 175, 275, 400, 550) # dis
deg = floor(alpha)

set.seed(244) # 244


m = length(p)
H = sapply(p, function(x)(rev(seq(0.2, 3/x, -0.005))))
H

cl = makeCluster(detectCores( ) - 1);
plan(multisession)
erg = lapply(1:m, function(j){
  print(j)
  return(h_est_BB(n, p[j], N, H[[j]], deg = deg )) 
})
stopCluster(cl);

med_bandwidth_from_cv = c(0.081245, 0.073190, 0.065435, 0.061080, 0.055635, 0.054645)
erg.df = data.frame(sup.err = unlist(erg), p = 
                      factor(unlist(sapply(1:m, function(j){rep(p[j], length(H[[j]]) + 1)}))),
                    h = c(0, H[[1]], 0, H[[2]], 0, H[[3]], 0, H[[ 4]], 0, H[[ 5]], 0, H[[ 6]])
) |>
  mutate(cv_h = unlist(sapply(1:m, function(j){rep(med_bandwidth_from_cv[j], length(H[[j]]) + 1)}))) |> 
  group_by(p) |>
  mutate(min_h = h[which.min(sup.err)]) |>
  ungroup()

# Figure 2.2
erg.df |> 
  ggplot() + 
  geom_point(aes(x = h, y = sup.err, color = p, pch = p)) +
  geom_point(aes(y = 0, x = cv_h, color = p, pch = p)) + 
  lims(y = c(0, 0.3)) +
  labs(subtitle = "n = 600") 

erg.df |>  
  group_by(p) |> 
  slice_min(sup.err)


#ggsave("mean/grafics/optimal_bw_comp_w_interp.png", device = "png", width = 5, height = 3.8, units = "in")



###### n = 100 ######
# not contained in dissertation
n = 100

set.seed(554)

m = length(p)
H = sapply(p, function(x)(rev(seq(0.2, 3/x, -0.005))))
H

cl = makeCluster(detectCores( ) - 1);
plan(cluster);
erg_n100 = lapply(1:m, function(j){
  print(j)
  return(h_est_BB(n, p[j], N, H[[j]], deg = deg, sigma = 0.5)) 
})
stopCluster(cl);


erg_df_n100 = data.frame(sup.err = unlist(erg_n100), p = 
                           factor(unlist(sapply(1:m, function(j){rep(p[j], length(H[[j]]) + 1)}))),
                         h = c(0, H[[1]], 0, H[[2]], 0, H[[3]], 0, H[[ 4]], 0, H[[ 5]], 0, H[[ 6]])
) |>
  group_by(p) |>
  mutate(min_h = h[which.min(sup.err)]) %>% 
  ungroup()

erg_df_n100 |> 
  filter(p %in% c(65, 115, 175, 275, 400, 550)) |>
  ggplot() + 
  geom_point(aes(x = h, y = sup.err, color = p, pch = p)) +
  lims(y = c(0, 0.3)) +
  labs(subtitle = "n = 100") 

#ggsave("mean/grafics/optimal_bw_comp_w_interp_n100.png", device = "png", width = 5, height = 3.8, units = "in")
