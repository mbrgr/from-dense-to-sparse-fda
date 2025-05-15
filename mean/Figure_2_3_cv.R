##### source codes #####
source("mean/functions_mean_derivative.R")
source("mean/functions_mean.R")
# Figure 2.3 results from the leave-one-curve-out cross validation

#### results ####
load("mean/data/bw_cross_validation.RData")

#### Figure 2.3 - Cross Validation ####
###### LOCOCV ######

N = 1000
n = 600
alpha = 2.1
p = c(65, 115, 175, 275, 400, 550)
deg = floor(alpha)
H = sapply(p, function(x)(rev(seq(0.2, 3/x, -0.005))))
set.seed(89)
cl = makeCluster(detectCores( ) - 1);
plan(multisession);
timestamp()
cv.bw = sapply(1:length(p), function(j){
  sim.lococv(n, p[j], N, H[[j]])
})


stopCluster(cl);

cv.bw %>% dim()


#save.image("mean/data/bw_cross_validation.RData")


cv.bw |> head()
colnames(cv.bw)  = p
cv.bw |> as_tibble() |> 
  pivot_longer(1:6, names_to = "p",values_to = "h") |> 
  mutate(p = factor(p, levels = c(65, 115, 175, 275, 400, 550))) |> 
  ggplot() + 
  geom_boxplot(aes(p, h, color = p), alpha = .9, show.legend = F) + 
  labs(y = "sup.err", subtitle = "n = 600") + 
  lims(y = c(0.00, 0.155)) + 
  deriv_est_theme


ggsave("mean/grafics/cv_bw_comparison.pdf", device = "pdf", width = 5, height = 3.8, units = "in")

cv.bw |> as_tibble() |> 
  pivot_longer(1:6, names_to = "p",values_to = "h") |> 
  mutate(p = factor(p, levels = c(65, 115, 175, 275, 400, 550))) |> 
  group_by(p) |> 
  summarise(avg = mean(h), med = median(h))


#### k-fold cross validation ####
##### results #####
load("mean/data/bw_cross_validation_n100.RData")

N = 1000
n = 600
alpha = 2.1
p = c(65, 115, 175, 275, 400, 550)
deg = floor(alpha)
H = sapply(p, function(x)(rev(seq(0.2, 3/x, -0.005))))
k = 60

set.seed(89)
cl = makeCluster(detectCores( ) - 1);
plan(multisession);
timestamp()
k_fold_cv_bw = sapply(1:length(p), function(j){
  sim_k_fold_cv(n, p[j], N, H[[j]], k = k)
})


k_fold_cv_bw %>% dim()

k_fold_cv_bw |> head()
colnames(k_fold_cv_bw)  = p
k_fold_cv_bw |> as_tibble() |> 
  pivot_longer(1:6, names_to = "p",values_to = "h") |> 
  mutate(p = factor(p, levels = c(65, 115, 175, 275, 400, 550))) |> 
  ggplot() + 
  geom_boxplot(aes(p, h, color = p), alpha = .9, show.legend = F) + 
  labs(y = "sup.err", subtitle = "n = 600") + 
  lims(y = c(0, 0.17)) + 
  deriv_est_theme

ggsave("mean/grafics/k_fold_cv_bw_comparison.pdf", width = 5, height = 3.8, units = "in")

k_fold_cv_bw |> as_tibble() |> 
  pivot_longer(1:6, names_to = "p",values_to = "h") |> 
  mutate(p = factor(p, levels = c(65, 115, 175, 275, 400, 550))) |> 
  group_by(p) |> 
  summarise(avg = mean(h), med = median(h))

## n = 100
n = 100
k = 10

set.seed(989)
cl = makeCluster(detectCores( ) - 1);
plan(multisession);
timestamp()
k_fold_cv_bw_n100 = sapply(1:length(p), function(j){
  sim_k_fold_cv(n, p[j], N, H[[j]], k = k)
})


colnames(k_fold_cv_bw_n100)  = p
k_fold_cv_bw_n100 |> as_tibble() |> 
  pivot_longer(1:6, names_to = "p",values_to = "h") |> 
  mutate(p = factor(p, levels = c(65, 115, 175, 275, 400, 550))) |> 
  ggplot() + 
  geom_boxplot(aes(p, h, color = p), alpha = .9, show.legend = F) + 
  labs(y = "sup.err", subtitle = "n = 600") + 
  lims(y = c(0, 0.2))  + 
  deriv_est_theme
#theme_grey(base_size = 15) + 
# theme(plot.title = element_text(size = 14))

ggsave("mean/grafics/k_fold_cv_bw_comparison_n100.pdf", device = "pdf", width = 5, height = 3.8, units = "in")

k_fold_cv_bw_n100 |> as_tibble() |> 
  pivot_longer(1:6, names_to = "p",values_to = "h") |> 
  mutate(p = factor(p, levels = c(65, 115, 175, 275, 400, 550))) |> 
  group_by(p) |> 
  summarise(avg = mean(h), med = median(h), sd(h))


erg_df_n100 %>% 
  group_by(p) %>% 
  slice_min(h)

save.image("mean/data/bw_cross_validation_n100.RData")
