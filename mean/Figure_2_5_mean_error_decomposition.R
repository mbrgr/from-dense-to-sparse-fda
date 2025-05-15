# Includes Figure 2.5 and 2.6

##### source code #####
source("mean/functions_mean.R")
source("mean/functions_mean_derivative.R")

##### results #####
load("mean/data/mean_error_decomposition.RData")

##### Estimator Comparison #####
N = 1000
alpha = 2.1
p = c(40, 60, 80, 120, 160, 200, 260, 320, 400, 480, 560)
deg = floor(alpha)
n = 600 # use fixed n

err.decomp = matrix(0, length(p), 4)
err.decomp.interp = matrix(0, length(p), 4)
err.decomp.spline = matrix(0, length(p), 4)



set.seed(423)
spar = sapply(1:length(p), function(j){spl.smth.para(n, p[j], N)})

h = c(rep(.09, 3), rep(.085, 2), .08, rep(.075, 3), rep(.07, 2))
rbind(p, h)

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

# Figure 2.5
ggplot(df.erg[df.erg$term != "sup.err",], aes(p, sup.error, lty = term, pch = method, col = term)) + 
  geom_point() + geom_line() + labs(subtitle = "n = 600")  + 
  lims(y = c(0, 0.16)) + 
  deriv_est_theme

ggsave("mean/grafics/supErr_comp_interp_lp_spline.pdf", device = "pdf",
       width = 5, height = 3.8, units = "in")

# Figure 2.6
ggplot(df.erg[df.erg$term == "sup.err",],
       aes(p, sup.error, pch = method, lty = method, col = method)) + 
  geom_point() + geom_line() + 
  labs(subtitle = "n = 600") + 
  lims(y = c(0, 0.16))  + 
  deriv_est_theme

ggsave("mean/grafics/supErrDecomp_comp_interp_lp_spline.pdf", device = "pdf",
       width = 5, height = 3.8, units = "in")

save.image("mean/data/mean_error_decomposition.RData")
