## from: Reporting FDR analogous confidence intervals for the log fold change of differentially expressed genes
## doi: 10.1186/1471-2105-12-288
pacman::p_load(
  limma,
  mvtnorm,
  tidyverse,
  glue
)
##################################
### Artificial microarray data ###
##################################
d = 500
n = 10
mu1 = rep(0, d)
mu2 = rep(0, d)
tau = 0.25
deg = sample(1:d, tau * d, replace=FALSE)
mu2[deg] = 1
Sigma = rbind(
  c(0.80, 0.60, 0.40, 0.20, 0.00),
  c(0.60, 0.80, 0.60, 0.40, 0.20),
  c(0.40, 0.60, 0.80, 0.60, 0.40),
  c(0.20, 0.40, 0.60, 0.80, 0.60),
  c(0.00, 0.20, 0.40, 0.60, 0.80))
Sigma = Sigma %x% matrix(1, d/5, d/5)
diag(Sigma) = seq(1, 2, length.out=d)
X1 = t(rmvnorm(n, mu1, Sigma))
X2 = t(rmvnorm(n, mu2, Sigma))
X = cbind(X1, X2)
#############
### Limma ###
#############
group = gl(2, n)
design = model.matrix(~ group)
fit1 = lmFit(X, design)
fit = eBayes(fit1)
p.un = fit$p.value[,2]
p.bh = p.adjust(p.un, method="BH")
p.by = p.adjust(p.un, method="BY")
#######################################
### Unadjusted confidence intervals ###
#######################################
alpha = 0.05
beta = fit$coefficients[,2]
std = sqrt(fit$s2.post) * sqrt(fit$cov.coefficients[2,2])
dof = fit$df.prior + fit$df.residual[1]
cl.un = 1 - alpha / 2
lower.un = beta - qt(cl.un, dof) * std
upper.un = beta + qt(cl.un, dof) * std
########################################
### BH-adjusted confidence intervals ###
########################################
R.deg = length(which(p.bh < alpha))
cl = 1 - (R.deg * alpha / d) / 2
lower.bh = beta - qt(cl, dof) * std
upper.bh = beta + qt(cl, dof) * std
########################################
### BY-adjusted confidence intervals ###
########################################
R.deg = length(which(p.by < alpha))
cl = 1 - (R.deg * alpha / d) / 2
lower.by = beta - qt(cl, dof) * std
upper.by = beta + qt(cl, dof) * std

