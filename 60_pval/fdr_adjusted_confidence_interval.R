## Compute adjusted confidence intervals for log fold change
#
# Based on jung et al paper
## doi: 10.1186/1471-2105-12-288

# # A tibble: 23 x 7
# marker      FC FC_gmean       p_bh lower_bh upper_bh lower_abs
# <chr>    <dbl>    <dbl>      <dbl>    <dbl>    <dbl>     <dbl>
#   1 CD62L  -3.68     -2.88  0.00000122   -4.49    -2.87      2.87 
# 2 CD39    3.57      3.75  0.00000115    2.81     4.33      2.81 
# 3 CD25   -2.70     -1.95  0.0000907    -3.69    -1.71      1.71 
# 4 CD57   -2.68     -3.14  0.000196     -3.78    -1.58      1.58 
# 5 CXCR3   2.62      2.71  0.00000215    1.99     3.25      1.99 
# 6 CD45RA -2.44     -1.82  0.0000215    -3.17    -1.70      1.70 
# 7 CD161   2.20      2.27  0.0000253     1.49     2.90      1.49 
# 8 CD73   -2.03     -1.40  0.000176     -2.85    -1.22      1.22 
# 9 CD38    1.99      2.06  0.0000253     1.35     2.62      1.35 
# 10 ICOS    1.88      1.92  0.0000253     1.29     2.46      1.29 
# 11 CCR7   -1.78     -1.45  0.000176     -2.50    -1.07      1.07 
# 12 CXCR4  -1.66     -0.815 0.00207      -2.59    -0.733     0.733
# 13 CD127  -1.59     -1.31  0.00408      -2.58    -0.603     0.603
# 14 OX40    1.54      1.71  0.000614      0.807    2.27      0.807
# 15 CD49d   1.43      1.46  0.000393      0.792    2.08      0.792
# 16 CCR6    1.09      1.10  0.00694       0.346    1.83      0.346
# 17 CD28    1.02      0.820 0.00408       0.382    1.65      0.382
# 18 CD27   -0.676    -0.794 0.163        -1.68     0.327     0.327
# 19 CD69   -0.632    -0.442 0.133        -1.49     0.222     0.222
# 20 KLRG1  -0.479    -0.724 0.163        -1.20     0.240     0.240
# 21 CCR5   -0.148    -0.684 0.764        -1.15     0.855     0.855
# 22 CXCR5   0.125     0.189 0.663        -0.440    0.690     0.440
# 23 CCR4    0.0542    0.297 0.868        -0.693    0.801     0.693
pacman::p_load(
  limma,
  mvtnorm,
  tidyverse,
  glue,
  statmod,
  assertthat
)

df <- read_csv('out/mean_marker_expression_per_donor.csv') %>% 
  filter(donor != 'CD1570')
df %>% 
  filter(biosource == 'PBMC', disease_category == 'Celiac - Untreated') %>% 
  dplyr::select(-biosource, -disease, -disease_category) %>% 
  mutate(sample_category = if_else(sample_category == 'tetramer+', 'tetramer_pos', 'tetramer_neg')) %>% 
  mutate(sample = glue('{donor}__{sample_category}')) %>% 
  identity() ->
  x
x %>% 
  distinct(donor, sample_category, sample) %>% 
  identity() ->
  pbmc_sample_to_donor_table

x %>% 
  select(-donor, -sample_category) %>% 
  spread(sample, value) %>% 
  identity() ->
  x
M <- select(x, -1) %>% 
  as.matrix() %>% 
  log2()
rownames(M) <- x$variable  
M

assert_that(all(pbmc_sample_to_donor_table$sample == colnames(M)))
# design with donors too complex for the low number of samples (since we are using means)
# possible to do this at the cell level ... 
design <- model.matrix(~ donor + sample_category, data = pbmc_sample_to_donor_table) 
design_simple <- model.matrix(~ sample_category, data = pbmc_sample_to_donor_table)
rescol <- 'sample_categorytetramer_pos'
assert_that(rescol %in% colnames(design))
# fit1 <- lmFit(log2(M), design)
fit1 <- lmFit(M, design)
# fit1 <- lmFit(M, design)
fit = eBayes(fit1)
p.un = fit$p.value[,rescol]
p.bh = p.adjust(p.un, method="BH")
p.by = p.adjust(p.un, method="BY")

#######################################
### Unadjusted confidence intervals ###
#######################################
alpha = 0.05
beta = fit$coefficients[,rescol]
std = sqrt(fit$s2.post) * sqrt(fit$cov.coefficients[rescol,rescol])
dof = fit$df.prior + fit$df.residual[1]
cl.un = 1 - alpha / 2
lower.un = beta - qt(cl.un, dof) * std
upper.un = beta + qt(cl.un, dof) * std
########################################
### BH-adjusted confidence intervals ###
########################################
R.deg = length(which(p.bh < alpha))
cl = 1 - (R.deg * alpha / nrow(M)) / 2
lower.bh = beta - qt(cl, dof) * std
upper.bh = beta + qt(cl, dof) * std
########################################
### BY-adjusted confidence intervals ###
########################################
R.deg = length(which(p.by < alpha))
cl = 1 - (R.deg * alpha / nrow(M)) / 2
lower.by = beta - qt(cl, dof) * std
upper.by = beta + qt(cl, dof) * std

assert_that(all(names(p.bh) == names(lower.bh)))
ldf <- tibble(marker = names(p.bh), p_bh = p.bh, lower_bh = lower.bh, upper_bh = upper.bh, FC=fit$coefficients[, rescol]) %>% 
  mutate(lower_abs = pmin(abs(lower_bh), abs(upper_bh))) %>% 
  arrange(p_bh)


## Testing
gmean <- read_csv('out/grand_mean_marker_fold_change.csv')
gmean %>% 
  filter(variable == 'PD-1')

gmean %>% 
  filter(disease == 'Ced' & biosource == 'PBMC') %>% 
  spread(sample_category, value) %>% 
  mutate(fold_change = log2(`tetramer+` / `tetramer-`)) %>% 
  select(biosource, variable, fold_change) %>% 
  mutate(category = 'log(CeD+ / CeD-)') %>% 
  identity() ->
  fc_tetpos_over_Ced
fc_tetpos_over_Ced %>% 
  dplyr::rename(FC_gmean = fold_change, marker = variable) %>% 
  dplyr::select(marker, FC_gmean) %>% 
  inner_join(ldf) %>% 
  select(marker, FC, FC_gmean, everything()) %>% 
  arrange(desc(abs(FC))) %>% 
  dput()
  print(n = Inf)

df %>% 
  filter(biosource == 'PBMC', disease_category == 'Celiac - Untreated') %>% 
  dplyr::select(-biosource, -disease, -disease_category) %>% 
  spread(sample_category, value) %>% 
  mutate(fold_change = log2(`tetramer+` / `tetramer-`)) %>% 
  filter(variable == 'PD_1')

###
### Limma fold change for PD-1 is larger than for any donor ... this must be wrong or an estimate of something else?
###
ldf %>% 
  print(n = Inf)
