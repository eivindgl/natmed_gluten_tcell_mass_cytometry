pacman::p_load(
  tidyverse
)

df <- read_csv('out/cytof_marker_asinh_mean_per_donor.csv') %>% 
  filter(variable != 'CD244')
adf <- read_csv('out/cytof_marker_mean_per_donor.csv') %>% 
  filter(variable != 'CD244')

df %>% 
  distinct(biosource, disease, donor) %>% 
  dplyr::count(biosource, disease)

df %>% 
  distinct(variable) %>% 
  arrange(variable) %>% 
  print(n = Inf)


df %>%
  filter(disease_status == 'UCD') %>% 
  spread(sample_category, value)

df %>% 
  filter(disease_status == 'UCD') %>% 
  filter(biosource == 'PBMC') %>% 
  spread(sample_category, value) %>% 
  na.omit() %>% 
  group_by(variable) %>% 
  summarize(pval = t.test(`tetramer+`, `tetramer-`, paired = T)$p.value) %>% 
  arrange(pval) %>% 
  mutate(FDR = p.adjust(pval, method = 'BH')) %>% 
  identity() ->
  pval_pbmc
skip_in_gut <- c('dummy')
# skip_in_gut <- c('CD45RA', 'CD62L', 'CD49d', 'Integrin_b7')
df %>% 
  filter(disease_status == 'UCD') %>% 
  filter(biosource == 'SCS') %>% 
  anti_join(tibble(variable = skip_in_gut)) %>% 
  # mutate(value = asinh(value) / 5.0) %>%
  # identity() ->
  # xdf
  spread(sample_category, value) %>% 
  na.omit() %>% 
  group_by(variable) %>% 
  summarize(pval = t.test(`tetramer+`, `tetramer-`, paired = T)$p.value) %>%
  # summarize(pval = wilcox.test(`tetramer+`, `tetramer-`, paired = T)$p.value) %>% 
  arrange(pval) %>% 
  mutate(FDR = p.adjust(pval, method = 'BH')) %>% 
  identity() ->
  pval_scs

pval_scs %>% 
  filter(FDR < 0.05)
pval_pbmc

list(SCS = pval_scs,
     PBMC = pval_pbmc) %>% 
  bind_rows(.id = 'biosource') %>% 
  identity() ->
  pval_df

adf %>% 
  filter(disease_status == 'UCD') %>% 
  spread(sample_category, value) %>% 
  na.omit() %>% 
  group_by(biosource, variable) %>% 
  summarize(tetpos = mean(`tetramer+`), tetneg = mean(`tetramer-`)) %>% 
  mutate(fold_change = log2(tetpos / tetneg)) %>% 
  dplyr::select(biosource, variable, fold_change) %>% 
  inner_join(pval_df) %>% 
  arrange(biosource, FDR) %>% 
  write_tsv('out/cytof_ucd_marker_paired_T_test_pval.tsv')

### tet + vs HC
extract_foldchanges <- function(x, n) {
  x %>% slice(n) %>% .$data %>% magrittr::extract2(1) %>% .$value 
}
df %>% 
  filter(biosource == 'PBMC') %>% 
  filter((disease_status == 'UCD' & sample_category == 'tetramer+') | 
           (disease == 'HC' & sample_category == 'tetramer-')) %>% 
  select(disease, variable, value) %>% 
  group_by(variable, disease) %>% 
  nest() %>% 
  split(., .$variable) %>% 
  map_dfr(function(x)t.test(extract_foldchanges(x, 1), extract_foldchanges(x, 2))$p.value) %>% 
  gather(marker, pval) %>% 
  mutate(FDR = p.adjust(pval)) %>% 
  arrange(FDR) %>% 
  identity() ->
  pval_df_pbmc_hc

df %>% 
  filter(biosource == 'SCS') %>% 
  filter((disease_status == 'UCD' & sample_category == 'tetramer+') | 
           (disease == 'HC' & sample_category == 'tetramer-')) %>% 
  select(disease, variable, value) %>% 
  # ggplot(aes(disease, value)) + geom_jitter(width = 0.2) + facet_wrap(~ variable) + theme_cowplot()
  group_by(variable, disease) %>% 
  nest() %>% 
  split(., .$variable) %>% 
  map_dfr(function(x)t.test(extract_foldchanges(x, 1), extract_foldchanges(x, 2))$p.value) %>% 
  gather(marker, pval) %>% 
  mutate(FDR = p.adjust(pval)) %>% 
  arrange(FDR) %>% 
  identity() ->
  pval_df_scs_hc
      
list(PBMC = pval_df_pbmc_hc,
     SCS = pval_df_scs_hc) %>% 
  bind_rows(.id = 'biosource') %>% 
  identity() ->
  CeD_vs_HC_pval
  
adf %>% 
  filter(disease == 'HC' | disease_status == 'UCD') %>% 
  filter(str_detect(sample_category, 'tetramer')) %>% 
  # distinct(biosource, disease, donor) %>% 
  # count(biosource, disease)
  group_by(biosource, disease, variable) %>% 
  summarize(expression = mean(value)) %>% 
  ungroup() %>% 
  spread(disease, expression) %>% 
  mutate(fold_change = log2(Ced / HC)) %>% 
  dplyr::select(biosource, variable, fold_change) %>% 
  identity() ->
  CeD_vs_HC_fc

CeD_vs_HC_fc %>% 
  full_join(CeD_vs_HC_pval, by = c('biosource', 'variable' = 'marker')) %>% 
  arrange(biosource, FDR) %>% 
  write_csv('out/cytof_ucd_vs_HC_marker_T_test_pval.tsv')


### tet - vs HC
df %>% 
  filter(biosource == 'PBMC' & sample_category == 'tetramer-') %>% 
  filter(disease_status == 'UCD' | disease == 'HC') %>% 
  select(disease, variable, value) %>% 
  group_by(variable, disease) %>% 
  nest() %>% 
  split(., .$variable) %>% 
  map_dfr(function(x)t.test(extract_foldchanges(x, 1), extract_foldchanges(x, 2))$p.value) %>% 
  gather(marker, pval) %>% 
  mutate(FDR = p.adjust(pval)) %>% 
  arrange(FDR) %>% 
  identity() ->
  pval_pbmc_hc_neg

df %>% 
  filter(biosource == 'SCS' & sample_category == 'tetramer-') %>% 
  filter(disease_status == 'UCD' | disease == 'HC') %>% 
  select(disease, variable, value) %>% 
  # ggplot(aes(disease, value)) + geom_jitter(width = 0.2) + facet_wrap(~ variable) + theme_cowplot()
  group_by(variable, disease) %>% 
  nest() %>% 
  split(., .$variable) %>% 
  map_dfr(function(x)t.test(extract_foldchanges(x, 1), extract_foldchanges(x, 2))$p.value) %>% 
  gather(marker, pval) %>% 
  mutate(FDR = p.adjust(pval)) %>% 
  arrange(FDR) %>% 
  identity() ->
  pval_scs_hc_neg
      
list(PBMC = pval_pbmc_hc_neg,
     SCS = pval_scs_hc_neg) %>% 
  bind_rows(.id = 'biosource') %>% 
  identity() ->
  CeDneg_vs_HC_pval
  
adf %>% 
  filter(disease == 'HC' | disease_status == 'UCD') %>% 
  filter(sample_category == 'tetramer-') %>% 
  # distinct(biosource, disease, donor) %>% 
  # count(biosource, disease)
  group_by(biosource, disease, variable) %>% 
  summarize(expression = mean(value)) %>% 
  ungroup() %>% 
  spread(disease, expression) %>% 
  mutate(fold_change = log2(Ced / HC)) %>% 
  dplyr::select(biosource, variable, fold_change) %>% 
  identity() ->
  CeDneg_vs_HC_fc

CeDneg_vs_HC_fc %>% 
  full_join(CeDneg_vs_HC_pval, by = c('biosource', 'variable' = 'marker')) %>% 
  arrange(biosource, FDR) %>% 
  write_csv('out/cytof_ucd-neg_vs_HC_marker_T_test_pval.tsv')
