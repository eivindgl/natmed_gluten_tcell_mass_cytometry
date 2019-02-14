pacman::p_load(
  tidyverse,
  glue,
  viridis
)

outdir <- '75_flu_vs_ai_v2/out/clustering'
adf <- read_csv('75_flu_vs_ai_v2/out/clust_subsampled_activated_tcells.csv.gz')
raw_clust_df <- read_csv('75_flu_vs_ai_v2/out/clustering/seed_1234567/cluster_c30_seed1234567.csv.gz')
# raw_clust_df <- read_csv('75_flu_vs_ai_v2/out/clustering/seed_123456/cluster_c30_seed123456.csv.gz')
tsne_idx <- read_csv('75_flu_vs_ai_v2/out/tsne_idx.csv')
samples <- read_csv('75_flu_vs_ai_v2/out/samples.csv')

cdf <- adf %>% 
  inner_join(tsne_idx) %>% 
  inner_join(raw_clust_df)

##
## 1. Plot each cluster faceted by sample_group
##
plots <- cdf %>%
  # filter(sample_category == 'full') %>% 
  mutate(sample_group = if_else(sample_category == 'tetramer+', str_c(sample_group, '_tet+'), sample_group)) %>% 
  group_by(cell_clustering) %>% 
  do(
    plot = ggplot(data = .) + aes(x = tSNE1, y = tSNE2, color = sample_group) +
      geom_point(size = 0.5) +
      facet_wrap(~ sample_group, ncol = 2)
  ) 
  
tsne_clust_dir <- file.path(outdir, 'tsne_clust_facet')
dir.create(tsne_clust_dir, showWarnings = F, recursive = T)
plots %>% 
  mutate(filename = glue('cluster_{cell_clustering}.png'),
         outpath = file.path(tsne_clust_dir, filename)) %>% 
  rowwise() %>% 
  {walk2(.$plot, .$outpath, function(x,y) {
    print(glue('saving to {y}'))
    ggsave(y, x)}
  )}

##
## 2. Plot Markers enriched per cluster
##
adf %>% 
  inner_join(raw_clust_df) %>% 
  filter(sample_category == 'full') %>% 
  select(-unique_name, -cell_n, -sample_category, -sample_group) %>% 
  gather(marker, expr, -cell_clustering) %>% 
  group_by(cell_clustering, marker) %>% 
  summarize(expr = median(expr)) %>% 
  ungroup() %>% 
  ggplot(aes(marker, factor(cell_clustering), fill = expr)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7))
ggsave(file.path(outdir, 'clusters_all.png'))

##
## 3. Plot ratio per subgroup per cluster
##
adf %>% 
  inner_join(raw_clust_df) %>% 
  filter(sample_category == 'full' | sample_group %in% c('ccd', 'ucd')) %>% 
  mutate(sample_group = if_else(sample_category == 'full', sample_group, str_c(sample_group, '_tet+'))) %>% 
  dplyr::count(sample_group, cell_clustering) %>% 
  add_count(sample_group, wt = n) %>% 
  dplyr::rename(clust_n = n, tot_n = nn) %>% 
  mutate(ratio = clust_n / tot_n) %>% 
  complete(sample_group, cell_clustering, fill = list(ratio = 0)) %>% 
  ggplot(aes(sample_group, cell_clustering)) +
  geom_tile(aes(fill = ratio)) +
  scale_fill_continuous(type = 'viridis', limits = c(0, 0.2), oob = scales::squish) +
  theme_minimal() +
  labs(caption = 'ratio capped at 0.2')
ggsave(file.path(outdir, 'clusters_enrichment_by_sample_group.png'))


##
## 4. plot tentative disease specfific clusters
##
cdf %>% 
  # filter(sample_category == 'full') %>% 
  ## Seed 1234567 cluster
  mutate(cluster = case_when(cell_clustering %in% c(9, 13) ~ 'AT',
                             cell_clustering == 6 ~ 'Flu',
                             # cell_clustering == 20 ~ 'Flu_CD57',
                             TRUE ~ 'undef'),
  # mutate(cluster = case_when(cell_clustering %in% c(16, 29) ~ 'AT',
  #                            cell_clustering == 24 ~ 'Flu',
  #                            cell_clustering == 26 ~ 'Flu_CD57',
                             # TRUE ~ 'undef'),
         sample_group = if_else(sample_category == 'tetramer+', str_c(sample_group, '_tet+'), sample_group)) %>% 
  filter(sample_category != 'tetramer+' | sample_group == 'ucd_tet+') %>% 
  mutate(sample_group = fct_relevel(sample_group, 'ucd_tet+', 'ucd', 'hc', 'ccd', 'gfd',  'sle')) %>% 
  identity() ->
  clust_df
  
p_act <- clust_df %>% 
  ggplot(aes(tSNE1, tSNE2)) +
  geom_point(data = filter(clust_df, cluster == 'undef'), color = 'gray', size = 1) +
  geom_point(data = filter(clust_df, cluster != 'undef'), mapping = aes(color = cluster), size = 1) +
  # geom_point(aes(color = cluster), size = 0.6) +
  scale_color_manual(values = c('AT' = '#f36d16', 'Flu' = '#3672bc')) +
  # scale_color_brewer(palette = 'Set1', type = 'qual') +
  facet_wrap(~ sample_group) +
  theme_minimal()
p_act <- p_act +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    legend.position = ''
  )
ggsave(file.path(outdir, 'act_cd4_tsne_viz_flu_and_at.png'), p_act)
ggsave(file.path(outdir, 'act_cd4_tsne_viz_flu_and_at.pdf'), p_act)
pacman::p_load(
  janitor
)

cdf %>% 
  filter(sample_category == 'tetramer+') %>% 
  tabyl(cell_clustering, sample_group) %>% 
  adorn_percentages('col') %>% 
  adorn_totals() %>% 
  adorn_pct_formatting() 

# cdf %>% 
adf %>% 
  inner_join(raw_clust_df) %>% 
  filter(sample_category != 'tetramer+') %>% 
  tabyl(cell_clustering, sample_group) %>% 
  adorn_percentages('col') %>% 
  adorn_totals('row') %>%
  adorn_pct_formatting() 

raw_clust_df %>% 
  inner_join(adf) %>% 
  dplyr::count(cell_clustering, sample_category) %>% 
  spread(sample_category, n, fill = 0) %>% 
  mutate(frac = `tetramer+` / (full + `tetramer+`)) %>% 
  arrange(desc(frac))

##
## Plot frequency of Flu cluster per disease state
##
adf %>% 
  inner_join(raw_clust_df, .) %>% 
  inner_join(samples) %>% 
  # mutate(cluster = case_when(cell_clustering %in% c(16, 29) ~ 'AT',
  #                            cell_clustering == 24 ~ 'Flu',
  #                            cell_clustering == 26 ~ 'Flu_CD57',
  #                            TRUE ~ 'undef')) %>% 
  mutate(cluster = case_when(cell_clustering %in% c(9, 13) ~ 'AT',
                             cell_clustering == 6 ~ 'Flu',
                             # cell_clustering == 20 ~ 'Flu_CD57',
                             TRUE ~ 'undef')) %>%
  filter(sample_category == 'full') %>% 
  dplyr::count(sample_group, unique_name, cluster) %>% 
  add_count(unique_name, wt = n) %>% 
  mutate(freq = n / nn) %>% 
  mutate(sample_group = fct_relevel(sample_group, 'ucd', 'ccd', 'gfd', 'sle', 'ssc', 'flu_inf', 'flu_rec', 'hc')) %>% 
  identity() ->
  freqdf
p <- freqdf %>% 
  filter(cluster != 'undef') %>% 
  ggplot(aes(x=sample_group, y = freq, fill = cluster)) + 
  geom_boxplot(outlier.color = NA, alpha = 0.7) +
  geom_point(aes(color = cluster), position = position_jitterdodge(dodge.width = 0.75, seed = 1234), size = 2.5) +
  scale_y_continuous(labels = scales::percent) + # , limits = c(0, 0.25), oob = scales::squish) +
  scale_fill_manual(values = c('Flu' = '#247b7b', 'AT' = '#80475e')) +
  theme_minimal() +
  labs(y = 'Percentage of activated CD4+ cells', x = 'Disease group')
p
ggsave(file.path(outdir, 'prct_flu_specific_response.png'))
ggsave(file.path(outdir, 'prct_flu_specific_response.pdf'))

p + theme(
  strip.background = element_blank(),
  strip.text.x = element_blank(),
  axis.text = element_blank(),
  axis.title = element_blank(),
  legend.position = '',
  strip.text = element_blank()
)
ggsave(file.path(outdir, 'prct_flu_specific_response_stripped.png'))
ggsave(file.path(outdir, 'prct_flu_specific_response_stripped.pdf'))


freqdf %>% 
  dplyr::rename(n_cluster = n, n_sample = nn) %>% 
  write_csv(file.path(outdir, 'number_of_activated_cells_per_cluster.csv'))


adf %>% 
  inner_join(raw_clust_df, .) %>% 
  inner_join(samples) %>% 
  mutate(cluster = case_when(cell_clustering %in% c(9, 13) ~ 'AT',
                             cell_clustering == 6 ~ 'Flu',
                             # cell_clustering == 20 ~ 'Flu_CD57',
                             TRUE ~ 'undef')) %>% 
  filter(sample_category == 'full') %>% 
  select(-cell_n, -cell_clustering, -sample_category, -disease, -disease_status, -ncells, -donor) %>% 
  gather(marker, expr, -unique_name, -sample_group, -cluster) %>% 
  group_by(sample_group, cluster, marker, unique_name) %>% 
  summarise(expr = mean(expr)) %>% 
  summarise(expr = mean(expr)) %>% 
  identity() ->
  mean_expr_by_clust_and_marker
    
mean_expr_by_clust_and_marker %>% 
  spread(marker, expr)

mean_expr_by_clust_and_marker %>% 
  ggplot(aes(marker, interaction(sample_group, cluster), fill = expr)) +
  geom_tile()+
  scale_fill_viridis_c()

mean_expr_by_clust_and_marker %>% 
  group_by(cluster, marker) %>% 
  summarise(expr = mean(expr)) %>% 
  ungroup() %>% 
  ggplot(aes(marker, cluster, fill = expr)) +
  geom_tile()+
  scale_fill_viridis_c()
