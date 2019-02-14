pacman::p_load(
  tidyverse,
  glue,
  viridis
)

adf <- read_csv('75_flu_vs_ai_v2/out/clust_subsampled_activated_tcells.csv.gz')
raw_clust_df <- read_csv('75_flu_vs_ai_v2/out/clustering/seed_1234567/cluster_c30_seed1234567.csv.gz')
# raw_clust_df <- read_csv('75_flu_vs_ai_v2/out/clustering/seed_123456/cluster_c30_seed123456.csv.gz')
tsne_idx <- read_csv('75_flu_vs_ai_v2/out/tsne_idx.csv')
samples <- read_csv('75_flu_vs_ai_v2/out/samples.csv')

cdf <- adf %>% 
  inner_join(tsne_idx) %>% 
  inner_join(raw_clust_df)

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
  filter(sample_group %in% c('ucd_tet+', 'ucd', 'flu_inf', 'flu_rec')) %>% 
  identity() ->
  clust_df

#m2

p_act <- clust_df %>% 
  ggplot(aes(tSNE1, tSNE2)) +
  geom_point(data = filter(clust_df, cluster == 'undef'), color = 'gray', size = 1) +
  geom_point(data = filter(clust_df, cluster != 'undef'), mapping = aes(color = cluster), size = 1) +
  # geom_point(aes(color = cluster), size = 0.6) +
  scale_color_manual(values = c('Flu' = '#247b7b', 'AT' = '#80475e')) +
  # scale_color_manual(values=wes_palette(n=2, name="Cavalcanti1")) +
  # scale_color_brewer(palette = 'Set1', type = 'qual') +
  facet_grid(sample_group ~ .) +
  theme_minimal() +
  coord_fixed()
p_act
p_act <- p_act +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    legend.position = '',
    strip.text = element_blank()
  )

p_act
ggsave('out/paper_figs/ai_vs_flu_tSNE_subset.png', plot = p_act, bg = 'transparent')
ggsave('out/paper_figs/ai_vs_flu_tSNE_subset.pdf', plot = p_act, bg = 'transparent')
