pacman::p_load(
  tidyverse,
  viridis,
  glue
)
outdir <- 'out/ai_vs_flu/paper_figs'
dir.create(outdir, showWarnings = F, recursive = T)
samples <- read_csv('out/ai_vs_flu/samples_used.csv') %>% 
  mutate(sample_group = case_when(
    sample_group == 'PBMC_CCD' ~ 'CeD Challenge',
    sample_group == 'PBMC_UCD' ~ 'CeD Untreated',
    sample_group == 'PBMC_Flu' ~ disease,
    sample_group == 'PBMC_AutoImmune' ~ disease
  ),
  disease = if_else(str_detect(disease, 'Flu'), 'Flu', disease),
  sample_group = if_else(sample_category == 'tetramer+', str_c(sample_group, '_tet+'), sample_group),
  sample_group = tools::toTitleCase((str_replace_all(sample_group, '_', ' '))),
  sample_group = fct_relevel(sample_group, "CeD Challenge", "CeD Challenge Tet+", "CeD Untreated", "CeD Untreated Tet+", 
                             "SLE", "SSc", "Flu Infected", "Flu Recovering")) %>% 
  dplyr::select(-path, -biosource, -disease_status, -disease, -sample)
  
adf <- read_csv('out/ai_vs_flu/activated_tcells.csv.gz') %>% 
  dplyr::select(-sample_group, -sample_category)
clusterdf <- read_csv("out/ai_vs_flu/cluster/seed_123456/cluster_c30_seed123456.csv.gz")

clusterdf %>% janitor::tabyl(cell_clustering) %>% arrange(n)
clusterdf %>% 
  mutate(cell_clustering = as.integer(cell_clustering),
         cluster = case_when(cell_clustering %in% c(1, 2, 5, 7, 8) ~ 'CM',
                             # cell_clustering %in% c(3, 4, 16, 13) ~ 'sharedCM',
                             # cell_clustering %in% c(6, 11) ~ 'sharedCM',
                             # cell_clustering == 1 ~ 'CM CD27- CCR4-',
                             cell_clustering == 24 ~ 'AT',
                             # cell_clustering == 18 ~ 'CeD',
                             # cell_clustering == 5 ~ 'CD28+CXCR3+HLA-DR-',
                             # cell_clustering == 3 ~ 'Exhausted',
                             cell_clustering == 21 ~ 'Flu', # keep it simple, stupid
                             # cell_clustering %in% c(21, 25, 30) ~ 'Flu', # all these clusters are prevalent in flu
                             TRUE ~ 'undef'),
         cell_clustering = factor(cell_clustering)) %>% 
  identity() ->
  clusterdf
adf %>% 
  dplyr::select(cell_n, unique_name) %>% 
  inner_join(clusterdf) %>% 
  dplyr::count(unique_name, cell_clustering) %>% 
  add_count(unique_name, wt = n) %>% 
  dplyr::rename(n_cluster_sample = n,
                tot_sample = nn) %>% 
  mutate(ratio_cluster_sample = n_cluster_sample / tot_sample,
         log_ratio_cluster_sample = log(n_cluster_sample) - log(tot_sample)) %>% 
  identity() ->
  cells_cluster_sample

cells_cluster_sample %>% 
  dplyr::count(cell_clustering, wt = n_cluster_sample) %>% 
  dplyr::rename(tot_cluster = n) %>% 
  add_count(wt = tot_cluster) %>% 
  transmute(cell_clustering, 
            avg_log_ratio = log(tot_cluster) - log(n),
            avg_ratio = tot_cluster / n) %>% 
  identity() ->
  cluster_freq

cluster_freq %>% 
  arrange(desc(avg_log_ratio))

clusterdf %>% 
  inner_join(adf) %>% 
  inner_join(samples) %>% 
  dplyr::count(cell_clustering, sample_category) %>% 
  spread(sample_category, n, fill = 0) %>% 
  mutate(tot = full + `tetramer+`,
         ratio_tetp = `tetramer+` / tot) %>% 
  arrange(desc(ratio_tetp))

tsnedf <- read_csv('out/ai_vs_flu/tsne_tweak_df.csv') %>% 
  dplyr::select(-sample_group, -sample_category)

clusterdf %>% 
  inner_join(tsnedf, by = 'cell_n') %>% 
  mutate(cell_clustering = factor(cell_clustering)) %>% 
  inner_join(samples) %>% 
  filter(sample_category == 'full') %>% 
  identity() ->
  sdf
  
set.seed(SEED)
tsnedf %>% 
  select(cell_n, unique_name, tSNE1, tSNE2) %>% 
  inner_join(clusterdf) %>% 
  select(-cluster) %>% 
  inner_join(samples) %>% 
  filter(sample_category == 'full') %>% 
  group_by(sample_group) %>% 
  sample_n(3228) %>% 
  ungroup() %>% 
  identity() ->
  df
  
plots <- df %>% 
  mutate(cell_clustering = factor(cell_clustering)) %>% 
  group_by(cell_clustering) %>% 
  nest() %>% 
  mutate(plot = map2(data, cell_clustering, ~ggplot(data = .x, mapping = aes(tSNE1, tSNE2)) + 
                       geom_point(mapping = aes(color = sample_group), size = 0.6) +
                       facet_wrap(~ sample_group, ncol = 3) +
                       ggtitle(.y) +
                       theme_minimal()
                       ))
dir.create('out/ai_vs_flu/cluster_plots', showWarnings = F)

plots %>% 
  mutate(outpath = glue('out/ai_vs_flu/cluster_plots/cluster_{cell_clustering}.png')) %>% 
  rowwise() %>% 
  {walk2(.$plot, .$outpath, function(x,y) {
    ggsave(y, x)}
    )}

df %>% 
  group_by(sample_group, cell_clustering) %>% 
  summarize(cluster_ratio= n()) %>% 
  ungroup() %>% 
  group_by(cell_clustering) %>% 
  mutate(tot_cluster = sum(cluster_ratio),
         zscore = scale(cluster_ratio)) %>% 
  ungroup() %>% 
  mutate(cluster_ratio = log(cluster_ratio / tot_cluster)) %>% 
  ungroup() %>% 
  ggplot(aes(sample_group, cell_clustering)) +
  geom_tile(aes(fill = cluster_ratio)) +
  scale_fill_viridis(limits = c(-3, -0.5), oob = scales::squish)

adf %>% 
  inner_join(clusterdf) %>% 
  select(-unique_name, -cell_n, -cluster) %>% 
  gather(marker, expr, -cell_clustering) %>% 
  group_by(cell_clustering, marker) %>% 
  summarize(expr = median(expr)) %>% 
  ungroup() %>% 
  ggplot(aes(marker, factor(cell_clustering), fill = expr)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7))
ggsave('out/ai_vs_flu/clusters_all.png')
adf %>% 
  inner_join(clusterdf) %>% 
  mutate(cell_clustering = if_else(cluster == 'undef', as.character(cell_clustering), cluster)) %>% 
  select(-unique_name, -cell_n, -cluster) %>% 
  gather(marker, expr, -cell_clustering) %>% 
  group_by(cell_clustering, marker) %>% 
  summarize(expr = median(expr)) %>% 
  ungroup() %>% 
  ggplot(aes(marker, factor(cell_clustering), fill = expr)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7))
ggsave('out/ai_vs_flu/clusters_condensed.png')
   



sdf %>% 
  ggplot(aes(tSNE1, tSNE2)) +
  # geom_density_2d() +
  geom_point(data = filter(sdf, cluster == 'undef'), color = 'gray', alpha = 0.4, size = 0.5) +
  geom_point(data = filter(sdf, cluster != 'undef'), mapping = aes(color = cluster), size = 0.6, alpha = 0.9) +
  # geom_density_2d() +
  facet_wrap(~ sample_group, ncol = 2) +
  theme_minimal() +
  scale_color_brewer(palette = 'Set1', type = 'qual')
ggsave(file.path(outdir, 'activated_tsne_main_groups.pdf'))
ggsave(file.path(outdir, 'activated_tsne_main_groups.png'))

clusterdf %>% 
  inner_join(dplyr::select(adf, cell_n, unique_name)) %>% 
  inner_join(samples) %>% 
  dplyr::select(-unique_name) %>% 
  filter(sample_category == 'full') %>% 
  dplyr::count(sample_group, donor, cluster) %>% 
  add_count(donor, wt = n) %>% 
  mutate(freq = n / nn) %>% 
  identity() ->
  freqdf
freqdf %>% 
  filter(cluster != 'undef') %>% 
  ggplot(aes(x=sample_group, y = freq, fill = cluster)) + 
  geom_boxplot() +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.6), oob = scales::squish) +
  scale_fill_brewer(palette = 'Set1', type = 'qual') +
  theme_minimal() 
ggsave(file.path(outdir, 'activated_main_groups_boxplot.pdf'))
ggsave(file.path(outdir, 'activated_main_groups_boxplot.png'))
freqdf %>% 
  dplyr::rename(n_cluster = n, tot_sample = nn) %>% 
  write_csv(file.path(outdir, 'activated_main_groups_counts.csv'))
  
