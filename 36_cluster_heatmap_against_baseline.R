##
## Heatmap of Cluster specific T-cells against average of all T-cells.
## A little complicated by the fact that have to read all of the cells for each sample and compute an average
## cluster averages are more simple. Read clustered cells and compute grand mean per marker / sample_cluster

pacman::p_load(
  tidyverse,
  glue,
  viridis
)

outdir <- '75_flu_vs_ai_v2/out/clustering/heatmaps'
dir.create(outdir, showWarnings = F)
adf <- read_csv('75_flu_vs_ai_v2/out/clust_subsampled_activated_tcells.csv.gz')
raw_clust_df <- read_csv('75_flu_vs_ai_v2/out/clustering/seed_1234567/cluster_c30_seed1234567.csv.gz')
samples <- read_csv('75_flu_vs_ai_v2/out/samples.csv')
cluster_names <- read_csv('75_flu_vs_ai_v2/out/clustering/seed_1234567/cluster_names.csv')
predmarkers <- read_csv('75_flu_vs_ai_v2/out/predmarkers.csv')
predorder <- read_csv('input_data/original/asbjorn_markers_specific.csv') %>% 
  filter(category == 'CyTOF')

reverse_scale <- function(x) {
  sinh(x) * 5
}

adf %>% 
  filter(sample_category == 'full') %>% 
  inner_join(raw_clust_df) %>% 
  inner_join(cluster_names) %>% 
  filter(cluster != 'undef') %>% 
  select(-sample_category, -sample_group, -cell_n, -cell_clustering) %>% 
  gather(marker, expression, -unique_name, -cluster) %>% 
  mutate(expression = reverse_scale(expression)) %>% 
  group_by(cluster, marker, unique_name) %>% 
  summarize(expression = mean(expression)) %>% 
  summarize(expression = mean(expression)) %>% 
  ungroup() %>% 
  identity() ->
  clust_mean_expr

p_abs <- clust_mean_expr %>% 
  inner_join(predmarkers, by = c('marker' = 'fcs_name')) %>% 
  mutate(antibody = fct_relevel(antibody, rev(predorder$common_name))) %>% 
  ggplot(aes(cluster, antibody, fill = expression)) + 
  geom_tile() + 
  # scale_fill_viridis(limits = c(0, 70), oob = scales::squish) +
  scale_fill_gradient(low = 'grey20', high = 'yellow2', limits = c(0, 75), oob = scales::squish) +
  theme_minimal() +
  theme(legend.position = 'bottom')
p_abs
##
## Next, read average expression for all samples involved in clustering
##
pacman::p_load(
  rhdf5
)
read_dataset <- function(path, h5_path) {
  M <- h5read(path, h5_path)
  colnames(M) <- h5readAttributes(path, h5_path)$colnames
  M
}

alldf <- samples %>% 
  mutate(h5path = glue('samples/{unique_name}')) %>%
  .$h5path %>%
  map(partial(read_dataset, 'out/cytof.h5')) %>%
  map(~ as_tibble(.x[, predmarkers$fcs_name])) %>% 
  set_names(samples$unique_name) %>% 
  bind_rows(.id = 'unique_name') 

baseline <- alldf %>% 
  gather(marker, expression, -unique_name) %>% 
  mutate(expression = reverse_scale(expression)) %>% 
  group_by(marker, unique_name) %>% 
  summarize(expression = mean(expression)) %>% 
  summarize(expression = mean(expression)) 

p_logfc <- baseline %>% 
  dplyr::rename(baseline_expr = expression) %>% 
  inner_join(clust_mean_expr) %>% 
  inner_join(predmarkers, by = c('marker' = 'fcs_name')) %>% 
  mutate(antibody = fct_relevel(antibody, rev(predorder$common_name))) %>% 
  mutate(logfc = log2(expression / baseline_expr)) %>% 
  ggplot(aes(cluster, antibody, fill = logfc)) + 
  geom_tile() + 
  scale_fill_gradient2(low = 'dodgerblue4', mid = 'grey20', high = 'yellow2', limits = c(-4, 4), oob=scales::squish) +
  theme_minimal() +
  theme(legend.position = 'bottom')
p_logfc

##

pacman::p_load(
  cowplot
)
p <- plot_grid(p_abs, p_logfc)
p
ggsave(file.path(outdir, 'cluster_heatmap_abs_and_fc.png'), p)
ggsave(file.path(outdir, 'cluster_heatmap_abs_and_fc.pdf'), p)
