#
## compute TSNE_plot
#
pacman::p_load(
  tidyverse,
  glue,
  viridis
)
devtools::install_github("jkrijthe/Rtsne", ref = 'openmp')
library(Rtsne)
outdir <- '75_flu_vs_ai_v2/out'
tsne_markers <- read_csv('75_flu_vs_ai_v2/out/predmarkers.csv') %>% 
  filter(fcs_name != 'CD45RA')
  

get_tsne <- function(sdf, SEED = 1234, perplexity = 30) {
  print(glue('computing tSNE for {nrow(sdf)} cells'))
  sdf %>% 
    select(-unique_name, -sample_category, -sample_group, -cell_n) %>% 
    as.matrix() %>% 
    identity() ->
    M
  set.seed(SEED)
  print(head(M))
  tsne_out <- Rtsne(M, check_duplicates = FALSE, pca = FALSE, num_threads = 6, perplexity = perplexity)
  
  sdf %>% 
    dplyr::select(cell_n) %>% 
    mutate(tSNE1 = tsne_out$Y[, 1], 
           tSNE2 = tsne_out$Y[, 2])
}
adf <- read_csv('75_flu_vs_ai_v2/out/tsne_input_subsampled_activated_tcells.csv.gz')
samples <- read_csv('75_flu_vs_ai_v2/out/samples.csv')

adf %>% 
  # dplyr::count(unique_name, sample_group, sample_category) %>% 
  # arrange(sample_group) %>% 
  # print(n = Inf)
  get_tsne(perplexity = 50) %>% 
  identity() ->
  tsne_tweak_df
write_csv(tsne_tweak_df, file.path(outdir, 'tsne_idx.csv'))

adf %>% 
  filter(sample_category == 'full') %>% 
  distinct(sample_group, unique_name) %>% 
  count(sample_group)

tsne_tweak_df %>%
  inner_join(adf) %>% 
  identity() ->
  xdf
xdf %>% 
  filter(sample_category == 'full') %>% 
  ggplot(aes(tSNE1, tSNE2, color = CD62L)) +
  geom_point() +
  # geom_point(data = filter(xdf, sample_category == 'tetramer+'), color = 'darkorange') +
  geom_density_2d(alpha = 0.4) +
  facet_wrap(~ sample_group, ncol = 2) +
  theme_minimal() +
  scale_color_viridis_c()
ggsave('out/ai_vs_flu/tSNE_by_group.png', height = 12)

xdf %>% 
  gather(marker, expression, -unique_name, -sample_group, -tSNE1, -tSNE2, -cell_n, -sample_category) %>% 
  group_by(marker) %>% 
  mutate(frac_of_max = expression / max(expression)) %>% 
  ungroup() %>% 
  mutate(marker = fct_relevel(marker, tsne_markers$fcs_name)) %>% 
  # filter(str_detect(marker, 'PD')) %>%
  # head(20000) %>%
  ggplot(aes(tSNE1, tSNE2)) +
  # geom_point(aes(color = scaled_intensity)) +
  # geom_point(aes(color = scaled_intensity), size = 1) +
  stat_summary_2d(aes(z = frac_of_max), fun = median, bins = 70) +
  # annotate("rect", xmin = tlim$PBMC$xmin, xmax = tlim$PBMC$xmax, 
  #            ymin = tlim$PBMC$ymin, ymax = tlim$PBMC$ymax,
  #            alpha = .1,
  #          color = 'black') +
  facet_wrap(~ marker, ncol = 6) +
  theme_minimal() +
  # scale_fill_viridis() +
  scale_fill_continuous(type = 'viridis', limits = c(0, 1)) +
  labs(fill = 'Scaled Intensity') +
  coord_fixed() %>% 
  identity() ->
  tsne_pbmc_by_marker
tsne_pbmc_by_marker
ggsave(file.path(outdir, 'tsne_by_marker.png'), tsne_pbmc_by_marker, width = 10, height = 10)
ggsave(file.path(outdir, 'tsne_by_marker.pdf'), tsne_pbmc_by_marker, width = 10, height = 10)
