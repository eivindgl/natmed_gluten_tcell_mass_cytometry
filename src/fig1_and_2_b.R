# t-SNE plotting fun
pacman::p_load(
  tidyverse,
  rhdf5,
  glue,
  tictoc,
  assertthat,
  svglite,
  devtools,
  viridis
)
devtools::install_github("jkrijthe/Rtsne", ref = 'openmp')
library(Rtsne)

predmarkers <- read_csv('meta/common_markers_for_pred.csv') %>%
  mutate(fcsmarker = str_replace_all(Antibody, '[ -]', '_'))

read_dataset <- function(unique_name, path = 'out/cytof.h5') {
  assert_that(is_bare_character(unique_name))
  h5_path <- glue('samples/{unique_name}') 
  print(glue("reading {h5_path}"))
  M <- h5read(path, h5_path)
  colnames(M) <- h5readAttributes(path, h5_path)$colnames
  M %>% 
    as_tibble() %>% 
    # gather(marker, value) %>% 
    mutate(filename = unique_name) %>% 
    select(filename , everything())
}

all_samples <- read_csv('meta/meta.csv')
filenames <- c("P3p_TetPos.fcs", "P3p_Pre.fcs", "P3d_TetPos.fcs", "P3d_TetNeg.fcs", 
               "P46d_TetPos.fcs", "P46d_TetNeg.fcs", "P50p_Pre.fcs")
all_samples %>% 
  inner_join(tibble(filename = filenames)) %>% 
  identity() ->
  sne_samples


sne_samples %>% 
  group_by(biosource) %>% 
  do(.$filename %>% 
       map_dfr(read_dataset)
  ) %>% 
  ungroup() %>% 
  identity() ->
  fcs_df

get_tsne <- function(sdf, SEED = 1234, perplexity = 30) {
  print(glue('computing tSNE for {nrow(sdf)} cells'))
  sdf %>% 
    select(-biosource, -filename) %>% 
    as.matrix() %>% 
    identity() ->
    M
  set.seed(SEED)
  tsne_out <- Rtsne(M, check_duplicates = FALSE, pca = FALSE, num_threads = 6, perplexity = perplexity)
  
  sdf %>% 
    mutate(tSNE1 = tsne_out$Y[, 1], 
           tSNE2 = tsne_out$Y[, 2])
}

fcs_df %>% 
  filter(biosource == 'SCS') %>% 
  get_tsne(SEED = 99933, perplexity = 90) %>% 
  inner_join(
    dplyr::select(all_samples, filename, disease, disease_status),
    .,
    by = 'filename') %>% 
  identity() ->
  tsne_scs

fcs_df %>% 
  filter(biosource == 'PBMC') %>% 
  get_tsne(SEED = 99933, perplexity = 90) %>% 
  inner_join(
    dplyr::select(all_samples, filename, disease, disease_status),
    .,
    by = 'filename') %>% 
  identity() ->
  tsne_pbmc
# PLOT SCS
dir.create('out/plots', recursive = T, showWarnings = F)
tsne_scs <- all_samples %>% 
  select(filename, sample_category) %>% 
  inner_join(tsne_scs)

tsne_scs %>% 
  ggplot(aes(x = tSNE1, y = tSNE2, color = sample_category)) +
  geom_point(data = filter(tsne_scs, sample_category != 'tetramer+'), size = 1) +
  geom_point(data = filter(tsne_scs, sample_category == 'tetramer+'), size = 1) +
  theme_minimal() +
  scale_color_manual(values = c('tetramer+' = '#f36d16', 'tetramer-' = '#3672bc')) +
  facet_grid(~ disease) +
  coord_fixed() %>% 
  identity() ->
  scs_plot

# PLOT PBMC
tsne_pbmc <- all_samples %>% 
  select(filename, sample_category) %>% 
  inner_join(tsne_pbmc)
tsne_pbmc %>% 
  filter(biosource == 'PBMC') %>% 
  mutate(disease_cat = str_c(disease, sample_category, sep = '_')) %>% 
  # filter(sample_category == 'full' & disease == 'Ced') %>% 
  # plot_group = fct_relevel('Ced_tetramer+', 'Ced_full')) %>% 
  ggplot(aes(x = tSNE1, y = tSNE2, color = sample_category)) +
  geom_point(data = filter(tsne_pbmc, sample_category != 'tetramer+'), size = 1) +
  geom_point(data = filter(tsne_pbmc, sample_category == 'tetramer+'), size = 1) +
  theme_minimal() +
  scale_color_manual(values = c('tetramer+' = '#f36d16', 'full' = '#3672bc')) +
  facet_grid(~ disease) +
  coord_fixed() %>% 
  identity() ->
  pbmc_plot

ggsave('out/plots/fig1b_tsne_scs.png', plot = scs_plot)
ggsave('out/plots/fig1b_tsne_scs.pdf', plot = scs_plot)
ggsave('out/plots/fig2b_tsne_pbmc.png', plot = pbmc_plot)
ggsave('out/plots/fig2b_tsne_pbmc.pdf', plot = pbmc_plot)
