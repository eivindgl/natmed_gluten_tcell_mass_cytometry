pacman::p_load(
  tidyverse,
  rhdf5
)
all_samples <- read_csv('input_data/renamed/sample_meta_full.csv')

samples_cd1300 <- all_samples %>% 
  filter(donor == 'CD1300' & biosource == 'PBMC') %>% 
  select(-path, -note, -instrument, -category, -sample_group, -sample)

read_dataset <- function(path, h5_path) {
  M <- h5read(path, h5_path)
  colnames(M) <- h5readAttributes(path, h5_path)$colnames
  M
}

read_sample <- function(unique_names, hdf_path = 'out/cytof.h5') {
  unique_names %>% 
    map(function(u) {
    h5_path <- glue("samples/{u}")
    read_dataset(hdf_path, h5_path)
  }) %>% 
    set_names(unique_names)
}

M <- read_sample(samples_cd1300$unique_name)

NCells <- 1e4
xs <- M %>% 
  map(function(x) {
    if (nrow(x) > NCells) {
      set.seed(2)
      x <- x[sample(nrow(x), NCells, replace = F), ]
    }
    x
  }) 
M_class <- rep(names(M), map_int(xs, nrow))
M <- do.call(rbind, xs)

devtools::install_github("jkrijthe/Rtsne", ref = 'openmp')
library(Rtsne)

get_tsne <- function(M, SEED = 1234, perplexity = 30) {
  set.seed(SEED)
  Rtsne(M, check_duplicates = FALSE, pca = FALSE, num_threads = 6, perplexity = perplexity)
}
tsne_raw <- get_tsne(M)
tsne_df <- as_tibble(tsne_raw$Y) %>% 
  set_names(c('x', 'y')) %>% 
  mutate(sample = M_class) %>% 
  bind_cols(as_tibble(M)) %>% 
  inner_join(samples_cd1300, by = c('sample' = 'unique_name'))

tsne_df %>% 
  filter(sample_category == 'full') %>% 
  group_by(sample) %>% 
  sample_n(5000) %>%
  identity() ->
  sub_pre
tsne_df %>% 
  filter(sample_category != 'full') %>% 
  bind_rows(sub_pre) %>% 
  identity() ->
  df

df %>% 
  mutate(sample_category = fct_relevel(sample_category, 'full')) %>% 
  ggplot(aes(x,y, color = sample_category)) +
  geom_point(size = 1, alpha = 0.5) +
  facet_grid(~ disease_status)

df %>% 
  # mutate(sample_category = fct_relevel(sample_category, 'full')) %>% 
  ggplot(aes(x,y)) +
  geom_point(color = 'grey',
             data = filter(df, sample_category != 'tetramer+'),
             size = 1, alpha = 0.5) +
  geom_point(mapping = aes(color = disease_status),
             data = filter(df, sample_category == 'tetramer+'),
             size = 3) +
  scale_color_brewer(type = 'qual', palette = 'Dark2') +
  theme_minimal()
adf %>% 
  filter(sample_category == 'full') %>% 
  group_by(unique_name) %>% 
  sample_n(min_group) %>% 
  ungroup() %>% 
  bind_rows(filter(adf, sample_category == 'tetramer+')) %>% 
  arrange(cell_n) %>% 
  identity() ->
  tsne_indf

tsne_indf %>% 
  distinct(unique_name) %>% 
  inner_join(ai_flu_samples) %>% 
  write_csv('out/ai_vs_flu/samples_used.csv')

tsne_indf %>% 
  get_tsne(perplexity = 50) %>% 
  identity() ->
  tsne_tweak_df
