pacman::p_load(
  tidyverse,
  rhdf5,
  assertthat,
  ggthemes,
  glue
)
all_samples <- read_csv('input_data/renamed/sample_meta_full.csv')

read_dataset <- function(path, h5_path) {
  M <- h5read(path, h5_path)
  colnames(M) <- h5readAttributes(path, h5_path)$colnames
  M
}

tetp_samples <- all_samples %>% 
  filter(sample_category == 'tetramer+' & disease == 'Ced' & biosource == 'PBMC') %>% 
  filter(!str_detect(unique_name, 'TCD1456P|TCD1469P'))

bg_sample <- all_samples %>% 
  filter(str_detect(unique_name, 'BC11_full_N2'))

samples <- bind_rows(
  bg_sample,
  tetp_samples
)

# samples %>% 
samples %>% 
  .$unique_name %>% 
  map(function(unique_name){
    print(glue('Reading {unique_name}'))
    h5_path <- glue("samples/{unique_name}")
    x <- read_dataset('out/cytof.h5', h5_path)
    if (nrow(x) == 0) {
      print(glue("Skipping sample {unique_name} - Data set is empty (0 cells)"))
      return (NULL)
    }
    x
  }) %>% 
  set_names(samples$unique_name) %>% 
  discard(is.null) %>% 
  identity() ->
  xs

n_ai_samples <- xs %>% 
  map(colnames) %>% 
  bind_rows() %>% 
  gather(sample, marker) %>% 
  filter(marker == 'TIGIT') # AI panel has TIGIT - so this should not be present
assert_that(are_equal(0, nrow(n_ai_samples)))

N_bg = 3e4
set.seed(1234)
xs$BC11_full_N2 <- xs$BC11_full_N2[sample(1:N_bg, replace=FALSE), ]

df <- xs %>% 
  map(as_tibble) %>% 
  bind_rows(.id = 'unique_name')

df <- all_samples %>% 
  select(unique_name, donor, disease_status) %>% 
  inner_join(df)

##
## Compute TSNE
##
devtools::install_github("jkrijthe/Rtsne", ref = 'openmp')
library(Rtsne)
get_tsne <- function(sdf, SEED = 1234, perplexity = 30) {
  print(glue('computing tSNE for {nrow(sdf)} cells'))
  sdf %>% 
    dplyr::select_if(is_bare_double) %>% 
    as.matrix() %>% 
    identity() ->
    M
  set.seed(SEED)
  print(head(M))
  tsne_out <- Rtsne(M, check_duplicates = FALSE, pca = FALSE, num_threads = 6, perplexity = perplexity)
  
  sdf %>% 
    mutate(tSNE1 = tsne_out$Y[, 1], 
           tSNE2 = tsne_out$Y[, 2])
}

sdf <- get_tsne(df)
dir.create('out/rechallenge', showWarnings = F)
p <- sdf %>% 
  # filter(donor != 'BC11') %>% 
  ggplot(aes(tSNE1, tSNE2, color = disease_status)) +
  geom_point() +
  theme_minimal()
p
ggsave('out/rechallenge/sketch_tetp_by_disease_status.png')
p + facet_wrap(~ disease_status)
ggsave('out/rechallenge/sketch_tetp_facet_by_disease_status.png')

sdf %>% 
  filter(disease_status == 'Challenge') %>% 
  distinct(donor) %>% 
  semi_join(sdf, .) %>% 
  ggplot(aes(tSNE1, tSNE2, color = disease_status)) +
  geom_point() +
  theme_minimal() +
  facet_grid(donor ~ disease_status)

sdf %>% 
  ggplot(aes(tSNE1, donor)) +
  stat_bin2d(aes(fill = stat(density)), binwidth = c(8,1)) +
  # stat_density(aes(fill = stat(density)), geom = "raster", position = "identity") +
  scale_fill_continuous(type = 'viridis') + #, limits = c(0, 0.3)) +
  facet_wrap(~ disease_status, ncol = 1, scales = 'free_y')
ggsave('out/rechallenge/sketch_tetp_density_by_disease_status.png')
sdf %>% 
  ggplot(aes(tSNE2, donor)) +
  stat_bin2d(aes(fill = stat(density)), binwidth = c(8,1)) +
  # stat_density(aes(fill = stat(density)), geom = "raster", position = "identity") +
  scale_fill_continuous(type = 'viridis') + #, limits = c(0, 0.3)) +
  facet_wrap(~ disease_status, ncol = 1, scales = 'free_y')

sdf %>% 
  filter(tSNE2 > 32) %>% 
  dplyr::count(disease_status, donor) %>% 
  arrange(desc(n))

sdf %>% 
  gather(marker, expression, -unique_name, -donor, -disease_status, -tSNE1, -tSNE2) %>% 
  group_by(marker) %>% 
  mutate(frac_of_max = expression / max(expression)) %>% 
  ungroup() %>% 
  ggplot(aes(tSNE1, tSNE2)) +
  stat_summary_2d(aes(z = frac_of_max), fun = median, bins = 70) +
  facet_wrap(~ marker, ncol = 6) +
  theme_minimal() +
  # scale_fill_viridis() +
  scale_fill_continuous(type = 'viridis', limits = c(0, 1)) +
  labs(fill = 'Scaled Intensity') +
  coord_fixed()
ggsave('out/rechallenge/sketch_marker_intensity.png')

sdf %>% 
  gather(marker, expression, -unique_name, -donor, -disease_status, -tSNE1, -tSNE2) %>% 
  ggplot(aes(tSNE1, tSNE2)) +
  stat_summary_2d(aes(z = expression), fun = median, bins = 70) +
  facet_wrap(~ marker, ncol = 6) +
  theme_minimal() +
  # scale_fill_viridis() +
  scale_fill_continuous(type = 'viridis', limits = c(0, 5), oob = scales::squish) +
  labs(fill = 'Absolute expression') +
  coord_fixed()
ggsave('out/rechallenge/sketch_marker_absolute_intensity.png')
  
df %>% 
  dplyr::count(disease_status, donor) %>% print(n = Inf)

samples %>% 
  dplyr::count(disease_status, donor) %>% print(n = Inf)
