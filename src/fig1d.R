pacman::p_load(
  tidyverse,
  assertthat,
  glue,
  rhdf5
)

select <- dplyr::select
filter <- dplyr::filter
count <- dplyr::count

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


## Read antibodies to generate plots for
markers <- read_csv('meta/cytof_markers_specific.csv')

samples <- read_csv('meta/meta.csv') %>% 
  filter(biosource == 'SCS' & disease_status == 'ucd' & sample_category == 'tetramer+')

fcs_df <- map_dfr(samples$filename, read_dataset)

gmean <- fcs_df %>% 
  gather(marker, expr, -filename) %>% 
  mutate(expr = sinh(expr)) %>% 
  group_by(marker, filename) %>% 
  summarize(expr = mean(expr)) %>% 
  summarize(expr = mean(expr))

gmean %>% 
  inner_join(markers) %>% 
  mutate(antibody = factor(antibody, levels = rev(markers$antibody))) %>% 
  ggplot(aes('x', antibody)) +
  geom_tile(aes(fill = expr)) +
  scale_fill_gradient(low = 'grey20', high = 'yellow2', limits = c(0, 75), oob = scales::squish) +
  coord_equal()

ggsave('out/plots/fig1d_gut_absexpr.png')
