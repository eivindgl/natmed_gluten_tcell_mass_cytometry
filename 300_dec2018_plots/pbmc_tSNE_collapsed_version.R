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

predmarkers <- read_csv('input_data/common_markers_for_pred.csv') %>%
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
    mutate(unique_name = unique_name) %>% 
    select(unique_name, everything())
}

all_samples <- read_csv('input_data/renamed/sample_meta_full.csv')
tibble(unique_name = c("UCD1467P.A_tetramer+", "UCD1467P.A_full", "UCD1467A.A_tetramer+", 
                       "UCD1467A.A_tetramer-", "CD1413D.A_tetramer+", "CD1413D.A_tetramer-",
                       'BC11_full_N2')) %>% 
  inner_join(all_samples) %>% 
  identity() ->
  sne_samples


sne_samples %>% 
  group_by(biosource) %>% 
  do(.$unique_name %>% 
       map_dfr(read_dataset)
  ) %>% 
  ungroup() %>% 
  identity() ->
  fcs_df

tsne_pbmc <- read_csv('out/paper_figs/main/tSNE/tSNE_PBMC_raw_data.csv.gz')
tsne_scs <- read_csv('out/paper_figs/main/tSNE/tSNE_SCS_raw_data.csv.gz')

##
## Plots with rectangle annotation
##
tetplim <- tribble(
  ~biosource, ~xmin, ~xmax, ~ymin, ~ymax,
  "PBMC", 16.4, 18.4, 5, 10,
  "SCS", -7, -2, 11.5, 21.8
  # "SCS", -7, -4, 15, 21.8 # 62% vs 2.5%
)
tlim <- split(tetplim, tetplim$biosource)

bind_rows(
  tsne_scs,
  tsne_pbmc
) %>% 
  inner_join(tetplim) %>% 
  mutate(within_rect = tSNE1 > xmin & tSNE1 < xmax &
           tSNE2 > ymin & tSNE2 < ymax) %>% 
  group_by(biosource, sample_category, disease) %>% 
  summarize(xmin = first(xmin), xmax = first(xmax), ymin = first(ymin), ymax = first(ymax),
            tot = n(), within = sum(within_rect), prct = within / tot * 100.0) %>% 
  ungroup() %>% 
  identity() ->
  rect_prct
rect_prct


z <- tsne_pbmc %>% 
  filter(biosource == 'PBMC') %>% 
  mutate(disease_cat = str_c(disease, sample_category, sep = '_'))
z %>% 
  ggplot(aes(x = tSNE1, y = tSNE2, color = sample_category)) +
  geom_point(data = filter(z, sample_category == 'full'), size = 1) +
  geom_point(data = filter(z, sample_category == 'tetramer+'), size = 1) +
  annotate("rect", xmin = tlim$PBMC$xmin, xmax = tlim$PBMC$xmax,
           ymin = tlim$PBMC$ymin, ymax = tlim$PBMC$ymax,
           alpha = .1,
           # fill = NULL,
           color = 'black') +
  # labs(title = 'PBMC - CD1467 - All CD4+',
  #      subtitle = glue('rectangle covers {format(filter(rect_prct, biosource == "PBMC" & sample_category == "full" & disease == "Ced")$prct, digits = 3)}%')) +
  # geom_density2d(data = z, aes(color = sample_category, group = sample_category)) + 
  theme_minimal() +
  scale_color_manual(values = c('tetramer+' = '#f36d16', 'full' = '#3672bc')) +
  facet_grid(disease ~ .) +
  theme(axis.text = element_blank(),
        legend.position = NaN,
        strip.text = element_blank()) +
  labs(x = '', y = '') +
  coord_fixed() %>% 
  identity() ->
  pbmc_cd1467_full
pbmc_cd1467_full
ggsave('out/paper_figs/tSNE_PBMC_collapsed_version.png', plot = pbmc_cd1467_full, bg = 'transparent')
ggsave('out/paper_figs/tSNE_PBMC_collapsed_version.pdf', plot = pbmc_cd1467_full, bg = 'transparent')
