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

get_tsne <- function(sdf, SEED = 1234, perplexity = 30) {
  print(glue('computing tSNE for {nrow(sdf)} cells'))
  sdf %>% 
    select(-biosource, -unique_name) %>% 
    as.matrix() %>% 
    identity() ->
    M
  set.seed(SEED)
  tsne_out <- Rtsne(M, check_duplicates = FALSE, pca = FALSE, num_threads = 6, perplexity = perplexity)
  
  sdf %>% 
    mutate(tSNE1 = tsne_out$Y[, 1], 
           tSNE2 = tsne_out$Y[, 2])
}
## RECOMPUTE SCS
fcs_df %>% 
  filter(biosource == 'SCS') %>% 
  get_tsne(SEED = 99933, perplexity = 90) %>% 
  identity() ->
  tsne_scs 
tsne_scs %>% 
  inner_join(
    dplyr::select(all_samples, unique_name, sample_category, disease),
    .,
    by = 'unique_name') %>% 
  identity() ->
  tsne_scs

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

# PLOT SCS
save_tsne_plot <- function(plot, bsource, out_folder = 'out/paper_figs/main/tSNE') {
  dir.create(out_folder, showWarnings = F, recursive = T)
  file_formats <- c('svg', 'pdf', 'png') 
  file_formats %>% 
    map_chr(~ glue("{out_folder}/tSNE_{bsource}.{.x}")) %>% 
    walk(~ ggsave(.x, plot))
  bind_rows(
    tsne_pbmc,
    tsne_scs
  ) %>% 
    filter(biosource == bsource) %>% 
    write_csv(glue('{out_folder}/tSNE_{bsource}_raw_data.csv.gz'))
  rect_prct %>% 
    filter(biosource == bsource) %>% 
    write_tsv(glue('{out_folder}/tSNE_{bsource}_rectangle_percentage.tsv'))
}

tsne_scs %>% 
  ggplot(aes(x = tSNE1, y = tSNE2, color = sample_category)) +
  geom_point(data = filter(tsne_scs, sample_category != 'tetramer+'), size = 1) +
  geom_point(data = filter(tsne_scs, sample_category == 'tetramer+'), size = 1) +
  annotate("rect", xmin = tlim$SCS$xmin, xmax = tlim$SCS$xmax, 
           ymin = tlim$SCS$ymin, ymax = tlim$SCS$ymax,
           alpha = .1,
           color = 'black') +
  theme_minimal() +
  scale_color_manual(values = c('tetramer+' = '#f36d16', 'tetramer-' = '#3672bc')) +
  facet_grid(~ disease) +
  coord_fixed() %>% 
  identity() ->
  scs_full

# PLOT PBMC
tsne_pbmc %>% 
  filter(biosource == 'PBMC') %>% 
  mutate(disease_cat = str_c(disease, sample_category, sep = '_')) %>% 
  # filter(sample_category == 'full' & disease == 'Ced') %>% 
  # plot_group = fct_relevel('Ced_tetramer+', 'Ced_full')) %>% 
  ggplot(aes(x = tSNE1, y = tSNE2, color = sample_category)) +
  geom_point(size = 1) +
  annotate("rect", xmin = tlim$PBMC$xmin, xmax = tlim$PBMC$xmax, 
           ymin = tlim$PBMC$ymin, ymax = tlim$PBMC$ymax,
           alpha = .1,
           # fill = NULL,
           color = 'black') +
  # labs(title = 'PBMC - CD1467 - All CD4+',
  #      subtitle = glue('rectangle covers {format(filter(rect_prct, biosource == "PBMC" & sample_category == "full" & disease == "Ced")$prct, digits = 3)}%')) +
  theme_minimal() +
  scale_color_manual(values = c('tetramer+' = '#f36d16', 'full' = '#3672bc')) +
  facet_grid(~ disease_cat) +
  coord_fixed() %>% 
  identity() ->
  pbmc_cd1467_full

save_tsne_plot(scs_full, 'SCS')
save_tsne_plot(pbmc_cd1467_full, 'PBMC')
##
## Supplemental marker intensity plots
##
ac_markers <- read_csv('input_data/original/asbjorn_markers_specific.csv') %>% 
  filter(category == 'CyTOF') %>% 
  dplyr::select(-gene_name, -category) %>% 
  mutate(marker = str_replace_all(common_name, '[- ]', '_'))

tsne_pbmc %>% 
  mutate(tetramer = if_else(sample_category == 'tetramer+', 1.0, 0.0)) %>% 
  # dplyr::select(-donor) %>%
  dplyr::select(-sample_category, -unique_name, -biosource) %>% 
  gather(marker, value, -tSNE1, -tSNE2, -disease) %>% 
  filter(!(marker == 'tetramer' & (value < 1.0 | disease == 'HC'))) %>% 
  dplyr::select(-disease) %>% 
  mutate(marker = fct_relevel(marker, 'tetramer')) %>% 
  filter(marker != 'CD244') %>% 
  group_by(marker) %>% 
  mutate(scaled_intensity = value / max(value)) %>% 
  ungroup() %>% 
  identity() ->
  x_pbmc

x_pbmc %>% 
  left_join(ac_markers, by = 'marker') %>% 
  mutate(marker = if_else(is.na(common_name), marker, common_name)) %>% 
  dplyr::select(-common_name) %>% 
  identity() ->
  x_pbmc

x_pbmc %>% 
  mutate(marker = fct_relevel(marker, ac_markers$common_name)) %>% 
  mutate(marker = fct_relevel(marker, 'tetramer')) %>% 
  filter(!(marker == 'tetramer' & value < 1.0)) %>% 
  # filter(str_detect(marker, 'PD')) %>%
  # head(20000) %>%
  ggplot(aes(tSNE1, tSNE2)) +
  # geom_point(aes(color = scaled_intensity)) +
  # geom_point(aes(color = scaled_intensity), size = 1) +
  stat_summary_2d(aes(z = scaled_intensity), fun = median, bins = 70) +
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
# tsne_pbmc_by_marker
ggsave('out/paper_figs/main/tSNE/tSNE_PBMC_by_marker.png', tsne_pbmc_by_marker)
ggsave('out/paper_figs/main/tSNE/tSNE_PBMC_by_marker.pdf', tsne_pbmc_by_marker)
tsne_scs %>% 
  mutate(tetramer = if_else(sample_category == 'tetramer+', 1.0, 0.0)) %>% 
  dplyr::select(-sample_category, -unique_name, -biosource) %>% 
  gather(marker, value, -tSNE1, -tSNE2, -disease) %>% 
  filter(!(marker == 'tetramer' & (value < 1.0 | disease == 'HC'))) %>% 
  dplyr::select(-disease) %>% 
  mutate(marker = fct_relevel(marker, 'tetramer')) %>% 
  filter(marker != 'CD244') %>% 
  group_by(marker) %>% 
  mutate(scaled_intensity = value / max(value)) %>% 
  ungroup() %>% 
  identity() ->
  x_scs

x_scs %>% 
  left_join(ac_markers, by = 'marker') %>% 
  mutate(marker = if_else(is.na(common_name), marker, common_name)) %>% 
  dplyr::select(-common_name) %>% 
  identity() ->
  x_scs

x_scs %>% 
  mutate(marker = fct_relevel(marker, ac_markers$common_name)) %>% 
  mutate(marker = fct_relevel(marker, 'tetramer')) %>% 
  ggplot(aes(tSNE1, tSNE2)) +
  # geom_point(aes(color = scaled_intensity)) +
  # geom_point(aes(color = scaled_intensity), size = 1) +
  stat_summary_2d(aes(z = scaled_intensity), fun = median, bins = 70) +
  # scale_color_viridis() +
  facet_wrap(~ marker, ncol = 6) +
  theme_minimal() +
  scale_fill_continuous(type = 'viridis', limits = c(0, 1)) +
  labs(fill = 'Scaled Intensity') +
  coord_fixed() %>% 
  identity() ->
  tsne_scs_by_marker
ggsave('out/paper_figs/main/tSNE/tSNE_SCS_by_marker.png', tsne_scs_by_marker)
ggsave('out/paper_figs/main/tSNE/tSNE_SCS_by_marker.pdf', tsne_scs_by_marker)
# x %>% 
#   filter(representative_of == 'PBMC_CeD') %>% 
#   filter(marker %in% c('tetramer', 'CXCR3')) %>%
#   ggplot(aes(tSNE1, tSNE2)) +
#   geom_point(aes(color = scaled_intensity)) +
#   facet_wrap(~ marker, ncol = 5)
p_pbmc_marker <- x %>% 
  filter(representative_of == 'PBMC_CeD') %>% 
  # filter(marker == 'CXCR3') %>% 
  ggplot(aes(tSNE1, tSNE2)) +
  # geom_point(aes(color = scaled_intensity)) +
  geom_point(aes(color = scaled_intensity, alpha = scaled_intensity), size = 1) +
  # stat_summary_2d(aes(z = value), fun = median, bins = 70) +
  scale_color_viridis() +
  facet_wrap(~ marker, ncol = 5)
p_scs_marker <- x %>% 
  filter(representative_of == 'SCS_CeD') %>% 
  # filter(marker == 'CXCR3') %>% 
  ggplot(aes(tSNE1, tSNE2)) +
  geom_point(aes(color = scaled_intensity, alpha = scaled_intensity), size = 1) +
  # stat_summary_2d(aes(z = value), fun = median, bins = 70) +
  scale_color_viridis() +
  facet_wrap(~ marker, ncol = 5)
dir.create('out/paper_figs/supplemental', showWarnings = FALSE, recursive = TRUE)
ggsave('out/paper_figs/supplemental/tSNE_CeD_PBMC_by_marker.svg', p_pbmc_marker)
ggsave('out/paper_figs/supplemental/tSNE_CeD_PBMC_by_marker.png', p_pbmc_marker)
ggsave('out/paper_figs/supplemental/tSNE_CeD_PBMC_by_marker.pdf', p_pbmc_marker)
ggsave('out/paper_figs/supplemental/tSNE_CeD_SCS_by_marker.svg', p_scs_marker)
ggsave('out/paper_figs/supplemental/tSNE_CeD_SCS_by_marker.png', p_scs_marker)
ggsave('out/paper_figs/supplemental/tSNE_CeD_SCS_by_marker.pdf', p_scs_marker)