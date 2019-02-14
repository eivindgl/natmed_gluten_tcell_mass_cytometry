# t-SNE plotting fun
pacman::p_load(
  tidyverse,
  rhdf5,
 glue,
 tictoc,
 assertthat,
 svglite
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

pacman::p_load(
  devtools,
  viridis
)
devtools::install_github("jkrijthe/Rtsne", ref = 'openmp')
library(Rtsne)

sne_samples %>% 
  group_by(biosource) %>% 
  do(.$unique_name %>% 
       map_dfr(read_dataset)
  ) %>% 
  ungroup() %>% 
  identity() ->
  fcs_df

get_tsne <- function(sdf) {
  print(glue('computing tSNE for {nrow(sdf)} cells'))
  sdf %>% 
    select(-biosource, -unique_name) %>% 
    as.matrix() %>% 
    identity() ->
    M
  set.seed(12345678)
  tsne_out <- Rtsne(M, check_duplicates = FALSE, pca = FALSE, num_threads = 6)
  
  sdf %>% 
    mutate(tSNE1 = tsne_out$Y[, 1], 
           tSNE2 = tsne_out$Y[, 2])
}
fcs_df %>% 
  group_by(biosource) %>% 
  do(get_tsne(.)) %>% 
  ungroup() %>% 
  identity() ->
  tsne_df

sne_samples %>% 
  select(unique_name, sample_category, donor, disease) %>% 
  inner_join(tsne_df) %>% 
  group_by(biosample) %>% 
  do(
    plot = ggplot(., aes(x = tSNE1, y = tSNE2, color = sample_category)) +
      geom_point() +
      facet_wrap(~ interaction(disease, sample_category), nrow = 1)
  ) %>% 
  identity() ->
  sne_plots

##
## MAIN TETRAMER PLOTS
##
p_pbmc_main <- sne_plots %>% 
  filter(representative_of == 'PBMC_CeD') %>% 
  .$plot %>% 
  magrittr::extract2(1) + 
  theme_minimal() +
  scale_color_manual(values = c('tetramer+' = '#f36d16', 'tetramer-' = '#3672bc'))
p_pbmc_main
p_scs_main <- sne_plots %>% 
  filter(representative_of == 'SCS_CeD') %>% 
  .$plot %>% 
  magrittr::extract2(1) + 
  theme_minimal() +
  scale_color_manual(values = c('tetramer+' = '#f36d16', 'tetramer-' = '#3672bc'))
p_scs_main
ggsave('out/paper_figs/main/tSNE_CeD_PBMC.svg', p_pbmc_main)
ggsave('out/paper_figs/main/tSNE_CeD_SCS.svg', p_scs_main)

##
## Plots with rectangle annotation
##
tetplim <- tribble(
  ~biosource, ~xmin, ~xmax, ~ymin, ~ymax,
  "PBMC", -10.2, -7, -18.5, -15.2,
  "SCS", 1, 15, -30, -13
)
tlim <- split(tetplim, tetplim$biosource)

tsne_df %>% 
  inner_join(select(all_samples, unique_name, sample_category, biosource, disease)) %>% 
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
  
sne_samples %>% 
  select(unique_name, sample_category, donor, disease) %>% 
  inner_join(tsne_df) %>% 
  mutate(plot_group = str_c(disease, sample_category, sep = '_')) %>% # ,
  mutate(plot_group = fct_relevel(plot_group, 'Ced_tetramer+', 'Ced_full')) %>% 
  identity() ->
  xdf
xdf %>% 
  filter(biosource == 'PBMC') %>% 
  filter(sample_category == 'tetramer+') %>% 
  ggplot(aes(x = tSNE1, y = tSNE2, color = sample_category)) +
  geom_point() +
  annotate("rect", xmin = tlim$PBMC$xmin, xmax = tlim$PBMC$xmax, 
             ymin = tlim$PBMC$ymin, ymax = tlim$PBMC$ymax,
             alpha = .2,
           color = 'black') +
  labs(title = 'PBMC - CD1467 - tetramer+',
       subtitle = glue('rectangle covers {format(filter(rect_prct, biosource == "PBMC" & sample_category == "tetramer+")$prct, digits = 3)}%')) +
  theme_minimal() +
  scale_color_manual(values = c('tetramer+' = '#f36d16', 'tetramer-' = '#3672bc')) %>% 
  identity() ->
  pbmc_cd1467_tp
pbmc_cd1467_tp
ggsave('out/paper_figs/main/tSNE_CeD_PBMC_CD1467_tetp.svg', pbmc_cd1467_tp)
xdf %>% 
  filter(biosource == 'PBMC') %>% 
  # filter(sample_category == 'full' & disease == 'Ced') %>% 
         # plot_group = fct_relevel('Ced_tetramer+', 'Ced_full')) %>% 
  ggplot(aes(x = tSNE1, y = tSNE2, color = sample_category)) +
  geom_point() +
  annotate("rect", xmin = tlim$PBMC$xmin, xmax = tlim$PBMC$xmax, 
             ymin = tlim$PBMC$ymin, ymax = tlim$PBMC$ymax,
             alpha = .3,
           color = 'black') +
  # labs(title = 'PBMC - CD1467 - All CD4+',
  #      subtitle = glue('rectangle covers {format(filter(rect_prct, biosource == "PBMC" & sample_category == "full" & disease == "Ced")$prct, digits = 3)}%')) +
  theme_minimal() +
  scale_color_manual(values = c('tetramer+' = '#f36d16', 'full' = '#3672bc')) +
  facet_grid(~ plot_group) +
  coord_fixed() %>% 
  identity() ->
  pbmc_cd1467_full
pbmc_cd1467_full
ggsave('out/paper_figs/main/tSNE_CeD_PBMC_CD1467_BC13.svg', pbmc_cd1467_full)
ggsave('out/paper_figs/main/tSNE_CeD_PBMC_CD1467_BC13.pdf', pbmc_cd1467_full)
ggsave('out/paper_figs/main/tSNE_CeD_PBMC_CD1467_BC13.png', pbmc_cd1467_full)
x %>% 
  filter(representative_of == 'PBMC_CeD') %>% 
  write_csv('out/paper_figs/main/tSNE_CeD_PBMC_CD1467_BC13_raw_data.csv.gz')
rect_prct %>% 
  filter(biosource == 'PBMC') %>% 
  write_csv('out/paper_figs/main/tSNE_CeD_PBMC_CD1467_BC13_rectange_prct.csv')

xdf %>% 
  filter(biosource == 'SCS') %>% 
  filter(!(disease == 'HC' & sample_category == 'tetramer+')) %>% 
  # filter(sample_category == 'full' & disease == 'Ced') %>% 
         # plot_group = fct_relevel('Ced_tetramer+', 'Ced_full')) %>% 
  ggplot(aes(x = tSNE1, y = tSNE2, color = sample_category)) +
  geom_point() +
  annotate("rect", xmin = tlim$SCS$xmin, xmax = tlim$SCS$xmax, 
             ymin = tlim$SCS$ymin, ymax = tlim$SCS$ymax,
             alpha = .1,
           color = 'black') +
  theme_minimal() +
  scale_color_manual(values = c('tetramer+' = '#f36d16', 'tetramer-' = '#3672bc')) +
  facet_grid(~ plot_group) +
  coord_fixed() %>% 
  identity() ->
  scs_cd1467_full
scs_cd1467_full

xdf %>% 
  write_csv('out/tmp/tsne.csv.gz')

x %>% 
  filter(representative_of == 'PBMC_CeD') %>% 
  filter(sample_category == 'tetramer-') %>% 
  ggplot(aes(x = tSNE1, y = tSNE2, color = sample_category)) +
  geom_point() +
  annotate("rect", xmin = tlim$PBMC$xmin, xmax = tlim$PBMC$xmax, 
             ymin = tlim$PBMC$ymin, ymax = tlim$PBMC$ymax,
             alpha = .2) +
  labs(title = 'PBMC - CD1467 - tetramer+',
       subtitle = glue('rectangle covers {format(filter(rect_prct, biosource == "PBMC" & sample_category == "tetramer+")$prct, digits = 4)}%')) +
  theme_minimal() +
  scale_color_manual(values = c('tetramer+' = '#f36d16', 'tetramer-' = '#3672bc')) %>% 
  identity() ->
  pbmc_cd1467_tp
ggsave('out/paper_figs/main/tSNE_CeD_PBMC_CD1467_tetp.svg', pbmc_cd1467_tp)
##
## Supplemental marker intensity plots
##
tsne_df %>% 
  inner_join(
    select(all_samples, unique_name, sample_category)
  ) %>% 
  mutate(tetramer = if_else(sample_category == 'tetramer+', 1.0, 0.0)) %>% 
  dplyr::select(-sample_category) %>% 
  gather(marker, value, -representative_of, -unique_name, -tSNE1, -tSNE2) %>% 
  mutate(marker = fct_relevel(marker, 'tetramer')) %>% 
  inner_join(
    select(all_samples, unique_name, disease, disease_status)
  ) %>% 
  group_by(marker) %>% 
  mutate(scaled_intensity = value / max(value)) %>% 
  ungroup() %>% 
  identity() ->
  x

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


