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

get_tsne <- function(sdf, SEED = 1234, perplexity = 30) {
  print(glue('computing tSNE for {nrow(sdf)} cells'))
  sdf %>% 
    select(-unique_name) %>% 
    as.matrix() %>% 
    identity() ->
    M
  set.seed(SEED)
  tsne_out <- Rtsne(M, check_duplicates = FALSE, pca = FALSE, num_threads = 6, perplexity = perplexity)
  
  sdf %>% 
    mutate(tSNE1 = tsne_out$Y[, 1], 
           tSNE2 = tsne_out$Y[, 2])
}

ac_markers <- read_csv('input_data/original/asbjorn_markers_specific.csv') %>% 
  filter(category == 'CyTOF') %>% 
  dplyr::select(-gene_name, -category) %>% 
  mutate(marker = str_replace_all(common_name, '[- ]', '_'))

all_samples <- read_csv('input_data/renamed/sample_meta_full.csv') %>% 
  filter(biosource == 'PBMC')

list(
  ucd_tetp = all_samples %>% 
    filter(disease_status == 'UCD' & sample_category == 'tetramer+'),
  gfd_tetp = all_samples %>% 
    filter(disease_status == 'GFD' & sample_category == 'tetramer+'),
  autoimmune = all_samples %>% 
    filter(sample_category == 'AutoPhenotype') %>% 
    filter(sample_group == 'PBMC_AutoImmune') %>% 
    mutate(disease = case_when(str_detect(disease, '^SSc') ~ 'SSc',
                               disease == 'SLE' ~ 'SLE',
                               TRUE ~ 'Other Autoimmune')),
  flu = all_samples %>% 
    filter(sample_category == 'AutoPhenotype') %>% 
    filter(disease == 'Flu'),
  HC = all_samples %>% 
    filter(unique_name == 'BC11_full_N2')
) %>% 
  bind_rows(.id = 'tsneCat') %>% 
  identity() ->
  sne_samples

sne_samples %>% 
  .$unique_name %>% 
  map_dfr(read_dataset) %>% 
  select(unique_name, ac_markers$marker) %>% 
  select_if(function(x) !any(is.na(x))) %>% 
  select(-Integrin_b7, -CD49d) %>% 
  identity() ->
  fcs_df
  
## RECOMPUTE SCS
fcs_df %>% 
  get_tsne(SEED = 1234, perplexity = 30) %>% 
  identity() ->
  tsne_df 

tsne_df %>% 
  select(-starts_with('sample_category')) %>%
  select(-starts_with('disease')) %>%
  inner_join(
    dplyr::select(sne_samples, unique_name, sample_category, disease, disease_status),
    .,
    by = 'unique_name') %>% 
  mutate(disease = case_when(disease_status == 'UCD' ~ 'CeD_UCD',
                             disease_status == 'GFD' ~ 'CeD_TCD',
                             TRUE ~ disease)) %>% 
  mutate(disease= fct_relevel(disease, 'HC', 'CeD_UCD', 'CeD_TCD', 'SLE', 'SSc')) %>% 
  identity() ->
  tsne_df

if (!file.exists('out/paper_figs/supplemental/tsne/tSNE_tetp_and_autophenotype_against_BC11_raw_data.csv.gz')) {
  tsne_df %>% 
    write_csv('out/paper_figs/supplemental/tsne/tSNE_tetp_and_autophenotype_against_BC11_raw_data.csv.gz')
} else {
  print("file exists. skipping auto overwrite of raw data")
}
# tsne_df <- read_csv('out/paper_figs/supplemental/tsne/tSNE_tetp_and_autophenotype_against_BC11_raw_data.csv.gz') %>%
#   mutate(disease= fct_relevel(disease, 'HC', 'CeD_UCD', 'CeD_TCD', 'SLE', 'SSc'))

tsne_df %>%  # this donor had cancer as well as flu (and tons of AT+ cells)
  filter(!str_detect(unique_name, '22-002-01-12')) %>% 
  identity() ->
  tsne_df

tlim = list(xmin = 19.5, xmax = 26,
            ymin = -10, ymax = -3.5)


tsne_df %>% 
  filter(disease != 'CeD_TCD') %>% 
  identity() ->
  xdf
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette <- setNames(cbPalette, levels(xdf$disease))
p_tsne <- xdf %>%
  ggplot(aes(x = tSNE1, y = tSNE2, color = disease)) +
  # geom_point(data = head(filter(xdf, disease == 'HC'), 30000), size = 1, alpha = 1) +
  geom_point(data = filter(xdf, disease == 'HC'), size = 1, alpha = 0.1) +
  geom_point(data = filter(xdf, disease != 'HC'), size = 1) +
  annotate("rect", xmin = tlim$xmin, xmax = tlim$xmax,
             ymin = tlim$ymin, ymax = tlim$ymax,
             alpha = .1,
           color = 'black') +
  theme_minimal() +
  # scale_color_manual(values = c('tetramer+' = '#f36d16', 'tetramer-' = '#3672bc')) +
  # facet_grid(~ sample_category) +
  coord_fixed() +
  scale_colour_manual(values=cbPalette)
# PLOT SCS
p_tsne
dir.create('out/paper_figs/supplemental/tsne', recursive = T, showWarnings = F)
ggsave('out/paper_figs/supplemental/tsne/tSNE_tetp_and_autophenotype_against_BC11.pdf', p_tsne)
ggsave('out/paper_figs/supplemental/tsne/tSNE_tetp_and_autophenotype_against_BC11.png', p_tsne)

xdf %>% 
  crossing(as_data_frame(tlim)) %>% 
  mutate(within_rect = tSNE1 > xmin & tSNE1 < xmax &
           tSNE2 > ymin & tSNE2 < ymax) %>% 
  group_by(sample_category, disease) %>% 
  summarize(xmin = first(xmin), xmax = first(xmax), ymin = first(ymin), ymax = first(ymax),
    tot = n(), within = sum(within_rect), prct = within / tot * 100.0) %>% 
  ungroup() %>% 
  identity() ->
  rect_prct
rect_prct
rect_prct %>% 
  write_csv('out/paper_figs/supplemental/tsne/tSNE_tetp_and_autophenotype_against_BC11.csv')
