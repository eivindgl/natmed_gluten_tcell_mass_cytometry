pacman::p_load(
  tidyverse,
  rhdf5,
  glue,
  tictoc,
  ggrepel,
  viridis
)
dir.create('out/ai_vs_flu', showWarnings = F)

# predmarkers <- 
tibble(antibody = c("CCR4", "CCR6", "CD127", "CD161", "CD25", "CD27", "CD28", "CD38", 
                    "CD39", "CD57", "CD69", "CD73", "CTLA-4", "CXCR3", "CXCR5", "HLA-DR", 
                    "ICOS", "KLRG1", "PD-1", "CD45RA", "CD62L", "CCR7")) %>% 
  mutate(fcs_name = str_replace_all(antibody, '[- ]', '_')) %>% 
  filter(!(antibody %in% c('CD73'))) %>% 
  identity() ->
  predmarkers

read_dataset <- function(path, h5_path) {
  M <- h5read(path, h5_path)
  colnames(M) <- h5readAttributes(path, h5_path)$colnames
  M
}

all_samples  <- read_csv('input_data/renamed/sample_meta_full.csv')

all_samples %>% 
  count(sample_group) # PBMC_AutoImmune | PBMC_Flu | PBMC_UCD

all_samples %>% 
  filter(sample_group %in% c('PBMC_AutoImmune', 'PBMC_Flu', 'PBMC_UCD') | disease_status == 'Challenge') %>% 
  filter(donor != 'CCD1535') %>% # no pre enriched
  filter(sample_category != 'tetramer-') %>% 
  filter(sample_category != 'AutoPhenotype') %>% 
  mutate(sample_group = if_else(disease_status == 'Challenge', 'PBMC_CCD', sample_group)) %>% 
  # filter(disease_status != 'recovering') %>% 
  dplyr::select(-category, -instrument, -date, -note) %>% 
  mutate(disease = case_when(sample_group == 'PBMC_AutoImmune' ~
                               case_when(str_detect(disease, '^SSc') ~ 'SSc',
                                         disease == 'SLE' ~ 'SLE',
                                         TRUE ~ 'Other Autoimmune'), 
                             disease == 'Flu' ~ str_c(disease, disease_status, sep = '_'),
                             TRUE ~ disease)) %>% 
  filter(disease != 'Other Autoimmune') %>% 
  filter(donor != '22-002') %>% # Patient has comorbidities... exclude
  identity() ->
  ai_flu_samples

alldf <- ai_flu_samples %>% 
  mutate(h5path = glue('samples/{unique_name}')) %>%
  .$h5path %>%
  map(partial(read_dataset, 'out/cytof.h5')) %>%
  map(~ as_tibble(.x[, predmarkers$fcs_name])) %>% 
  set_names(ai_flu_samples$unique_name) %>% 
  bind_rows(.id = 'unique_name') 

ai_flu_samples %>% 
  group_by(sample_category) %>% 
  sample_n(7) %>% 
  ungroup() %>% 
  inner_join(alldf %>% 
               dplyr::select(unique_name, CD38, HLA_DR) %>% 
               gather(marker, expression, -unique_name)) %>% 
  ggplot(aes(sample_category, expression)) +
  geom_boxplot(aes(color = marker))

ai_flu_samples %>% 
  filter(sample_category == 'tetramer+') %>% 
  semi_join(alldf, .) %>% 
  identity() ->
  gdf
ai_flu_samples %>% 
  filter(sample_category == 'tetramer+') %>% 
  anti_join(alldf, .) %>% 
  identity() ->
  rdf

## Based on these initial results I define activated T-cells
## to be HLA-DR>0 and CD38 > 1.5 (20th percentile for gsT-cells)
##
quantile(gdf$CD38, probs = 0.2)    
quantile(rdf$CD38, probs = 0.65)  

quantile(gdf$HLA_DR, probs = 0.3)    
quantile(rdf$HLA_DR, probs = 0.7)    

gdf %>% 
  mutate(has_HLA_DR = HLA_DR > 0,
         strong_CD38 = CD38 > 1.5,
         mature = CD45RA < 3,
         activated = has_HLA_DR & strong_CD38) %>% 
  dplyr::count(has_HLA_DR, strong_CD38, mature)
rdf %>% 
  mutate(has_HLA_DR = HLA_DR > 0,
         strong_CD38 = CD38 > 1.5,
         mature = CD45RA < 1,
         activated = has_HLA_DR & strong_CD38) %>% 
dplyr::count(has_HLA_DR, strong_CD38)

alldf %>% 
  inner_join(ai_flu_samples) %>% 
  filter(sample_category == 'full') %>% 
  mutate(activated = if_else(CD38 > 1.5 & CD45RA < 3, 'activated', 'naive_or_dormat')) %>% # & HLA_DR > 0) %>% 
  dplyr::count(disease, activated) %>% 
  spread(activated, n) %>% 
  mutate(total = activated + naive_or_dormat,
         frac = activated / total) %>% 
  arrange(frac)
  

alldf %>% 
  # bind_rows(gdf) %>% 
  inner_join(dplyr::select(ai_flu_samples, unique_name, sample_group, sample_category), .) %>% 
  filter(CD38 > 1.5 & CD45RA < 3) %>% # & HLA_DR > 0) %>% 
  add_count(unique_name) %>% 
  filter(sample_category != 'tetramer+' & n < 800) %>% 
  distinct(unique_name)
alldf %>% 
  # bind_rows(gdf) %>% 
  inner_join(dplyr::select(ai_flu_samples, unique_name, sample_group, sample_category), .) %>% 
  filter(CD38 > 1.5 & CD45RA < 3) %>% # & HLA_DR > 0) %>% 
  add_count(unique_name) %>% 
  filter(sample_category == 'tetramer+' | n > 800) %>% 
  dplyr::select(-n) %>% 
  mutate(cell_n = 1:n()) %>% 
  select(cell_n, everything()) %>% 
  select(-CD45RA) %>% 
  identity() ->
  adf

adf %>% 
  distinct(unique_name) %>% 
  inner_join(ai_flu_samples) %>% 
  distinct(donor, disease, disease_status) %>% 
  dplyr::count(disease, disease_status)

adf %>% 
  distinct(unique_name) %>% 
  inner_join(ai_flu_samples) %>% 
  dplyr::count(donor, disease, disease_status) %>% 
  filter(disease_status == 'Challenge')

adf %>%
  write_csv('out/ai_vs_flu/activated_tcells.csv.gz')
adf %>% 
  inner_join(ai_flu_samples) %>% 
  dplyr::count(unique_name, donor, disease) %>% 
  arrange(desc(n)) %>% print(n = Inf)
## not correct as long as I remove samples with low counts
# print(glue("{nrow(adf)} out of {nrow(rdf)} are activated. ({round((nrow(adf) / nrow(rdf)) * 100, digits = 1)}%)"))

#
## compute TSNE_plot
#
devtools::install_github("jkrijthe/Rtsne", ref = 'openmp')
library(Rtsne)

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
    mutate(tSNE1 = tsne_out$Y[, 1], 
           tSNE2 = tsne_out$Y[, 2])
}
min_group <- adf %>% 
   filter(sample_category == 'full') %>% 
  count(unique_name) %>% 
  .$n %>% 
  min()

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
write_csv(tsne_tweak_df, 'out/ai_vs_flu/tsne_tweak_df.csv')

tsne_tweak_df %>% 
  inner_join(ai_flu_samples) %>% 
  mutate(disease = if_else(disease == 'Ced', str_c('CeD_', disease_status), disease),
         disease = if_else(sample_category == 'tetramer+', str_c(disease, '_tet+'), disease)) %>% 
  # mutate(sample_group = if_else(sample_category == 'tetramer+', 'PBMC_UCD_tetp', sample_group)) %>% 
  # mutate(disease = if_else(sample_category == 'tetramer+', 'CeD_tet+', disease)) %>% 
  ggplot(aes(tSNE1, tSNE2)) +
  geom_point(alpha = 0.5) +
  geom_density_2d() +
  facet_wrap(~ disease, ncol = 2) +
  theme_minimal()
ggsave('out/ai_vs_flu/tSNE_by_group.png', height = 12)

plot_df <-tsne_tweak_df %>% 
  inner_join(ai_flu_samples) %>% 
  filter(sample_category == 'full') %>% 
  group_by(disease) %>%
  do(
    plots = ggplot(data = .) + aes(tSNE1, tSNE2) +
      geom_point(alpha = 0.5) +
      geom_density_2d() + facet_wrap(~ unique_name) + ggtitle(.$disease)
  ) 
plot_df %>%
  rowwise() %>%
  do(ggsave(filename = glue('out/ai_vs_flu/tSNE_by_donor_PBMC_{.$disease}.png'), plot = .$plots ))
     
tsne_markers <- keep(predmarkers$fcs_name, function(x) x != 'CD45RA')

tsne_tweak_df %>% 
  gather(marker, expression, -unique_name, -sample_group, -tSNE1, -tSNE2, -cell_n, -sample_category) %>% 
  group_by(marker) %>% 
  mutate(frac_of_max = expression / max(expression)) %>% 
  ungroup() %>% 
  mutate(marker = fct_relevel(marker, tsne_markers)) %>% 
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
ggsave('out/ai_vs_flu/tsne_by_marker.png', tsne_pbmc_by_marker, width = 10, height = 10)
ggsave('out/ai_vs_flu/tsne_by_marker.pdf', tsne_pbmc_by_marker, width = 10, height = 10)

