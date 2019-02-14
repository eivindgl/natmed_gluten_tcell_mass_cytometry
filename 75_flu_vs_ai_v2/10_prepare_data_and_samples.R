#
# 1. Selects samples and writes out/samples.csv
# 2. Extracts activated T-cells and cytof values for selected markers to out/activated_tcells.csv
#
pacman::p_load(
  tidyverse,
  rhdf5,
  glue,
  tictoc,
  ggrepel,
  viridis,
  assertthat
)
SEED <- 12345
outdir <- '75_flu_vs_ai_v2/out'
dir.create(outdir, showWarnings = F)

# predmarkers <- 
tibble(antibody = c("CCR4", "CCR6", "CD127", "CD161", "CD25", "CD27", "CD28", "CD38", 
                    "CD39", "CD57", "CD69", "CD73", "CTLA-4", "CXCR3", "CXCR5", "HLA-DR", 
                    "ICOS", "KLRG1", "PD-1", "CD45RA", "CD62L", "CCR7")) %>% 
  mutate(fcs_name = str_replace_all(antibody, '[- ]', '_')) %>% 
  filter(!(antibody %in% c('CD73'))) %>% 
  identity() ->
  predmarkers
predmarkers %>% 
  write_csv(file.path(outdir, 'predmarkers.csv'))

read_dataset <- function(path, h5_path) {
  M <- h5read(path, h5_path)
  colnames(M) <- h5readAttributes(path, h5_path)$colnames
  M
}

all_samples  <- read_csv('input_data/renamed/sample_meta_full.csv')

all_samples %>% 
  filter(sample_group %in% c('PBMC_AutoImmune', 'PBMC_Flu', 'PBMC_UCD', 'PBMC_HC', 'PBMC_GFD') | disease_status == 'Challenge') %>% 
  filter(donor != 'CD1535') %>% # no pre enriched
  filter(donor != 'CD1570') %>% # CeD no tet+
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
  filter(sample_category != 'tetramer-') %>% 
  identity() ->
  pot_samples

challenge_donors <- pot_samples %>% 
  filter(sample_group == 'PBMC_CCD') %>% 
  distinct(donor) %>% 
  .$donor

pot_samples %>% 
  mutate(sample_group = if_else(disease_status == 'GFD' & donor %in% challenge_donors, 'PBMC_GFD', sample_group)) %>% 
  identity() ->
  pot_samples

pot_samples %>% 
  dplyr::select(-sample, -biosource) %>% 
  mutate(disease = if_else(sample_group == 'PBMC_Flu', 'Flu', disease),
         disease = tolower(disease),
         disease_status = tolower(disease_status),
         sample_group = case_when(sample_group == 'PBMC_Flu' ~ str_c(sample_group, str_sub(disease_status, 1 ,3), sep = '_'), 
                                  sample_group == 'PBMC_AutoImmune' ~ str_c('PBMC', disease, sep = '_'),
                                  TRUE ~ sample_group),
         sample_group = tolower(str_sub(sample_group, 6))) %>% 
  identity() ->
  pot_samples

## READ all samples to ensure each selected sample has enough cells
alldf <- pot_samples %>% 
  mutate(h5path = glue('samples/{unique_name}')) %>%
  .$h5path %>%
  map(partial(read_dataset, 'out/cytof.h5')) %>%
  map(~ as_tibble(.x[, predmarkers$fcs_name])) %>% 
  set_names(pot_samples$unique_name) %>% 
  bind_rows(.id = 'unique_name') 

alldf %>% 
  # bind_rows(gdf) %>% 
  inner_join(dplyr::select(pot_samples, unique_name, sample_group, sample_category), .) %>% 
  filter(CD38 > 1.5 & CD45RA < 3) %>% # & HLA_DR > 0) %>% 
  add_count(unique_name) %>% 
  filter(sample_category == 'tetramer+' | n > 800) %>% 
  dplyr::select(-n) %>% 
  mutate(cell_n = 1:n()) %>% 
  select(cell_n, everything()) %>% 
  select(-CD45RA) %>% 
  identity() ->
  a_raw_df


a_raw_df %>% 
  dplyr::count(unique_name) %>% 
  dplyr::rename(ncells = n) %>% 
  inner_join(pot_samples) %>% 
  identity() ->
  pot_samples

pot_samples %>%
  filter(sample_category == 'full' & ncells > 807) %>% 
  dplyr::count(sample_group)
  

##
## Subset to 4-5 samples per sample group
##

## Use same donors for Challenge and GFD
ccd_samples <-  pot_samples %>% 
  filter(sample_group == 'ccd')
pot_samples %>% 
  filter(sample_group == 'ccd') %>% 
  distinct(donor) %>% 
  inner_join(pot_samples) %>% 
  filter(sample_group == 'gfd') %>% 
  identity() ->
  gfd_samples
assert_that(nrow(distinct(gfd_samples, donor)) == 4)
assert_that(nrow(gfd_samples) <= 8)
set.seed(SEED) # add 1 extra to get 5 samples
gfd_samples <- pot_samples %>% 
  filter(sample_group == 'gfd') %>% 
  anti_join(gfd_samples) %>% 
  distinct(donor) %>% 
  sample_n(1) %>% 
  inner_join(pot_samples) %>% 
  bind_rows(gfd_samples)
assert_that(nrow(distinct(gfd_samples, donor)) == 5)
assert_that(nrow(distinct(filter(gfd_samples, ncells >= 807), donor)) == 5)

# ## Add matching flu samples
# ##
# set.seed(SEED)
# flu_samples <- samples %>% 
#   filter(sample_group == 'flu_inf') %>% 
#   distinct(donor) %>% 
#   inner_join(samples) %>% 
#   select(donor, sample_group, ncells) %>% 
#   spread(sample_group, ncells) %>% 
#   na.omit() %>% 
#   distinct(donor) %>% 
#   sample_n(5) %>% 
#   inner_join(samples)

## sample HC controls from the bloodbank and not the gastrolab
set.seed(SEED)
pot_samples %>% 
  filter(sample_group != 'hc' | str_detect(donor, '^BC')) %>% 
  filter(sample_group != 'gfd' & sample_group != 'ccd') %>% 
  filter(sample_category == 'full') %>% 
  filter(ncells >= 807) %>% 
  group_by(sample_group) %>% 
  sample_n(5) %>% 
  ungroup() %>% 
  distinct(donor) %>% 
  inner_join(pot_samples) %>% # include tetramer+ samples
  # dplyr::count(sample_group, sample_category) %>% 
  bind_rows(ccd_samples) %>% 
  bind_rows(gfd_samples) %>% 
  identity() ->
  samples
samples %>% 
  dplyr::count(sample_group, sample_category)
  
samples %>% 
  dplyr::select(-path) %>% 
  write_csv(file.path(outdir, 'samples.csv'))

samples %>% 
  filter(sample_category == 'tetramer+') %>% 
  semi_join(alldf, .) %>% 
  identity() ->
  gdf
samples %>% 
  filter(sample_category == 'tetramer+') %>% 
  anti_join(alldf, .) %>% 
  identity() ->
  rdf

## Based on these initial results I define activated T-cells
## to be HLA-DR>0 and CD38 > 1.5 (20th percentile for gsT-cells)
##

rdf %>% 
  inner_join(samples) %>% 
  filter(sample_category == 'full') %>% 
  mutate(activated = if_else(CD38 > 1.5 & CD45RA < 3, 'activated', 'naive_or_dormat')) %>% # & HLA_DR > 0) %>% 
  dplyr::count(sample_group, donor, activated) %>% 
  spread(activated, n) %>% 
  mutate(total = activated + naive_or_dormat,
         frac = activated / total) %>% 
  arrange(frac) %>% 
  identity() ->
  x
x %>% 
  write_csv(file.path(outdir, 'fraction_activated_per_donor_and_state.csv'))
x %>% 
  dplyr::select(-total, -frac) %>% 
  gather(state, n, -sample_group, -donor) %>% 
  dplyr::count(sample_group, state, wt = n) %>% 
  spread(state, nn) %>% 
  mutate(total = activated + naive_or_dormat,
         frac = activated / total) %>% 
  arrange(frac) %>% 
  identity() ->
  xx
xx %>% 
  write_csv(file.path(outdir, 'fraction_activated_per_disease_and_state.csv'))

# These samples are dropped due to a low number of activated T-cells
#
dropped_samples <- alldf %>% 
  # bind_rows(gdf) %>%
  inner_join(dplyr::select(samples, unique_name, sample_group, sample_category), .) %>% 
  filter(CD38 > 1.5 & CD45RA < 3) %>% # & HLA_DR > 0) %>% 
  add_count(unique_name) %>% 
  filter(sample_category != 'tetramer+' & n < 800) %>% 
  distinct(unique_name, sample_group)
assert_that(nrow(dropped_samples) == 0)

alldf %>% 
  # bind_rows(gdf) %>% 
  inner_join(dplyr::select(samples, unique_name, sample_group, sample_category), .) %>% 
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
  write_csv(file.path(outdir, 'all_activated_tcells.csv.gz'))
## not correct as long as I remove samples with low counts
# print(glue("{nrow(adf)} out of {nrow(rdf)} are activated. ({round((nrow(adf) / nrow(rdf)) * 100, digits = 1)}%)"))
adf %>% 
  filter(sample_category == 'full') %>% 
  dplyr::count(sample_group, unique_name) %>% 
  identity() ->
  act_cells_per_sample
  
act_cells_per_sample %>% 
  summarize(median = median(n), mean = mean(n))

act_cells_per_sample %>% 
  group_by(sample_group) %>% 
  summarize(median = median(n), mean = mean(n))
median_cells <- median(act_cells_per_sample$n)
act_cells_per_sample %>% 
  ggplot(aes(sample_group, n)) +
  geom_boxplot()
ggsave(file.path(outdir, 'activated_cells_per_condition.png'))

## For clustering:
## sampling at most median(n) (3683) cells per sample
##
sample_at_most <- function (x, ncells) {
  set.seed(1234)
  sample_n(x, min(nrow(x), ncells))
}
sub_adf <- adf %>% 
  group_by(unique_name) %>%
  do(sample_at_most(., median_cells))

sub_adf %>% 
  write_csv(file.path(outdir, 'clust_subsampled_activated_tcells.csv.gz'))

## For t-SNE visualization
## sampling at most 807 cells per sample (except tetp)
##
tsne_adf <- sub_adf %>% 
  group_by(unique_name) %>%
  do(sample_at_most(., 807))
tsne_adf %>% 
  write_csv(file.path(outdir, 'tsne_input_subsampled_activated_tcells.csv.gz'))
