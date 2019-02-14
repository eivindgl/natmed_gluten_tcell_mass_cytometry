pacman::p_load(
  tidyverse,
  assertthat
)

all_samples <- read_csv('input_data/renamed/sample_meta_full.csv')

all_samples %>% 
  distinct(sample_group, donor) %>% 
  count(sample_group)

sample_groups <- list(PBMC = list(), SCS = list())

## PBMC

sample_groups$PBMC$HC <- all_samples %>% 
  filter(biosource == 'PBMC') %>% 
  filter(disease == 'HC') %>% 
  filter(donor != 'BC11' | str_detect(path, '170419 UCDpat 1 PBM'))

sample_groups$PBMC$Flu <- all_samples %>% 
  filter(biosource == 'PBMC') %>% 
  filter(disease == 'Flu')

sample_groups$PBMC$UCD <- all_samples %>% 
  filter(biosource == 'PBMC') %>% 
  filter(disease == 'Ced') %>% 
  filter(disease_status %in% c('UCD')) %>% 
  filter(donor != 'CD1535') %>% # challenge patient I should have annotated better ...
  filter(!(biosource == 'PBMC' & sample_category == 'tetramer-')) %>% 
  filter(donor != 'BC11' | str_detect(path, '170419 UCDpat 1 PBM'))

sample_groups$PBMC$GFD <- all_samples %>% 
  filter(biosource == 'PBMC') %>% 
  filter(disease == 'Ced') %>% 
  filter(disease_status %in% c('GFD')) %>% 
  filter(donor != 'CD1535') %>% # challenge patient I should have annotated better ...
  filter(!(biosource == 'PBMC' & sample_category == 'tetramer-')) %>% 
  filter(donor != 'BC11' | str_detect(path, '170419 UCDpat 1 PBM'))

sample_groups$PBMC$AutoImmune <- all_samples %>% 
  filter(biosource == 'PBMC') %>% 
  filter(disease %in% c('SLE', 'DM1', 'SSc'))

## GUT

sample_groups$SCS$HC <- all_samples %>% 
  filter(biosource == 'SCS') %>% 
  filter(disease == 'HC')

sample_groups$SCS$UCD <- all_samples %>% 
  filter(biosource == 'SCS') %>% 
  filter(disease == 'Ced') %>% 
  filter(disease_status == 'UCD')

sample_groups$SCS$GFD<- all_samples %>% 
  filter(biosource == 'SCS') %>% 
  filter(disease == 'Ced') %>% 
  filter(disease_status == 'GFD')

unnested_sample_groups <- sample_groups %>% 
  unlist(use.names = TRUE, recursive = FALSE)

unnested_sample_groups %>% 
  bind_rows(.id = 'sample_group') %>% 
  identity() ->
  sgdf

unnested_sample_groups$ungrouped <-  all_samples %>% 
  anti_join(sgdf)

unnested_sample_groups %>% 
  bind_rows(.id = 'sample_group') %>% 
  mutate(sample_group = str_replace_all(sample_group, '\\.', '_')) %>% 
  dplyr::select(unique_name, sample_group, everything()) %>% 
  identity() ->
  sgdf


# nrow(sgdf)
# nrow(all_samples)
# sgdf %>%
#   add_count(unique_name) %>% 
#   filter(n > 1) %>% 
#   arrange(unique_name)
assert_that(nrow(sgdf) == nrow(all_samples))

all_samples %>% 
  distinct(sample_group, donor) %>% 
  count(sample_group)
