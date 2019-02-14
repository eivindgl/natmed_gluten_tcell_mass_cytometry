pacman::p_load(
  tidyverse,
  assertthat,
  glue
)

#' Adds a column of publication relevant group name to input sample meta data frame.
add_sample_groups <- function(all_samples) {
  assert_that(has_name(all_samples, 'unique_name'))
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
    filter(category == 'autoimmune')
  
  ## GUT
  
  sample_groups$SCS$HC <- all_samples %>% 
    filter(biosource == 'SCS') %>% 
    filter(disease == 'HC')
  
  sample_groups$SCS$UCD <- all_samples %>% 
    filter(biosource == 'SCS') %>% 
    filter(disease == 'Ced') %>% 
    filter(disease_status == 'UCD')
  
  sample_groups$SCS$GFD <- all_samples %>% 
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
  sgdf
}


mdf <- data_frame(metapath=list.files('input_data/original', pattern = 'samples_meta.csv', 
                               recursive = TRUE, full.names = TRUE)) %>% 
  mutate(src_dir = dirname(metapath))
mdf$metapath %>% 
  set_names(., .) %>% 
  map_dfr(read_csv, .id = 'metapath') %>% 
  inner_join(mdf) %>% 
  select(-metapath) %>% 
  identity() ->
  df
df %>% 
  mutate(path = file.path(src_dir, filename),
         unique_name = str_c(sample, sample_category, sep = '_')) %>% 
  select(-src_dir, -filename) %>% 
  identity() ->
  fdf

fdf %>% 
  group_by(unique_name) %>% 
  filter(n() == 1) %>% 
  ungroup() %>% 
  identity() ->
  unique_samples

fdf %>% 
  group_by(unique_name) %>% 
  filter(n() > 1) %>% 
  ungroup() %>% 
  identity() ->
  dup_samples
print(glue("{nrow(fdf)} samples in total. {nrow(unique_samples)} are already unique. {nrow(fdf) - nrow(unique_samples)} remaining..."))
dup_samples %>% 
  mutate(name = basename(path)) %>% 
  select(sample, unique_name, name) %>% 
  print(n = Inf)

dup_samples %>% 
  group_by(unique_name) %>% 
  rownames_to_column() %>% 
  mutate(tmp = str_c(unique_name, '_N', rowname)) %>% 
  ungroup() %>% 
  select(-rowname, -unique_name) %>% 
  rename(unique_name = tmp) %>% 
  identity() ->
  dup_samples
dup_samples %>% 
  select(unique_name, everything()) %>% 
  print(n = Inf)

assert_that(nrow(dup_samples) == length(unique(dup_samples$unique_name)))
bind_rows(
  unique_samples,
  dup_samples
) %>% 
  add_sample_groups() %>% 
  identity() ->
  dfinal
dfinal %>% 
  write_csv('input_data/renamed/sample_meta_full.csv')
