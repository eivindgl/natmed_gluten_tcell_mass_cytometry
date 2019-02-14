# identifies missing donors between my curated list and asbjorns notes

pacman::p_load(
  tidyverse,
  readxl,
  janitor,
  glue
)

edf <- read_csv('input_data/renamed/sample_meta_full.csv')
adf <- read_csv('input_data/original/asbjorn_PBMC list Celiac disease.csv') %>% 
  clean_names() %>% 
  dplyr::rename(orig_name = sample) %>% 
  mutate(donor = if_else(str_detect(orig_name, 'CD'),
                         str_extract(orig_name, 'CD[[:number:]]+'),
                         str_extract(orig_name, '[[:alpha:]]+[^ ]+'))) %>% 
  mutate(donor= str_replace(donor, 'LRS', 'BC'))
  
# adf <- read_excel('input_data/original/asbjorn_PBMC list Celiac disease.xlsx') %>% 
#   clean_names() %>% 
#   select(1,2) %>% 
#   na.omit() %>% 
#   select(-1) %>% 
#   rename(orig_name = patients_included_in_the_primary_cy_tof_study) %>% 
#   mutate(donor = if_else(str_detect(orig_name, 'CD'),
#                         str_extract(orig_name, 'CD[[:number:]]+'),
#                         str_extract(orig_name, '[[:alpha:]]+[^ ]+'))) %>% 
#   mutate(donor= str_replace(donor, 'LRS', 'BC')) %>% 
#   identity() ->
#   adf

edf %>% 
  distinct(donor, biosource, disease) %>% 
  full_join(adf) %>% 
  print(n = Inf) %>% 
  filter(is.na(biosource)) %>% # donors not in my panel
  glue_data("donor = {donor} ({orig_name} - {sample_type} - {celiac_or_other})")
edf
