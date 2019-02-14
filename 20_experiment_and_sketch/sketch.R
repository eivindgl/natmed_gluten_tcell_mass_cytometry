pacman::p_load(
  tidyverse,
  flowCore
)
filter <- dplyr::filter

samplesdf <- read_csv('../input_data/renamed/fcs_matched_samples.csv')
samplesdf

samplesdf %>% 
  filter(Disease == 'Flu') %>% 
  head(1) %>% 
  magrittr::extract2('path') %>% 
  identity() ->
  test_fcs_path

cdf

cdf %>% 
  distinct(Antibody) %>% 
  write_csv('../input_data/original/antibody_meta.csv')

cdf <- read_csv('../input_data/original/antibody_meta.csv')
cdf



afcs <- extract_metals_and_ab_df(test_fcs_path)
cfcs <- extract_metals_and_ab_df('../input_data/original/SCS/170718 SCS Multiple pat for R/HCD1204D.A Tet POS.fcs')
x %>% 
  mutate(ml = st)
  str_extract()


adf <- samplesdf$path %>% 
  # head() %>% 
  map_df(extract_metals_and_ab_df)

adf %>% 
  count(met_ab) %>% 
  arrange(n) %>% 
  print(n = Inf)

cdf <- adf %>% 
  count(met_ab) %>% 
  arrange(n)

filter(cdf, n < 12)
cdf %>% 
  filter(str_detect(met_ab, 'CTLA'))

cdf %>% 
  filter(n < 77)

pdf <- read_csv('../input_data/original/cytof_panels.csv')
pdf %>% 
 # group_by(Metal) %>% 
  # mutate(idx = row_number()) %>% 
  spread(Panel, Antibody) %>% 
  filter(CeD != Autoimmune) %>% 
  identity() ->
  metals_that_varies

adf <- afcs %>% 
  inner_join(metals_that_varies, by = c('metal_a' = 'Metal'))

cdf <- cfcs %>% 
  inner_join(metals_that_varies, by = c('metal_a' = 'Metal'))


str_detect(adf$met_ab, adf$Autoimmune)
str_detect(cdf$met_ab, cdf$CeD)

df %>% 
  mutate(in_CeD = map_lgl(), 
         in_auto)


pdf %>% 
  add_rownames() %>% 
  filter(rowname %in% 30:32)
pdf


pdf %>% 
  distinct(Antibody) %>% 
  arrange(Antibody) %>% 
  mutate(lineage = 0, functional = 0) %>% 
  write_csv('../input_data/original/antibody_meta.csv')

adf <- read_csv('../input_data/original/antibody_meta.csv')
