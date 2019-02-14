# Unsure how to best construct training set.
# I would ideally like to test on donors not used for training...
# however, we have few donors and some variablity.
# It might make more sense to train on subset of cells / bootstrapping..
# anyway, first problem is to get the training data set.

pacman::p_load(
  tidyverse
)
common_markers <- read_csv('input_data/common_markers_for_pred.csv') %>% 
  mutate(Rname = str_replace_all(Antibody, "[- ]", "_"))

all_samples <- read_csv('input_data/renamed/sample_meta_full.csv')
all_samples %>% 
  filter(sample_group == 'PBMC_UCD') %>% 
  distinct(donor) %>% 
  inner_join(all_samples) %>% 
  filter(str_detect(sample_category, 'tetramer')) %>% # skip PRE-enriched fraction
  identity() ->
  pbmc_ucd_tetramer_samples 
  
cytof_panels <- read_csv('input_data/original/cytof_panels.csv')
cytof_panels %>% 
  distinct(Panel, Antibody) %>% 
  add_count(Antibody) %>% 
  distinct(Antibody, n) %>% 
  filter(n == 2) %>% 
  select(-n) %>% 
  arrange(Antibody) %>% 
  write_csv('input_data/common_markers_for_pred.csv')

source('20_experiment_and_sketch/load_fcs_and_rename_according_to_panel.R')
fcs <- read_all_fcs_files(sample_meta = pbmc_ucd_tetramer_samples, fcs_trans = list(cofactor = 5, func = asinh))
x <- flowset_of_subgroup(fcs$df, markers = common_markers$Rname)

df <- data.frame(unique_name = rep(pbmc_ucd_tetramer_samples$unique_name, fsApply(x, nrow)),
           fsApply(x, exprs), stringsAsFactors = FALSE) %>% 
  as_tibble()

pbmc_ucd_tetramer_samples %>% 
  select(unique_name, sample_category) %>% 
  inner_join(df) %>% 
  identity() ->
  df
dir.create('out/predict_tetp', showWarnings = FALSE)  
df %>% 
  write_csv('out/predict_tetp/pbmc_ucd_tetp_and_tetn.csv')

## write full data set
all_fcs <- read_all_fcs_files(fcs_trans = list(cofactor = 5, func = asinh))
allx <- flowset_of_subgroup(all_fcs$df, markers = common_markers$Rname)

alldf <- data.frame(unique_name = rep(all_samples$unique_name, fsApply(allx, nrow)),
           fsApply(allx, exprs), stringsAsFactors = FALSE) %>% 
  as_tibble()

all_samples %>% 
  select(unique_name, sample_category) %>% 
  inner_join(alldf) %>% 
  identity() ->
  alldf

alldf %>% 
  write_csv('out/predict_tetp/all.csv.gz')
