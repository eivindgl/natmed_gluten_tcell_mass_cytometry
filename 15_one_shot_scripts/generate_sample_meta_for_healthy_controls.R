pacman::p_load(
  tidyverse,
  fs,
  fuzzyjoin,
  lubridate,
  assertthat,
  glue
)
sdir <- 'input_data/original/PBMC/all_healthy_controls'
meta_path <- file.path(sdir, file.path('samples_meta.csv'))

panno <- read_csv('input_data/original/Autoimmune control and Flu CD4/parsed_asbjorn_annotation.csv') %>% 
  filter(category == 'healthy_control')

filenames <- list.files(sdir, pattern = '*.fcs')
tibble(filename = filenames) %>% 
  fuzzyjoin::fuzzy_inner_join(panno, ., by = c('donor' = 'filename'), match_fun = function(xs, ys) {
  # print(xs)
  # print(ys)
  map2_lgl(xs, ys, partial(grepl, fixed = TRUE))
}) %>% 
  mutate(sample = path_ext_remove(filename),
         sample = as.character(sample),
         sample = tidy_names(sample, syntactic = TRUE),
         date = lubridate::as_date(date),
         disease = 'HC',
         disease_status = 'normal',
         sample_category = if_else(str_detect(filename, regex('phenotype', ignore_case = TRUE)), 'AutoPhenotype', 'full')
  ) %>% 
  identity() ->
  parsed_df

# some donors are substrings of others......
# e.g. C4 BC40
parsed_df %>% 
  add_count(filename) %>% 
  filter(n > 1) %>% 
  arrange(filename)
# try to filter out missmatches
parsed_df %>% 
  add_count(filename) %>% 
  mutate(inferred_donor = str_detect(filename, glue("^{donor}"))) %>% 
  filter(n == 1 | inferred_donor) %>% 
  identity() ->
  pf_df

nrow(parsed_df)
nrow(pf_df)
length(filenames)
tibble(filename = filenames) %>% 
  anti_join(pf_df) 
# currently only missing donor `908 Control``
# manually add later
# assert_that(nrow(parsed_df) == length(filenames))

pf_df %>% 
  select(-inferred_donor, -n) %>% 
  write_csv(meta_path)

