pacman::p_load(
  tidyverse,
  fs,
  fuzzyjoin,
  lubridate
)
sdir <- 'input_data/original/Autoimmune control and Flu CD4/Autoimmune and control'
meta_path <- file.path(sdir, file.path('samples_meta.csv'))

panno <- read_csv('input_data/original/Autoimmune control and Flu CD4/parsed_asbjorn_annotation.csv')
panno
tibble(filename = list.files(sdir, pattern = '*.fcs')) %>% 
  fuzzyjoin::fuzzy_inner_join(panno, ., by = c('donor' = 'filename'), match_fun = function(xs, ys) {
  # print(xs)
  # print(ys)
  map2_lgl(xs, ys, partial(grepl, fixed = TRUE))
}) %>% 
  mutate(sample = path_ext_remove(filename),
         sample = as.character(sample),
         sample = tidy_names(sample, syntactic = TRUE),
         date = lubridate::as_date(date),
         disease_status = 'normal',
         sample_category = if_else(str_detect(filename, regex('phenotype', ignore_case = TRUE)), 'AutoPhenotype', 'full')
  ) %>% 
  identity() ->
  parsed_df
df %>% 
  anti_join(parsed_df, by = 'filename')
parsed_df %>% 
  anti_join(df, by = 'filename')
parsed_df %>% 
  write_csv(meta_path)

parsed_df %>% 
  mutate(sample_category = if_else(str_detect(filename, regex('phenotype', ignore_case = TRUE)), 'AutoPhenotype', 'full'))
  count(sample_category)


head(parsed_df$filename)  
