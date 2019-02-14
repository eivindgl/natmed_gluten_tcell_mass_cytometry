# Used to migrate from old monolith system where filenames where parsed and matched to sample database.
# Should not be exectued for normal usage
#
pacman::p_load(
  tidyverse,
  lubridate,
  janitor,
  glue
)
df <- read_csv('input_data/renamed/fcs_matched_samples.csv') %>%
  mutate(outdir = dirname(path),
         filename = basename(path)) %>%
  select(-path) %>% 
  mutate(date = lubridate::ymd(`Experimet Date`)) %>% 
  mutate(Instrument = if_else(Instrument == 'Bergen', 'Helios Bergen', Instrument)) %>% 
  select(-`Experimet Date`) %>% 
  clean_names() %>% 
  select(-note, note)

write_meta_csv = function(x) {
  path <- file.path(unique(x$outdir), 'samples_meta.csv')
  x %>% 
    select(-outdir) %>% 
    write_csv(path)
  return(x)
}
  
df %>% 
  group_by(outdir) %>% 
  do(write_meta_csv(.))
