# We have two panels of antibodies and metal combinations.
# However, it turns out that the metal naming within each panel type (between FCS files) are sometimes inconsistent.
# I want to get an overview of the differences.
pacman::p_load(
  tidyverse,
  flowCore,
  assertthat,
  glue,
  stringr,
  assertthat
)

filter <- dplyr::filter
select <- dplyr::select

sample_meta <- read_csv('input_data/renamed/fcs_matched_samples.csv')
cytof_panels <- read_csv('input_data/original/cytof_panels.csv')
ab_type <- read_csv('input_data/original/antibody_meta.csv')


all_fcs <- sample_meta$path %>% 
  map(~ read.FCS(., transformation = FALSE, truncate_max_range = FALSE))
names(all_fcs) <- sample_meta$path

extract_metals_and_ab_df <- function(fcs) {
  desc <- fcs %>% 
    parameters() %>% 
    magrittr::extract2('desc') %>% 
    as.character()
  metals <- colnames(fcs)
  metmatch <- str_match(metals, '^([[:alpha:]]+)([[:digit:]]+)')
  abmatch <- str_match(desc, '^[[:digit:]]+[[:alpha:]]+_(.+)')
  
  metal_number <- metmatch[, 3] # E.g. Dy160 -> 160 should be a unique identifier for this metal
  asbjorn_version <- str_c(metmatch[, 3], metmatch[, 2])
  assert_that(all(!duplicated(na.omit(asbjorn_version)))) # All metal numbers should be unique
  tibble(metal_r = metals,
         metal_number = metal_number,
         metal_asbjorn = asbjorn_version,
         met_ab = desc,
         ab = abmatch[, 2]) %>% 
    na.omit()
}

dir.create('out/tmp', showWarnings = FALSE, recursive = TRUE)
all_fcs %>% 
  map_dfr(extract_metals_and_ab_df, .id = 'path') %>% 
  filter(ab == 'empty')
  count(metal_asbjorn) %>% 
  arrange(n) %>% 
  print(n = 20)
all_fcs %>% 
  map_df(extract_metals_and_ab_df) %>% 
  count(metal_asbjorn, metal_r, ab) %>% 
  arrange(n) %>% 
  identity() ->
  ab_counts

ab_counts %>% 
  distinct(ab) %>% 
  left_join(ab_type, by = c('ab' = 'Antibody')) %>% 
  mutate(Antibody = if_else(is.na(lineage), '', ab)) %>% 
  select(fcs_antibody = ab, Antibody) %>% 
  write_csv('input_data/original/fcs_AB_to_Antibody.csv')

ab_counts %>% 
  filter(ab == 'LD')
