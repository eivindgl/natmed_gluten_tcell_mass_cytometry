pacman::p_load(
  tidyverse,
  assertthat,
  glue,
  rhdf5
)

select <- dplyr::select
filter <- dplyr::filter
count <- dplyr::count

read_dataset <- function(unique_name, path = 'out/cytof.h5') {
  assert_that(is_bare_character(unique_name))
  h5_path <- glue('samples/{unique_name}') 
  print(glue("reading {h5_path}"))
  M <- h5read(path, h5_path)
  colnames(M) <- h5readAttributes(path, h5_path)$colnames
  M %>% 
    as_tibble() %>% 
    # gather(marker, value) %>% 
    mutate(filename = unique_name) %>% 
    select(filename , everything())
}


markers <- read_csv('meta/cytof_markers_specific.csv')

samples <- read_csv('meta/meta.csv') %>% 
  filter(biosource == 'PBMC' & panel == 'CeD_panel') %>% 
  filter(sample_category != 'tetramer-') %>% 
  filter(disease_status == 'ucd' | (disease == 'HC' & sample_category == 'full'))

fcs_df <- map_dfr(samples$filename, read_dataset)

groups <- samples %>% 
  select(filename, disease, sample_category) %>% 
  transmute(filename, 
            group = str_c(disease, sample_category, sep = '_'),
            group = tolower(group))

gmean <- fcs_df %>% 
  inner_join(groups) %>% 
  gather(marker, expr, -filename, -group) %>% 
  na.omit() %>% 
  mutate(expr = sinh(expr)) %>% 
  group_by(group, marker, filename) %>% 
  summarize(expr = mean(expr)) %>% 
  summarize(expr = mean(expr)) %>% 
  ungroup()

ced_logfc <- gmean %>% 
  filter(str_detect(group, '^ced')) %>% 
  inner_join(markers) %>% 
  select(antibody, group, expr) %>% 
  spread(group, expr) %>% 
  transmute(antibody, group = 'UCeD Tet+ vs UCeD CD4+', logfc = log2(`ced_tetramer+` / `ced_full`))

tphc_logfc <- gmean %>% 
  filter(group != 'ced_full') %>% 
  inner_join(markers) %>% 
  select(antibody, group, expr) %>% 
  spread(group, expr) %>% 
  transmute(antibody, group = 'UCeD Tet+ vs Ctr.CD4+', logfc = log2(`ced_tetramer+` / `hc_full`))


bind_rows(ced_logfc, tphc_logfc) %>% 
  mutate(antibody = factor(antibody, levels = rev(markers$antibody)),
         group = fct_relevel(group, 'UCeD Tet+ vs UCeD CD4+')) %>% 
  ggplot(aes(group, antibody)) +
  geom_tile(aes(fill = logfc)) +
  scale_fill_gradient2(low = 'dodgerblue4', mid = 'grey20', high = 'yellow2', limits = c(-4, 4), oob=scales::squish) +
  coord_equal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave('out/plots/fig2e_pbmc_logfc.png')
