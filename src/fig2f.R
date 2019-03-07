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

ced_tp_samples <- read_csv('meta/meta.csv') %>% 
  filter(biosource == 'PBMC' & panel == 'CeD_panel') %>% 
  filter(disease == 'CeD' & sample_category == 'tetramer+')

gfd_ch_samples <- ced_tp_samples %>% 
  filter(disease_status == 'challenge') %>% 
  distinct(participant) %>% 
  inner_join(ced_tp_samples)

ucd_samples <- ced_tp_samples %>% 
  filter(disease_status == 'ucd')

samples <- bind_rows(gfd_ch_samples, ucd_samples)

fcs_df <- map_dfr(samples$filename, read_dataset)

groups <- samples %>% 
  select(filename, disease_status) %>% 
  transmute(filename, 
            group = tolower(disease_status))

gmean <- fcs_df %>% 
  inner_join(groups) %>% 
  gather(marker, expr, -filename, -group) %>% 
  na.omit() %>% 
  mutate(expr = sinh(expr)) %>% 
  group_by(group, marker, filename) %>% 
  summarize(expr = mean(expr)) %>% 
  summarize(expr = mean(expr)) %>% 
  ungroup()

challenge_logfc <- gmean %>% 
  filter(group != 'ucd') %>% 
  inner_join(markers) %>% 
  select(antibody, group, expr) %>% 
  spread(group, expr) %>% 
  transmute(antibody, group = 'Tet+ challenge vs Tet+ TCeD', logfc = log2(`challenge` / `gfd`))

ch_ucd_logfc <- gmean %>% 
  filter(group != 'gfd') %>% 
  inner_join(markers) %>% 
  select(antibody, group, expr) %>% 
  spread(group, expr) %>% 
  transmute(antibody, group = 'Tet+ challenge vs Tet+ UCeD', logfc = log2(`challenge` / `ucd`))


bind_rows(challenge_logfc, ch_ucd_logfc) %>% 
  mutate(antibody = factor(antibody, levels = rev(markers$antibody)),
         group = fct_relevel(group, 'Tet+ challenge vs Tet+ TCeD')) %>% 
  ggplot(aes(group, antibody)) +
  geom_tile(aes(fill = logfc)) +
  scale_fill_gradient2(low = 'dodgerblue4', mid = 'grey20', high = 'yellow2', limits = c(-4, 4), oob=scales::squish) +
  coord_equal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave('out/plots/fig2f_pbmc_logfc.png')
