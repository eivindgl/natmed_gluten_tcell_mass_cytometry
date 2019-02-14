# Display donor specific plots of absoulte expression and tet+ / tet- fold change
pacman::p_load(
  tidyverse,
  devtools,
  viridis
)

source('20_experiment_and_sketch/load_fcs_and_rename_according_to_panel.R')

ac_markers <- read_csv('input_data/original/asbjorn_markers_specific.csv') %>% 
  dplyr::filter(category == 'CyTOF') %>% 
  dplyr::select(-category, -gene_name) %>% 
  dplyr::rename(marker = common_name)
ac_markers <- ac_markers %>% 
  mutate(name = str_replace_all(marker, ' ', '_'),
         name = str_replace_all(name, '-', '_'))

all_samples <- read_csv('input_data/renamed/sample_meta_full.csv')
hc_ucd_samples_meta <- all_samples %>% 
  filter(disease == 'Ced' | disease == 'HC') %>% 
  filter(disease_status %in% c('UCD', 'Normal')) %>% 
  filter(sample_category != 'AutoPhenotype' ) %>% 
  filter(donor != 'CD1535') %>% # challenge patient I should have annotated better ...
  filter(!(biosource == 'PBMC' & sample_category == 'tetramer-')) %>% 
  filter(donor != 'BC11' | str_detect(path, '170419 UCDpat 1 PBM')) %>% 
  filter(is.na(note) | note != 'exclude_ac') %>% 
  filter(!(biosource == 'SCS' & donor == 'CD1414')) %>% # no tetpos
  filter(!(biosource == 'PBMC' & donor == '1570')) %>%  # no tetpos
  mutate(sample_category = if_else(sample_category == 'full', 'tetramer-', sample_category)) %>% 
  arrange(disease, donor, sample_category)
all_fcs <- read_all_fcs_files(sample_meta = hc_ucd_samples_meta, fcs_trans = list(cofactor = 1, func = identity))
fcs <- all_fcs$df
fs <- flowset_of_subgroup(fcs)
sample_ids <- rep(sampleNames(fs), fsApply(fs, nrow))
expr <- fsApply(fs, exprs)

dr <- data.frame(
  sample = sample_ids,
  expr, stringsAsFactors = FALSE
  ) %>% 
  as_data_frame() %>% 
  inner_join(
    select(fcs, sample=unique_name, disease, disease_status, donor, sample_category, biosource)
  )
# dr %>% 
dr %>% 
  rownames_to_column('cell_id') %>% 
  mutate(disease_category = case_when(disease == 'HC' ~ 'Healthy Control',
                                      disease_status == 'GFD' ~ "Celiac - GFD",
                                      disease_status == 'UCD' ~ "Celiac - Untreated"),
         disease_category = factor(disease_category, c('Celiac - Untreated', 'Celiac - GFD', 'Healthy Control'))) %>% 
  gather(variable, value, -cell_id, -disease, -disease_status, -sample_category, -sample, -donor, -biosource, -disease_category) %>% 
  identity() ->
  tdr


tdr %>% 
  # filter(disease_status == 'UCD') %>% 
  group_by(biosource, disease, disease_category, donor, sample_category, variable) %>% 
  summarize(value = mean(value)) %>% 
  ungroup() %>% 
  identity() ->
  tdr_mean

tdr_mean %>% 
  write_csv('out/mean_marker_expression_per_donor.csv')

tdr_mean %>% 
  inner_join(ac_markers, by = c('variable' = 'name')) %>% 
  select(-variable) %>% 
  dplyr::rename(variable = marker) %>% 
  mutate(variable = factor(variable, rev(ac_markers$marker))) %>% 
  identity() ->
  tdr_mean

tdr_mean %>% 
  filter(biosource == 'PBMC') %>% 
  ggplot(aes(donor, variable)) +
  geom_tile(aes(fill = value)) +
  facet_grid(sample_category ~ disease_category, scales = 'free') +
  scale_fill_gradient(low = 'grey20', high = 'yellow2', limits = c(0, 75), oob = scales::squish) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = -45)) +
  labs(title = 'Donor specific absolute expression - blood')
ggsave('out/reviewer_questions/donor_specific_abs_expression_blood.png')
ggsave('out/reviewer_questions/donor_specific_abs_expression_blood.pdf')

tdr_mean %>% 
  filter(biosource == 'SCS') %>% 
  ggplot(aes(donor, variable)) +
  geom_tile(aes(fill = value)) +
  facet_grid(sample_category ~ disease_category, scales = 'free') +
  scale_fill_gradient(low = 'grey20', high = 'yellow2', limits = c(0, 75), oob = scales::squish) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = -45, hjust = -.1)) +
  labs(title = 'Donor specific absolute expression - gut')
ggsave('out/reviewer_questions/donor_specific_abs_expression_gut.png')
ggsave('out/reviewer_questions/donor_specific_abs_expression_gut.pdf')
  
tdr_mean %>% 
  filter(biosource == 'PBMC') %>% 
  filter(disease == 'Ced') %>% 
  spread(sample_category, value) %>% 
  na.omit() %>% 
  mutate(fold_change = log2(`tetramer+` / `tetramer-`)) %>% 
  identity() ->
  x

x %>% 
  select(donor, variable, fold_change) %>% 
  spread(variable, fold_change) %>% 
  as.data.frame() %>% 
  column_to_rownames('donor') %>% 
  as.matrix %>% 
  dist() %>% 
  head()
  hclust() %>% 
  # (function(z)  hclust(dist(z))) %>% 
  (function(z)  z$labels[z$order]) %>% 
  identity() ->
  donor_order

x %>% 
  mutate(donor = factor(donor, rev(donor_order))) %>%
  ggplot(aes(donor, variable)) +
  geom_tile(aes(fill = fold_change)) +
  scale_fill_gradient2(low = 'dodgerblue4', mid = 'grey20', high = 'yellow2', limits = c(-3, 3), oob=scales::squish) +
  theme_minimal() +
  labs(title = 'Donor specific log2 fold change - blood')
ggsave('out/reviewer_questions/donor_specific_fold_change_blood.png')
ggsave('out/reviewer_questions/donor_specific_fold_change_blood.pdf')

tdr_mean %>% 
  filter(biosource == 'SCS') %>% 
  filter(disease == 'Ced') %>% 
  mutate(value = pmax(value, 1e-7)) %>% 
  spread(sample_category, value) %>% 
  mutate(fold_change = log2(`tetramer+` / `tetramer-`)) %>% 
  identity() ->
  x

x %>% 
  # filter(is.infinite(fold_change))
  mutate(fold_change = if_else(is.infinite(fold_change), NaN, fold_change)) %>% 
  select(donor, variable, fold_change) %>% 
  spread(variable, fold_change) %>% 
  as.data.frame() %>% 
  column_to_rownames('donor') %>% 
  as.matrix %>% 
  dist() %>% 
  hclust() %>% 
  {.$labels[.$order]} %>% 
  identity() ->
  donor_order

x %>% 
  mutate(donor = factor(donor, rev(donor_order))) %>%
  ggplot(aes(donor, variable)) +
  geom_tile(aes(fill = fold_change)) +
  scale_fill_gradient2(low = 'dodgerblue4', mid = 'grey20', high = 'yellow2', limits = c(-3, 3), oob=scales::squish) +
  theme_minimal() +
  facet_grid(~ disease_category, scales = 'free') +
  theme(axis.text.x = element_text(angle = -45, hjust = -.1)) +
  labs(title = 'Donor specific log2 fold change - gut')
ggsave('out/reviewer_questions/donor_specific_fold_change_gut.png')
ggsave('out/reviewer_questions/donor_specific_fold_change_gut.png')

##
## Rev 3 question 2 - fold change between tetramer+ and healthy controls
##
tdr %>% 
  filter(sample_category == 'full' | sample_category == 'tetramer-') %>% 
  filter(disease_status != 'GFD') %>% 
  group_by(biosource, variable, disease, donor) %>% 
  summarize(value = mean(value)) %>% 
  summarize(value = mean(value)) %>% 
  ungroup() %>% 
  spread(disease, value) %>% 
  mutate(fold_change = log2(Ced / HC), dummy = '') %>% 
  identity() ->
  bg_grand_mean_df

bg_grand_mean_df %>% 
  filter(variable %in% ac_markers$name) %>% 
  mutate(variable = factor(variable, rev(ac_markers$name))) %>% 
  identity() ->
  bg_grand_mean_df

p_bg_grand_mean <- bg_grand_mean_df %>% 
  filter(biosource == 'PBMC') %>% 
  ggplot(aes(biosource, variable)) +
  geom_tile(aes(fill = fold_change)) +
  scale_fill_gradient2(low = 'dodgerblue4', mid = 'grey20', high = 'yellow2', limits = c(-4, 4), oob=scales::squish) +
  theme_minimal() +
  labs(title = glue::glue('Log2 fold change between tetramer- cells from CeD patients
       and all cells from healthy controls'))
  # facet_grid(~ biosource, scales = 'free') +
p_bg_grand_mean  
ggsave(filename = 'out/reviewer_questions/ced_vs_ctrl_tetneg_fold_change.png', plot = p_bg_grand_mean)
ggsave(filename = 'out/reviewer_questions/ced_vs_ctrl_tetneg_fold_change.pdf', plot = p_bg_grand_mean)
  

p_bg_grand_mean_gut <- bg_grand_mean_df %>% 
  filter(biosource == 'SCS') %>% 
  ggplot(aes(biosource, variable)) +
  geom_tile(aes(fill = fold_change)) +
  scale_fill_gradient2(low = 'dodgerblue4', mid = 'grey20', high = 'yellow2', limits = c(-3, 3), oob=scales::squish) +
  theme_minimal() +
  labs(title = glue::glue('Log2 fold change between tetramer- cells from CeD patients
       and all cells from healthy controls'))
  # facet_grid(~ biosource, scales = 'free') +
p_bg_grand_mean_gut
ggsave(filename = 'out/reviewer_questions/ced_vs_ctrl_tetneg_fold_change_gut.png', plot = p_bg_grand_mean_gut)
ggsave(filename = 'out/reviewer_questions/ced_vs_ctrl_tetneg_fold_change_gut.pdf', plot = p_bg_grand_mean_gut)
