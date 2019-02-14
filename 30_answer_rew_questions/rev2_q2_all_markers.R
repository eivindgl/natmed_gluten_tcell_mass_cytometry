##
## Currently only for PBMCs
##
#      Boxplots or table displaying individual frequency data for all markers used in main figures
# among total CD4+ T cells, Tet+ and Tet- CD4+ T cells should be shown for untreated and treated celiac disease patients.
# The analysis as presented does not allow to understand whether particular surface protein
# drive the TA phenotype and what the variance in expression is for these
# proteins among different individuals.
pacman::p_load(
  tidyverse,
  devtools,
  viridis,
  assertr
)
##
## Ensure we use the right dplyr functions
##
select <- dplyr::select
filter <- dplyr::filter
count <- dplyr::count

source('20_experiment_and_sketch/load_fcs_and_rename_according_to_panel.R')

select_fcs_subset <- function(xs, markers) {
  f <- function(x) {
    expr <- exprs(x)
    exprs(x) <- expr[, markers, drop = F]
    x
  }
  xs %>% 
    map(f)
}

all_samples <- read_csv('input_data/renamed/sample_meta_full.csv')
ced_samples_meta <- all_samples %>% 
  filter(disease == 'Ced') %>% 
  # filter(disease_status %in% c('UCD', 'GFD')) %>% # Skipping GFD for by marker
  filter(disease_status == 'UCD') %>% 
  filter(donor != 'CD1535') %>% # challenge patient I should have annotated better ...
  filter(!(biosource == 'PBMC' & sample_category == 'tetramer-')) %>% 
  filter(donor != 'BC11' | str_detect(path, '170419 UCDpat 1 PBM')) %>% 
  filter(!(biosource == 'SCS' & donor == 'CD1414')) %>% # no tetpos
  filter(!(biosource == 'PBMC' & donor == 'CD1570')) %>%  # no tetpos
  mutate(sample_category = if_else(sample_category == 'full', 'tetramer-', sample_category)) %>% 
  arrange(disease, donor, sample_category)


ac_markers <- read_csv('input_data/original/asbjorn_markers_specific.csv') %>% 
  dplyr::filter(category == 'CyTOF') %>% 
  dplyr::select(-category, -gene_name) %>% 
  dplyr::rename(marker = common_name)
ac_markers <- ac_markers %>% 
  mutate(name = str_replace_all(marker, ' ', '_'),
         name = str_replace_all(name, '-', '_'))

all_fcs <- read_all_fcs_files(sample_meta = ced_samples_meta)
ced_fcs <- all_fcs$df  
  # special case (missing CD28 and CD57. Needs to be processed in isolation and joined afterwards to tdr)
  # filter(!(biosource == 'SCS' & donor == 'CD1414'))

fcs_as_tdf <- function(fcs_df) {
  fsced <- flowset_of_subgroup(fcs_df)#, markers = ac_markers$name)
  assert_that('OX40' %in% colnames(fsced)) # ensure that we only use CeD panel
  assert_that('CD137' %in% colnames(fsced)) # this marker was left out for some reason
  sample_ids <- rep(fcs_df$unique_name, fsApply(fsced, nrow))
  expr <- fsApply(fsced, exprs)

  dr <- data.frame(
    sample = sample_ids,
    expr, stringsAsFactors = FALSE
  ) %>% 
    as_data_frame() %>% 
    inner_join(
      select(fcs_df, sample=unique_name, disease, disease_status, donor, sample_category, biosource)
    ) %>% 
    rownames_to_column('cell_id') %>% 
    gather(variable, value, -cell_id, -disease, -disease_status, -sample_category, -sample, -donor, -biosource) %>% 
    filter(variable %in% ac_markers$name) %>% 
    mutate(variable = factor(variable, ac_markers$name, ac_markers$marker),
           sample_category = if_else(sample_category == 'tetramer+', 'tetramer+', 'All CD4'))
}
tdr <- fcs_as_tdf(ced_fcs)

tdr %>% 
  distinct(variable) %>%  
  rename(marker = variable) %>% 
  anti_join(ac_markers, .) %>% 
  verify(nrow(.) == 0)

###
### Boxplots by markers with donors on the x-axis
###
# p_blood_ucd <- 
pacman::p_load(
  ggforce
)

# xdf %>% 
p_pbmc_ucd_all_markers <- tdr %>% 
  filter(disease_status == 'UCD' & biosource == 'PBMC') %>% 
  ggplot(aes(donor, value, fill = sample_category)) +
  ## hollow outliers
  # geom_boxplot(aes(fill = sample_category), position = position_dodge2(preserve = 'total'), outlier.alpha = 0.1, outlier.shape = 1, outlier.fill = NULL) +
  geom_boxplot(aes(fill = sample_category), position = position_dodge2(preserve = 'total'), outlier.shape = NA) +
  # filter(donor == 'CD1349') %>%
  # theme_minimal() +
  theme_minimal() +
  theme(axis.text.x=element_text(angle = -45)) +
  labs(fill = 'Experiment',
       title = 'Boxplot of all CD4 markers',
       subtitle = 'PBMC UCD') +
  scale_fill_manual(values = c('tetramer+' = '#f36d16', 'All CD4' = '#3672bc')) +
  scale_y_continuous(limits = c(0, 5)) +
  facet_wrap(~ variable, ncol = 5, nrow = 6)
p_pbmc_ucd_all_markers

p_scs_ucd_all_markers <- tdr %>% 
  filter(disease_status == 'UCD' & biosource == 'SCS') %>% 
  ggplot(aes(donor, value, fill = sample_category)) +
  ## hollow outliers
  # geom_boxplot(aes(fill = sample_category), position = position_dodge2(preserve = 'total'), outlier.alpha = 0.1, outlier.shape = 1, outlier.fill = NULL) +
  geom_boxplot(aes(fill = sample_category), position = position_dodge2(preserve = 'total'), outlier.shape = NA) +
  # filter(donor == 'CD1349') %>%
  # theme_minimal() +
  theme_minimal() +
  theme(axis.text.x=element_text(angle = -45)) +
  labs(fill = 'Experiment',
       title = 'Boxplot of all CD4 markers',
       subtitle = 'SCS UCD') +
  scale_fill_manual(values = c('tetramer+' = '#f36d16', 'All CD4' = '#3672bc')) +
  scale_y_continuous(limits = c(0, 5)) +
  facet_wrap(~ variable, ncol = 5, nrow = 6)
ggsave(filename = 'out/reviewer_questions/boxplot_by_marker_gut_ucd.png', plot = p_scs_ucd_all_markers)
ggsave(filename = 'out/reviewer_questions/boxplot_by_marker_gut_ucd.pdf', plot = p_scs_ucd_all_markers)
ggsave(filename = 'out/reviewer_questions/boxplot_by_marker_blood_ucd.pdf', plot = p_pbmc_ucd_all_markers)
ggsave(filename = 'out/reviewer_questions/boxplot_by_marker_blood_ucd.png', plot = p_pbmc_ucd_all_markers)
