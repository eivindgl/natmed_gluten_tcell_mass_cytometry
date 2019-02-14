##
## Currently only for PBMCs
##
#      Boxplots or table displaying individual frequency data for all TA markers
# (CD45RA-, CD62L-, CXCR3+, CD39+, CD38+, PD-1, CD127low, CD25-, ICOS+, CD161+)
# among total CD4+ T cells, Tet+ and Tet- CD4+ T cells should be shown for untreated and treated celiac disease patients.
# The analysis as presented does not allow to understand whether particular surface protein
# drive the TA phenotype and what the variance in expression is for these
# proteins among different individuals.
pacman::p_load(
  tidyverse,
  devtools,
  viridis
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
  filter(disease_status %in% c('UCD', 'GFD')) %>% 
  filter(donor != 'CD1535') %>% # challenge patient I should have annotated better ...
  filter(!(biosource == 'PBMC' & sample_category == 'tetramer-')) %>% 
  filter(donor != 'BC11' | str_detect(path, '170419 UCDpat 1 PBM')) %>% 
  mutate(sample_category = if_else(sample_category == 'full', 'tetramer-', sample_category)) %>% 
  arrange(disease, donor, sample_category)

all_fcs <- read_all_fcs_files(sample_meta = ced_samples_meta)
ced_fcs <- all_fcs$df

fsced <- flowset_of_subgroup(ced_fcs)
assert_that('OX40' %in% colnames(fsced)) # ensure that we only use CeD panel
sample_ids <- rep(ced_fcs$unique_name, fsApply(fsced, nrow))
expr <- fsApply(fsced, exprs)

dr <- data.frame(
  sample = sample_ids,
  expr, stringsAsFactors = FALSE
  ) %>% 
  as_data_frame() %>% 
  inner_join(
    select(ced_fcs, sample=unique_name, disease, disease_status, donor, sample_category, biosource)
  )
# dr %>% 
dr %>% 
  rownames_to_column('cell_id') %>% 
  gather(variable, value, -cell_id, -disease, -disease_status, -sample_category, -sample, -donor, -biosource) %>% 
  identity() ->
  tdr

TA_markers <- c('CD45RA', 'CD62L', 'CXCR3', 'CD39', 'CD38', 'PD_1', 'CD127', 'CD25', 'ICOS', 'CD161')
p_blood_ucd <- tdr %>% 
  filter(disease_status == 'UCD' & biosource == 'PBMC') %>% 
  filter(variable %in% TA_markers) %>% 
  mutate(variable = factor(variable, TA_markers),
         sample_category = if_else(sample_category == 'tetramer+', 'tetramer+', 'All CD4')) %>% 
  # filter(donor == 'CD1349') %>%
  ggplot(aes(variable, value)) +
  geom_boxplot(aes(fill = sample_category), position = position_dodge2(preserve = 'total')) +
  facet_wrap(~ donor, ncol = 2) +
  # theme_minimal() +
  theme_minimal() +
  theme(axis.text.x=element_text(angle = -45)) +
  labs(fill = 'Experiment',
       title = 'Boxplot of all TA markers',
       subtitle = 'PBMC UCD')

p_blood_gfd <- tdr %>% 
  filter(disease_status == 'GFD' & biosource == 'PBMC') %>% 
  filter(variable %in% TA_markers) %>% 
  mutate(variable = factor(variable, TA_markers),
         sample_category = if_else(sample_category == 'tetramer+', 'tetramer+', 'All CD4')) %>% 
  # filter(donor == 'CD1349') %>%
  ggplot(aes(variable, value)) +
  geom_boxplot(aes(fill = sample_category), position = position_dodge2(preserve = 'total')) +
  facet_wrap(~ donor, ncol = 2) +
  # theme_minimal() +
  theme_minimal() +
  theme(axis.text.x=element_text(angle = -45)) +
  labs(fill = 'Experiment',
       title = 'Boxplot of all TA markers',
       subtitle = 'PBMC GFD')

dir.create('out/reviewer_questions', recursive = TRUE, showWarnings = FALSE)
ggsave(filename = 'out/reviewer_questions/boxplot_of_TA_markers_blood_ucd.png', plot = p_blood_ucd +
         scale_fill_manual(values = c('tetramer+' = '#f36d16', 'All CD4' = '#3672bc')), height = 11)
ggsave(filename = 'out/reviewer_questions/boxplot_of_TA_markers_blood_ucd.pdf', plot = p_blood_ucd +
         scale_fill_manual(values = c('tetramer+' = '#f36d16', 'All CD4' = '#3672bc')), height = 11)
ggsave(filename = 'out/reviewer_questions/boxplot_of_TA_markers_blood_gfd.png', plot = p_blood_gfd +
         scale_fill_manual(values = c('tetramer+' = '#f36d16', 'All CD4' = '#3672bc')), height = 11)
ggsave(filename = 'out/reviewer_questions/boxplot_of_TA_markers_blood_gfd.pdf', plot = p_blood_gfd +
         scale_fill_manual(values = c('tetramer+' = '#f36d16', 'All CD4' = '#3672bc')), height = 11)

p_gut_gfd <- tdr %>% 
  filter(variable %in% TA_markers) %>% 
  filter(biosource == 'SCS' & disease_status == 'GFD') %>% 
  mutate(variable = factor(variable, TA_markers),
         sample_category = if_else(sample_category == 'tetramer+', 'tetramer+', 'All CD4')) %>% 
  # filter(donor == 'CD1349') %>%
  ggplot(aes(variable, value)) +
  geom_boxplot(aes(fill = sample_category), position = position_dodge2(preserve = 'total')) +
  # theme_minimal() +
  facet_wrap(~ donor, ncol = 2) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle = -45)) +
  labs(fill = 'Experiment',
       title = 'Boxplot of all TA markers',
       subtitle = 'Gut GFD')

p_gut_ucd <- tdr %>% 
  filter(variable %in% TA_markers) %>% 
  filter(biosource == 'SCS' & disease_status == 'UCD') %>% 
  mutate(variable = factor(variable, TA_markers),
         sample_category = if_else(sample_category == 'tetramer+', 'tetramer+', 'All CD4')) %>% 
  # filter(donor == 'CD1349') %>%
  ggplot(aes(variable, value)) +
  geom_boxplot(aes(fill = sample_category), position = position_dodge2(preserve = 'total')) +
  # theme_minimal() +
  facet_wrap(~ donor, ncol = 2) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle = -45)) +
  labs(fill = 'Experiment',
       title = 'Boxplot of all TA markers',
       subtitle = 'Gut UCD')

dir.create('out/reviewer_questions', recursive = TRUE, showWarnings = FALSE)
ggsave(filename = 'out/reviewer_questions/boxplot_of_TA_markers_gut_gfd.png', plot = p_gut_gfd +
         scale_fill_manual(values = c('tetramer+' = '#f36d16', 'All CD4' = '#3672bc')), height = 11)
ggsave(filename = 'out/reviewer_questions/boxplot_of_TA_markers_gut_gfd.pdf', plot = p_gut_gfd +
         scale_fill_manual(values = c('tetramer+' = '#f36d16', 'All CD4' = '#3672bc')), height = 11)
ggsave(filename = 'out/reviewer_questions/boxplot_of_TA_markers_gut_ucd.png', plot = p_gut_ucd +
         scale_fill_manual(values = c('tetramer+' = '#f36d16', 'All CD4' = '#3672bc')), height = 11)
ggsave(filename = 'out/reviewer_questions/boxplot_of_TA_markers_gut_ucd.pdf', plot = p_gut_ucd +
         scale_fill_manual(values = c('tetramer+' = '#f36d16', 'All CD4' = '#3672bc')), height = 11)

