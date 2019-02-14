source('20_experiment_and_sketch/load_fcs_and_rename_according_to_panel.R')
# generate log fold change heatmaps for asbjorn's paper
# 
# Two experiments in blood and gut:
# tet+_vs_tet-
# tet+_vs_ctrl
#
# In addition, absolute expression of tet+ (either raw values or asinh)
pacman::p_load(
  flowCore,
  magrittr,
  tidyverse,
  stringr,
  matrixStats,
  assertthat,
  magrittr,
  forcats,
  cowplot,
  scales,
  svglite
)
egl_theme <- theme_minimal() +
  theme(legend.position='bottom') #, 
# axis.text.x=element_blank())
old_theme <- theme_set(egl_theme)

##
## Ensure we use the right dplyr functions
##
select <- dplyr::select
filter <- dplyr::filter
count <- dplyr::count

## Read antibodies to generate plots for
asbjorn_markers_all <- read_csv('input_data/original/asbjorn_markers_specific.csv')
asbjorn_markers <- asbjorn_markers_all %>% 
  filter(category == 'CyTOF') %>% 
  rename(antibody = common_name)



# ced_fcs <- all_fcs$df %>%
all_samples <- read_csv('input_data/renamed/sample_meta_full.csv')
ced_samples_meta <- all_samples %>% 
  filter(disease == 'Ced' | disease == 'HC') %>% 
  filter(disease_status %in% c('UCD', 'Normal')) %>% 
  filter(donor != 'CD1535') %>% # challenge patient I should have annotated better ...
  filter(! (donor == 'CD1414' & biosource == 'SCS')) %>% # no tetramer+ and no CD28 or CD57
  filter(!(biosource == 'PBMC' & sample_category == 'tetramer-')) %>% 
  filter(donor != 'BC11' | str_detect(path, '170421 UCDpat 2 PBM')) %>% 
  mutate(sample_category = if_else(sample_category == 'full', 'tetramer-', sample_category)) %>% 
  arrange(disease, donor, sample_category)

all_fcs <- read_all_fcs_files(sample_meta = ced_samples_meta,
  only_functional_AB = TRUE, fcs_trans = list(func = identity, cofactor = 1))

### DEBUG - some have AI-panel, one is missing CD28
bad_samples <- all_fcs$df$fcs %>% map(colnames) %>% 
  set_names(all_fcs$df$unique_name) %>% 
  keep(~'TIGIT' %in% .x) # AI panel
bad_samples
assert_that(length(bad_samples) == 0)
all_fcs$df$fcs %>% map(colnames) %>% 
  set_names(all_fcs$df$unique_name) %>% 
  purrr::discard(~ 'CD28' %in% .x) %>% 
  identity() ->
  bad_samples # missing CD28
bad_samples
assert_that(length(bad_samples) == 0)
###

ced_fcs <- all_fcs$df
fs <- flowset_of_subgroup(ced_fcs)
sample_ids <- rep(ced_fcs$unique_name, fsApply(fs, nrow))

dr <- data.frame(
  sample = sample_ids,
  fsApply(fs, flowCore::exprs),
  stringsAsFactors = FALSE
  ) %>% 
  as_data_frame() %>% 
  inner_join(
    select(ced_fcs, sample=unique_name, disease, disease_status, donor, sample_category, biosource)
  )

# missing_markers <- asbjorn_markers %>% 
#   anti_join(tibble(antibody= colnames(dr)))
# assert_that(nrow(missing_markers) == 0)
# dr %>% 
dr %>% 
  filter(disease_status != 'GFD') %>% #Drop GFD for now
  rownames_to_column('cell_id') %>% 
  # mutate(disease_category = case_when(disease == 'HC' ~ 'Healthy Control',
  #                                     disease_status == 'GFD' ~ "Celiac - GFD",
  #                                     disease_status == 'UCD' ~ "Celiac - Untreated"),
  #        disease_category = factor(disease_category, c('Celiac - Untreated', 'Celiac - GFD', 'Healthy Control'))) %>% 
  gather(variable, value, -cell_id, -disease, -disease_status, -sample_category, -sample, -donor, -biosource) %>% 
  identity() ->
  tdr

tdr %>% 
  group_by(biosource, disease, disease_status, sample_category, variable, donor) %>% 
  summarize(value = mean(value)) %>% 
  summarize(value = mean(value)) %>% 
  ungroup() %>% 
  # spread(sample_category, value) %>% 
  # mutate(log2_ratio = log2(`tetramer+` / `tetramer-`)) %>% 
  mutate(variable = str_replace(variable, '_', '-')) %>% 
  identity() ->
  tdr_gmean
  
tdr_gmean %>% 
  distinct(variable) %>% 
  anti_join(asbjorn_markers, by = c('variable' = 'antibody')) %>% 
  print(n = Inf)
missing_markers <- asbjorn_markers %>% 
  anti_join(distinct(tdr_gmean, variable), by = c('antibody' = 'variable' ))
assert_that(nrow(missing_markers) == 0)

tdr_gmean <- tdr_gmean %>% 
  inner_join(asbjorn_markers, by = c('variable' = 'antibody')) %>% 
  mutate(variable = factor(variable, rev(asbjorn_markers$antibody)))

##
## Save donor means to file for p-value computation elsewhere
##
tdr %>% 
  group_by(biosource, disease, disease_status, sample_category, variable, donor) %>% 
  summarize(value = mean(value)) %>% 
  ungroup() %>% 
  identity() ->
  donor_mean

tdr %>% 
  mutate(value = asinh(value) / 5.0) %>% 
  group_by(biosource, disease, disease_status, sample_category, variable, donor) %>% 
  summarize(value = mean(value)) %>% 
  ungroup() %>% 
  identity() ->
  donor_asinh_mean
  
donor_asinh_mean %>% 
  write_csv('out/cytof_marker_asinh_mean_per_donor.csv')

###
### Actual plotting
###

##
## get absolute expression
##

tdr_gmean %>% 
  filter(disease_status == 'UCD' & sample_category == 'tetramer+') %>% 
  filter(biosource == 'PBMC') %>%
  ggplot(aes(biosource, variable)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = 'grey20', high = 'yellow2', limits = c(0, 75), oob = squish) %>% 
  identity() ->
  p_blood_expr
p_blood_expr

tdr_gmean %>% 
  filter(disease_status == 'UCD' & sample_category == 'tetramer+') %>% 
  filter(biosource == 'SCS') %>%
  ggplot(aes(biosource, variable)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = 'grey20', high = 'yellow2', limits = c(0, 75), oob = squish) %>% 
  identity() ->
  p_gut_expr
p_gut_expr

tdr_gmean %>% 
  filter((disease_status == 'UCD' & sample_category == 'tetramer+') |
           (disease == 'HC' & sample_category == 'tetramer-')) %>% 
  select(biosource, sample_category, variable, value) %>% 
  spread(sample_category, value) %>% 
  mutate(fold_change = log2(`tetramer+` / `tetramer-`)) %>% 
  select(biosource, variable, fold_change) %>% 
  mutate(category = 'log(CeD+ / HC)') %>% 
  identity() ->
  fc_tetpos_over_HC

tdr_gmean %>% 
  filter(disease == 'Ced') %>% 
  spread(sample_category, value) %>% 
  mutate(fold_change = log2(`tetramer+` / `tetramer-`)) %>% 
  select(biosource, variable, fold_change) %>% 
  mutate(category = 'log(CeD+ / CeD-)') %>% 
  identity() ->
  fc_tetpos_over_Ced

bind_rows(fc_tetpos_over_Ced, fc_tetpos_over_HC) %>% 
  mutate(variable = factor(variable, rev(asbjorn_markers$antibody))) %>% 
  identity() ->
  fc_df
  

fc_df %>% 
  filter(biosource == 'PBMC') %>% 
  ggplot(aes(category, variable)) +
  geom_tile(aes(fill = fold_change)) +
  scale_fill_gradient2(low = 'dodgerblue4', mid = 'grey20', high = 'yellow2', limits = c(-4, 4), oob=squish) %>% 
  identity() ->
  p_fc_PBMC
p_fc_PBMC

fc_df %>% 
  filter(biosource == 'SCS') %>% 
  ggplot(aes(category, variable)) +
  geom_tile(aes(fill = fold_change)) +
  scale_fill_gradient2(low = 'dodgerblue4', mid = 'grey20', high = 'yellow2', limits = c(-3, 3), oob=squish) %>% 
  identity() ->
  p_fc_gut
p_fc_gut  

rnadf <- read_csv('input_data/renamed/cytof_log_ratio.csv') # TODO copy to this project
rnadf %>% 
  # ENSG00000227993 This HLA-DR gene has most reads -- but is not on a normal chromosome
  # To use a HLA-DR gene on primary assembly I would have to realign with a normal mapper to a normal
  # reference genome
  filter(common_name != 'HLA-DR' | gene_id == 'ENSG00000227993') %>%  # This HLA-DR gene has most reads -- use this
  filter(common_name != 'PD-1' | gene_id == 'ENSG00000188389') %>% # Either PD-1 gene has identical scores - this is the canonical one
  filter(common_name != 'CD45RA' | gene_id == 'ENSG00000081237') %>% # This is the canonical one and has the most reads
  select(-common_name) %>% 
  inner_join(select(asbjorn_markers_all, common_name, gene_name)) %>% 
  mutate(common_name = factor(common_name, levels = rev(asbjorn_markers_all$common_name)),
         comparison_group = if_else(str_detect(comparison_group, 'Control'), 'log(Ced+ / HC)', 'log(Ced+ / CeD-)'),
         comparison_group = fct_relevel(comparison_group, 'log(Ced+ / CeD-)') ) %>% 
  replace_na(replace = list(padj = 1)) %>% 
  identity() ->
  rna_seq_fc

rna_seq_fc %>% 
  filter(category == 'RNA-seq') %>% 
  ggplot(aes(comparison_group, common_name)) +
  geom_tile(aes(fill = log2FoldChange)) +
  geom_point(data=filter(rna_seq_fc, category == 'RNA-seq' & padj >= 0.05), shape = 4, size = 4, color = 'black') +
  scale_fill_gradient2(low = 'dodgerblue4', mid = 'grey20', high = 'yellow2', limits = c(-5, 5), oob = squish) +
  labs(y = '', x = '', title = '', fill = '') %>% 
  identity() ->
  p_rna_only
p_rna_only

##
## Per mail med asbjorn 8 desember. Split i to figurer. Fig1 med GARP og CCL22. Fig2 med CD84, CD200, CXCL13
##
# Fig A
rna_seq_fc %>% 
  filter(category == 'RNA-seq') %>% 
  filter(common_name %in% c('GARP', 'CCL22')) %>% 
  ggplot(aes(comparison_group, common_name)) +
  geom_tile(aes(fill = log2FoldChange)) +
  # geom_point(data=filter(rna_seq_fc, category == 'RNA-seq' & padj >= 0.05), shape = 4, size = 4, color = 'black') +
  scale_fill_gradient2(low = 'dodgerblue4', mid = 'grey20', high = 'yellow2', limits = c(-5, 5), oob = squish) +
  labs(y = '', x = '', title = '', fill = '') %>% 
  identity() ->
  p_rna_only_figA
p_rna_only_figA
# Fig B
rna_seq_fc %>% 
  filter(category == 'RNA-seq') %>% 
  filter(!(common_name %in% c('GARP', 'CCL22'))) %>% 
  ggplot(aes(comparison_group, common_name)) +
  geom_tile(aes(fill = log2FoldChange)) +
  # geom_point(data=filter(rna_seq_fc, category == 'RNA-seq' & padj >= 0.05), shape = 4, size = 4, color = 'black') +
  scale_fill_gradient2(low = 'dodgerblue4', mid = 'grey20', high = 'yellow2', limits = c(-5, 5), oob = squish) +
  labs(y = '', x = '', title = '', fill = '') %>% 
  identity() ->
  p_rna_only_figB
p_rna_only_figB


rna_seq_fc%>% 
  filter(category == 'CyTOF') %>% 
  ggplot(aes(comparison_group, common_name)) +
  geom_tile(aes(fill = log2FoldChange)) +
  scale_fill_gradient2(low = 'dodgerblue4', mid = 'grey20', high = 'yellow2', limits = c(-3, 3), oob = squish) %>% 
  identity() ->
  p_rna_cytof
p_rna_cytof 

list(p_blood_expr, p_fc_PBMC) %>% 
  map(~ . + theme(axis.title = element_blank(), legend.title = element_blank())) %>% 
  cowplot::plot_grid(plotlist = ., nrow = 1, rel_widths = c(1, 1.7)) %>% 
  identity() ->
  p_PBMC
list(p_gut_expr, p_fc_gut, p_rna_cytof) %>% 
  map(~ . + theme(axis.title = element_blank(), legend.title = element_blank())) %>% 
  cowplot::plot_grid(plotlist = ., nrow = 1, rel_widths = c(1, 1.7, 1.7)) %>% 
  identity() ->
  p_gut

dir.create('out/paper_figs/main', recursive = TRUE, showWarnings = FALSE)
ggsave('out/paper_figs/main/fig_rna-seq_only.svg', p_rna_only, width = 2.5, height = 4)
ggsave('out/paper_figs/main/fig_rna-seq_only_figA.svg', p_rna_only_figA, width = 3, height = 3)
ggsave('out/paper_figs/main/fig_rna-seq_only_figB.svg', p_rna_only_figB, width = 3, height = 3)
ggsave('out/paper_figs/main/cytof_expr_and_fc_PBMC.png', plot = p_PBMC, width = 6)
ggsave('out/paper_figs/main/cytof_expr_and_fc_PBMC.svg', plot = p_PBMC, width = 6)
ggsave('out/paper_figs/main/cytof_expr_and_fc_GUT.png', plot = p_gut, width = 8)
ggsave('out/paper_figs/main/cytof_expr_and_fc_GUT.svg', plot = p_gut, width = 8)


theme_set(old_theme)


###
### Make note of which samples and donors that is used
###
sample_summary_count <- ced_samples_meta %>% 
  distinct(donor, biosource, disease, disease_status) %>% 
  count(biosource, disease, disease_status) %>% 
  arrange(biosource, disease, n)

fp <- file('out/paper_figs/main/summary_stats_and_donors_used.txt', open = 'w')
write(glue('{date()}\n'), file = fp)
write(sample_summary_count %>% format_csv(), file = fp)
write(ced_samples_meta %>% format_csv(), file = fp)
close(fp)

## Suppl plot to reviewer for illustrating correspondance between RNA-seq and CyTOF
list(RNA_seq = rna_seq_fc %>% 
       filter(comparison_group == 'log(Ced+ / CeD-)') %>% 
       dplyr::select(variable = common_name, fold_change = log2FoldChange),
     CyTOF = fc_df %>% 
       filter(biosource == 'SCS' & category == 'log(CeD+ / CeD-)') %>% 
       dplyr::select(variable, fold_change)) %>% 
  bind_rows(.id = 'instrument') %>% 
  spread(instrument, fold_change) %>% 
  na.omit() %>% 
  identity() ->
  xdf

rsq_ced <- lm(CyTOF ~ RNA_seq, data = xdf) %>% 
  summary() %>% 
  .$r.squared %>% 
  formatC(digits = 2, format = 'f')

pacman::p_load(
  ggrepel
)
theme_set(theme_classic())
xdf %>% 
  ggplot(aes(CyTOF, RNA_seq, label = variable)) +
  geom_point(aes(color = variable)) +
  geom_smooth(method = 'lm') + 
  geom_text_repel() +
  labs(title = 'CyTOF vs RNA-seq log2 fold change for tetramer+ vs tetramer- in gut',
       caption = glue('R^2 = {rsq_ced}')) +
  theme(legend.position = 'none') %>% 
  identity() ->
  p_vs_tetneg
p_vs_tetneg  
ggsave('out/paper_figs/supplemental/cytof_vs_rna_tetp_vs_tetneg.pdf', p_vs_tetneg)

list(RNA_seq = rna_seq_fc %>% 
       filter(comparison_group == 'log(Ced+ / HC)') %>% 
       dplyr::select(variable = common_name, fold_change = log2FoldChange),
     CyTOF = fc_df %>% 
       filter(biosource == 'SCS' & category == 'log(CeD+ / HC)') %>% 
       dplyr::select(variable, fold_change)) %>% 
  bind_rows(.id = 'instrument') %>% 
  spread(instrument, fold_change) %>% 
  na.omit() %>% 
  identity() ->
  xdf_vs_hc

rsq_ced_vs_hc<- lm(CyTOF ~ RNA_seq, data = xdf_vs_tetneg) %>% 
  summary() %>% 
  .$r.squared %>% 
  formatC(digits = 2, format = 'f')

xdf_vs_hc%>% 
  ggplot(aes(CyTOF, RNA_seq, label = variable)) +
  geom_point(aes(color = variable)) +
  geom_smooth(method = 'lm') + 
  geom_text_repel() +
  labs(title = 'CyTOF vs RNA-seq log2 fold change for tetramer+ vs healthy controls in gut',
       caption = glue('R^2 = {rsq_ced_vs_hc}')) +
  theme(legend.position="none") %>% 
  identity() ->
  p_vs_hc
p_vs_hc
ggsave('out/paper_figs/supplemental/cytof_vs_rna_tetp_vs_hc.pdf', p_vs_hc)
  
  wfc_df %>% 
       filter(biosource == 'SCS' & category == 'log(CeD+ / CeD-)')