pacman::p_load(
  tidyverse,
  caret,
  rhdf5,
  assertthat,
  ggthemes,
  glue,
  randomForest
)
all_samples <- read_csv('input_data/renamed/sample_meta_full.csv')

read_dataset <- function(path, h5_path) {
  M <- h5read(path, h5_path)
  colnames(M) <- h5readAttributes(path, h5_path)$colnames
  M
}

predict_sample <- function(unique_name, model_name, model) {
  print(glue('predicting [{model_name}] for {unique_name}'))
  h5_path <- glue("samples/{unique_name}")
  x <- read_dataset('out/cytof.h5', h5_path)
  missing_markers <- base::setdiff(model$coefnames, colnames(x))
  if (nrow(x) == 0) {
    print(glue("Skipping sample {unique_name} - Data set is empty (0 cells)"))
    return (NULL)
  }
  if (length(missing_markers) > 1) {
    missing_markers <- paste(missing_markers, collapse = ", ")
    print(glue('Skipping sample {unique_name} - Missing markers: {missing_markers}'))
    return (NULL)
  }
  x <- x[, model$coefnames, drop = FALSE]
  predict(model, x) %>% 
    table() %>% 
    as_tibble() %>% 
    mutate(unique_name = unique_name,
           model_name = model_name) %>% 
    dplyr::select(unique_name, model_name, everything()) %>% 
    set_names(c('unique_name', 'model_name', 'predicted_phenotype', 'num_cells'))
}

predict_set <- function(model_name, model, samples) {
  predfunc <- partial(predict_sample, model_name = model_name, model = model)
  samples %>% 
    .$unique_name %>% 
    map(predfunc) %>% 
    discard(is.null) %>% 
    bind_rows()
}
mod_ced <- readRDS(file = 'out/predict_tetp/model_rf_tetp_ucd.robj')
mod_ai <- readRDS(file = 'out/predict_tetp/model_rf_tetp_ucd_no_CD62L_or_CCR7.robj')
mod_ch <- readRDS(file = 'out/predict_tetp/model_rf_tetp_ucd_no_CD39.robj')
mod_ch_flubg <- readRDS(file = 'out/predict_tetp/model_rf_tetp_ucd_no_CD39_with_cancerflu.robj')
mod_cm_ch_flubg <- readRDS(file = 'out/predict_tetp/model_rf_tetp_ucd_no_CD39_CD62L_CCR7_with_cancerflu.robj')

ch_samples <- all_samples %>% 
  filter(disease_status == 'Challenge')

hc_samples <- all_samples %>% 
  filter(biosource == 'PBMC' & sample_category == 'full' & disease == 'HC')

# samples <- bind_rows(
#   ch_samples,
#   hc_samples
# )
samples <- all_samples

bind_rows(
  # predict_set(model_name = 'autoimmune', model = mod_ai, samples = samples),
  predict_set(model_name = 'ced', model = mod_ced, samples = samples),
  predict_set(model_name = 'ch_friendly', model = mod_ch, samples = samples)
  # predict_set(model_name = 'ch_friendly_flubg', model = mod_ch_flubg, samples = samples),
  # predict_set(model_name = 'ch_cm_friendly_flubg', model = mod_cm_ch_flubg, samples = samples)
) %>% 
  identity() ->
  pred_results

pred_results %>% 
  inner_join(samples) %>% 
  filter(sample_category == 'full' & biosource == 'PBMC') %>% 
  spread(predicted_phenotype, num_cells) %>% 
  mutate(num_tot = tetramer_neg + tetramer_pos,
         per_million = tetramer_pos * 1e6 / num_tot) %>% 
  arrange(desc(per_million)) %>% 
  dplyr::select(disease, disease_status, model_name, donor, per_million) %>% 
  mutate(disease_status = tolower(disease_status),
         disease_status = case_when(disease_status == 'normal' & disease == 'HC' ~ 'HC', 
                                    disease_status == 'normal' ~ 'Other',
                                    TRUE ~ disease_status)) %>% 
  # filter(model_name == 'ch_friendly') %>% 
  ggplot(aes(disease_status, per_million)) + 
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 1.5e4), oob = scales::squish) +
# scale_y_log10() +
  facet_wrap(~ model_name)

pred_results %>% 
  inner_join(samples) %>% 
  filter(sample_category == 'full' & biosource == 'PBMC') %>% 
  spread(predicted_phenotype, num_cells) %>% 
  mutate(num_tot = tetramer_neg + tetramer_pos,
         per_million = tetramer_pos * 1e6 / num_tot) %>% 
  arrange(desc(per_million)) %>% 
  dplyr::select(disease, disease_status, model_name, donor, per_million) %>% 
  mutate(disease_status = tolower(disease_status),
         disease_status = case_when(disease_status == 'normal' & disease == 'HC' ~ 'HC', 
                                    disease_status == 'normal' ~ 'Other',
                                    TRUE ~ disease_status)) %>% 
  ggplot(aes(disease_status, per_million)) + 
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 5000), oob = scales::squish) +
# scale_y_log10() +
  facet_wrap(~ model_name, ncol = 1)
samples %>% 
  filter(biosource == 'PBMC') %>% 
  dplyr::select(-category, -biosource, -date, -sample, -note, -path) %>% 
  inner_join(pred_results) %>% 
  dplyr::select(-instrument) %>% 
  # filter(donor == 'CD1575') %>% 
  spread(predicted_phenotype, num_cells) %>% 
  mutate(tot_cells = tetramer_pos + tetramer_neg,
         per_million = tetramer_pos * 1e6 / tot_cells) %>% 
  arrange(desc(per_million)) %>% 
  identity() ->
  results
  
pred_df <- results %>% 
  filter(sample_category == 'full') %>% 
  filter(donor != '22-002') # Flu donor with cancer++

## Ensure we use paper relevant donors and groupings
##
paper_subsets <- list()
pred_df %>% 
  filter(sample_group == 'PBMC_AutoImmune') %>% 
  mutate(sample_subgroup = case_when(str_detect(disease, '^SSc') ~ 'SSc',
                                     disease == 'SLE' ~ 'SLE',
                                     TRUE ~ 'Other')) %>% 
  identity() ->
  paper_subsets$ai 
paper_subsets$ai %>% 
  distinct(donor, sample_subgroup) %>%
  count(sample_subgroup) # should be 10, 10, 5

pred_df %>%  #PBMC UCD OK
  filter(sample_group == 'PBMC_UCD') %>% 
  filter(sample_category == 'full') %>% 
  arrange(donor)

pred_df %>%  # PBMC GFD OK
  filter(sample_group == 'PBMC_GFD') %>% 
  filter(sample_category == 'full') %>% 
  arrange(donor)

pred_df %>%  # PBMC HC OK -- 10 with CeD panel - 17 with AI panel (BC11 overlaps)
  filter(sample_group == 'PBMC_HC') %>% 
  filter(sample_category == 'full') %>% 
  distinct(donor) %>% 
  arrange(donor) %>% 
  print(n = Inf)

pred_df %>%  # PBMC HC -- 10 with CeD panel - 17 with AI panel 
  filter(sample_group == 'PBMC_HC') %>% 
  filter(sample_category == 'full') %>% 
  add_count(donor) %>% 
  filter(n > 2)

pred_df %>%  # PBMC GFD OK
  filter(sample_category == 'full') %>% 
  filter(sample_group == 'PBMC_GFD' |
           sample_group == 'PBMC_UCD' | 
           sample_group == 'PBMC_HC') %>% 
  mutate(sample_subgroup = sample_group) %>% 
  bind_rows(paper_subsets$ai) %>% 
  identity() ->
  x

pred_df %>% 
  filter(unique_name == 'BC11_full_N1') %>% 
  mutate(sample_group = 'PBMC_HC',
         sample_subgroup = 'PBMC_HC') %>% 
  bind_rows(x)  %>% 
  identity() ->
  x

pred_df %>% 
  filter(disease_status == 'Challenge') %>% 
  mutate(sample_group = 'PBMC_CeD_Challenge',
         sample_subgroup = 'PBMC_CeD_Challenge') %>% 
  bind_rows(x)  %>% 
  identity() ->
  x

pred_df %>% 
  filter(sample_category == 'full') %>% 
  filter(sample_group == 'PBMC_Flu') %>% 
  mutate(sample_subgroup = disease_status) %>% 
  bind_rows(x) %>% 
  arrange(sample_group, sample_subgroup) %>% 
  identity() ->
  x
pacman::p_load(ggrepel)
x %>% 
  ggplot((aes(sample_subgroup, per_million, color = sample_group, label = donor))) +
  geom_boxplot(outlier.colour = NA) +
  geom_jitter(aes(color = sample_group), size = 5, width = 0.3) +
  geom_text_repel(data = filter(x, per_million > 6000)) +
  scale_y_continuous(limits = c(0, 10000), oob = scales::squish) +
  theme(axis.text.x = element_text(angle = 45, vjust = .5)) +
  theme_minimal() +
  facet_wrap(~ model_name, ncol = 1)
ggsave('out/paper_figs/main/prediction_figure/plot_sketch.png')

dir.create('out/paper_figs/main/prediction_figure')
x %>% 
  arrange(sample_group, sample_subgroup, per_million) %>% 
  filter(model_name == 'ced') %>% 
  write_csv('out/paper_figs/main/prediction_figure/predicted_tetp_per_million_with_CD39.csv')
x %>% 
  arrange(sample_group, sample_subgroup, per_million) %>% 
  filter(model_name == 'ch_friendly') %>% 
  write_csv('out/paper_figs/main/prediction_figure/predicted_tetp_per_million_without_CD39.csv')

randomForest::importance(mod_ch$finalModel) %>% 
  as.data.frame() %>% 
  rownames_to_column('antibody') %>% 
  arrange(desc(MeanDecreaseGini))
randomForest::importance(mod_ch$finalModel) %>% 
  as.data.frame() %>% 
  rownames_to_column('antibody') %>% 
  arrange(desc(MeanDecreaseGini)) %>% 
  write_csv('out/paper_figs/main/prediction_figure/model_without_CD39_parameter_importance.csv')


randomForest::importance(mod_ch$finalModel) %>% 
  as.data.frame() %>% 
  rownames_to_column('antibody') %>% 
  as_tibble() %>% 
  identity() ->
  mdg
mdg %>% 
  mutate(antibody = fct_reorder(antibody, MeanDecreaseGini)) %>% 
  ggplot(aes(antibody, MeanDecreaseGini)) +
  geom_point(shape = 5, size = 4) +
  coord_flip() +
  theme_minimal() +
  labs(x = 'Marker') +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))
ggsave('out/paper_figs/main/prediction_figure/model_without_CD39_parameter_importance.png')
ggsave('out/paper_figs/main/prediction_figure/model_without_CD39_parameter_importance.pdf')
