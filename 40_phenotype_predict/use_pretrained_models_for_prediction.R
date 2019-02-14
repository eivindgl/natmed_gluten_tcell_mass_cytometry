pacman::p_load(
  tidyverse,
  caret,
  rhdf5,
  assertthat,
  ggthemes
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

predict_set <- function(model_name, model) {
  predfunc <- partial(predict_sample, model_name = model_name, model = model)
  all_samples %>% 
    filter(biosource == 'PBMC' & !str_detect(sample_category, 'tetramer')) %>% 
    .$unique_name %>% 
    map(predfunc) %>% 
    discard(is.null) %>% 
    bind_rows()
}

bind_rows(
  predict_set(model_name = 'autoimmune', model = mod_ai),
  predict_set(model_name = 'ced', model = mod_ced)
) %>% 
  identity() ->
  pred_results

pred_results %>% 
  spread(predicted_phenotype, num_cells) %>% 
  mutate(num_tot = tetramer_neg + tetramer_pos,
         per_million = tetramer_pos * 1e6 / num_tot) %>% 
  inner_join(all_samples) %>% 
  select( -instrument, -date, -sample, -path) %>% 
  identity() ->
  pred_df

## Ensure we use paper relevant donors and groupings
##
paper_subsets <- list()
pred_df %>% 
  filter(sample_category == 'full') %>% 
  filter(sample_group == 'PBMC_AutoImmune') %>% 
  mutate(sample_subgroup = case_when(str_detect(disease, '^SSc') ~ 'SSc',
                                     disease == 'SLE' ~ 'SLE',
                                     TRUE ~ 'Other')) %>% 
  identity() ->
  paper_subsets$ai 
paper_subsets$ai %>% 
  distinct(donor, sample_subgroup) %>%
  count(sample_subgroup) # should be 10, 10, 15

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
  select(sample_group, sample_subgroup, unique_name, predicted_neg=tetramer_neg, predicted_pos=tetramer_pos, num_cells=num_tot, per_million) %>% 
  arrange(sample_group, sample_subgroup) %>% 
  identity() ->
  x

x %>% 
  ggplot((aes(sample_subgroup, per_million, color = sample_group))) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 1e4), oob = scales::squish) +
  theme(axis.text.x = element_text(angle = 45, vjust = .5))
ggsave('out/paper_figs/main/prediction_figure/plot_sketch.png')

xs <- split(x, x$model_name)
dir.create('out/paper_figs/main/prediction_figure')
xs$ced %>% 
  write_csv('out/paper_figs/main/prediction_figure/predicted_tetp_per_million_with_CD62L_CCR7.csv')
xs$autoimmune %>% 
  write_csv('out/paper_figs/main/prediction_figure/predicted_tetp_per_million_without_CD62L_CCR7.csv')

randomForest::importance(mod_ai$finalModel) %>% 
  as.data.frame() %>% 
  rownames_to_column('antibody') %>% 
  arrange(desc(MeanDecreaseGini)) %>% 
  write_csv('out/paper_figs/main/prediction_figure/model_without_CD62L_CCR7_parameter_importance.csv')

randomForest::importance(mod_ced$finalModel) %>% 
  as.data.frame() %>% 
  rownames_to_column('antibody') %>% 
  arrange(desc(MeanDecreaseGini)) %>% 
  write_csv('out/paper_figs/main/prediction_figure/model_with_CD62L_CCR7_parameter_importance.csv')


## DEBUG
##
# predict_set_bc11 <- function(model_name, model) {
#   predfunc <- partial(predict_sample, model_name = model_name, model = model)
#   all_samples %>% 
#     filter(biosource == 'PBMC' & !str_detect(sample_category, 'tetramer')) %>% 
#     .$unique_name %>% 
#     map(predfunc) %>% 
#     discard(is.null) %>% 
#     bind_rows()
# }
# 
# 
# bind_rows(
#   predict_set(model_name = 'autoimmune', model = mod_ai),
#   predict_set(model_name = 'ced', model = mod_ced)
# ) %>% 
#   identity() ->
#   pred_results
