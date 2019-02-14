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
mod_ced <- readRDS(file = 'out/predict_tetp/model_rf_tetp_gfd_ucd.robj')
mod_ai <- readRDS(file = 'out/predict_tetp/model_rf_tetp_ucd_no_CD62L_or_CCR7.robj')

ch_samples <- all_samples %>% 
  filter(disease_status == 'Challenge' & !str_detect(sample_category, 'tetramer'))

some_pbmc_CeD_samples <- all_samples %>% 
  filter(biosource == 'PBMC' & disease_status == 'UCD' & !str_detect(sample_category, 'tetramer')) %>% 
  sample_n(5)
  
pbmc_tetp_samples <- all_samples %>% 
  filter(sample_category == 'tetramer+')

samples <- bind_rows(
  ch_samples,
  some_pbmc_CeD_samples,
  pbmc_tetp_samples
)
  

bind_rows(
  # predict_set(model_name = 'autoimmune', model = mod_ai, samples = samples),
  predict_set(model_name = 'ced', model = mod_ced, samples = samples)
) %>% 
  identity() ->
  pred_results

pred_results %>% 
  spread(predicted_phenotype, num_cells) %>% 
  mutate(num_tot = tetramer_neg + tetramer_pos,
         per_million = tetramer_pos * 1e6 / num_tot) %>% 
  inner_join(all_samples) %>% 
  filter(disease == 'Ced') %>% 
  filter(sample_category == 'tetramer+') %>% 
  arrange(desc(per_million)) %>% 
  select(-date, -sample, -biosource) %>% 
  print(n = Inf)
  