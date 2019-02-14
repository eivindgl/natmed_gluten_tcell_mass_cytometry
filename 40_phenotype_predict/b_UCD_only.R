pacman::p_load(
  tidyverse,
  rhdf5,
  caret,
  doMC,
  glue,
  kknn,
  MLmetrics,
  tictoc
)

predmarkers <- read_csv('input_data/common_markers_for_pred.csv') %>%
  mutate(fcsmarker = str_replace_all(Antibody, '[ -]', '_'))

read_dataset <- function(path, h5_path) {
  M <- h5read(path, h5_path)
  colnames(M) <- h5readAttributes(path, h5_path)$colnames
  M
}

##
## Extract training data
##
all_samples  <- read_csv('input_data/renamed/sample_meta_full.csv')
ced_tet_samples <- all_samples %>%
  filter(disease_status == 'UCD') %>% 
  filter(biosource == 'PBMC' & disease == 'Ced' & disease_status != 'Challenge') %>%
  filter(sample_category != 'full' & (is.na(note) | note != 'exclude_ac'))

full_train_df <- ced_tet_samples %>%
  mutate(h5path = glue('samples/{unique_name}')) %>%
  .$h5path %>%
  map(partial(read_dataset, 'out/cytof.h5')) %>%
  map(~ as_tibble(.x[, predmarkers$fcsmarker])) %>% 
  set_names(ced_tet_samples$unique_name) %>% 
  bind_rows(.id = 'unique_name')

ced_tet_samples %>% 
  distinct(unique_name, sample_category, disease_status) %>% 
  mutate(tetramer = if_else(sample_category == 'tetramer+', 'tetramer_pos', 'tetramer_neg')) %>% 
  dplyr::select(-sample_category) %>% 
  inner_join(full_train_df) %>% 
  identity() ->
  full_train_df

full_train_meta <- full_train_df %>% 
  select(unique_name, disease_status, tetramer) %>% 
  mutate(tetramer_status = str_c(disease_status, tetramer, sep = '__'))
full_train_df <- full_train_df %>% 
  select(-unique_name, -disease_status) %>% 
  mutate(tetramer = fct_relevel(tetramer, 'tetramer_neg'))

##
## prepare training set
##
registerDoMC(cores = 6)
SEED <- 1234

## Not needed since I to repeated k-fold cross validation
## see https://machinelearningmastery.com/how-to-estimate-model-accuracy-in-r-using-the-caret-package/
# set.seed(SEED)
# trainingRows <- createDataPartition(full_train_meta$tetramer_status, p = .8, list = FALSE)
# traindf <- full_train_df %>% 
#   slice(trainingRows)
# testdf <- full_train_df %>% 
#   slice(-trainingRows)

set.seed(SEED)
control <- trainControl(method="repeatedcv", number=10, repeats=3, classProbs= TRUE, summaryFunction = multiClassSummary)
# can be "Accuracy",   "logLoss", "ROC",   "Kappa"
metric <- "logLoss"

tic('Training random forest')
set.seed(SEED)
m_rf <- train(tetramer ~ ., data=full_train_df, method="rf", metric=metric, 
              trControl=control) #, preProcess = c("center", "scale") )
                                # no scaling for random forest - https://stackoverflow.com/questions/8961586/do-i-need-to-normalize-or-scale-data-for-randomforest-r-package
toc()
print(m_rf)
# saveRDS(m_rf, file = 'out/predict_tetp/model_rf_tetp_ucd_no_CCR7_or_CD62L.robj')

###
### Ensure that model is reasonable
###
predict_sample <- function(unique_name) {
  print(glue('predicting for {unique_name}'))
  h5_path <- glue("samples/{unique_name}")
  x <- read_dataset('out/cytof.h5', h5_path)
  missing_markers <- base::setdiff(predmarkers$fcsmarker, colnames(x))
  if (nrow(x) == 0) {
    print(glue("Skipping sample {unique_name} - Data set is empty (0 cells)"))
    return (NULL)
  }
  if (length(missing_markers) > 1) {
    missing_markers <- paste(missing_markers, collapse = ", ")
    print(glue('Skipping sample {unique_name} - Missing markers: {missing_markers}'))
    return (NULL)
  }
  x <- x[, predmarkers$fcsmarker, drop = FALSE]
  predict(m_rf, x) %>% 
    table() %>% 
    as_tibble() %>% 
    mutate(unique_name = unique_name) %>% 
    dplyr::select(unique_name, everything()) %>% 
    set_names(c('unique_name', 'tetramer', 'num_cells'))
}

all_samples %>% 
  anti_join(ced_tet_samples) %>% 
  filter(biosource == 'PBMC') %>% 
  .$unique_name %>% 
  map(predict_sample) %>% 
  discard(is.null) %>% 
  bind_rows() %>% 
  identity() ->
  pred_results_raw

pred_results_raw %>% 
  spread(tetramer, num_cells) %>% 
  mutate(num_tot = tetramer_neg + tetramer_pos,
         per_million = tetramer_pos * 1e6 / num_tot) %>% 
  inner_join(all_samples) %>% 
  select(-path, -instrument, -date, -sample) %>% 
  identity() ->
  pred_results
pred_results %>% 
  write_csv('out/predict_tetp/pred_CD62L_UCD_results.csv')

pred_results %>% 
  filter(sample_category == 'full') %>% 
  arrange(desc(per_million)) %>% 
  print(n = Inf)

pred_results %>% 
  filter(sample_category == 'full') %>% 
  filter(sample_group != 'ungrouped') %>% 
  mutate(disease_status = if_else(sample_group == 'PBMC_AutoImmune', 'autoimmune', disease_status),
         disease_status = tolower(disease_status),
         disease_status = fct_relevel(disease_status, 'normal', 'recovering', 'gfd'),
         sample_group = fct_relevel(sample_group, 'PBMC_HC', 'PBMC_Flu', 'PBMC_GFD')) %>% 
  ggplot(aes(sample_group, per_million, color = disease_status)) + 
  geom_boxplot() +
  # geom_jitter(width = 0.1, alpha = 0.5) +
  scale_y_continuous(limits = c(0, 6e3), oob = scales::squish) +
  theme_minimal()
ggsave('out/predict_tetp/predicted_AT_cells_UCD_no_CD62L_or_CCR7.png')
    # + 
  # coord_flip()

randomForest::importance(m_rf$finalModel) %>% 
  as.data.frame() %>% 
  rownames_to_column() %>%  
  arrange(MeanDecreaseGini)
