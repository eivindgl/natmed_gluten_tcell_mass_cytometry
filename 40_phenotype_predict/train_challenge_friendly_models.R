pacman::p_load(
  tidyverse,
  rhdf5,
  caret,
  doMC,
  glue,
  MLmetrics,
  e1071,
  tictoc,
  ggrepel
)

# predmarkers <- 
tibble(antibody = c("CCR4", "CCR6", "CD127", "CD161", "CD25", "CD27", "CD28", "CD38", 
                    "CD39", "CD57", "CD69", "CD73", "CTLA-4", "CXCR3", "CXCR5", "HLA-DR", 
                    "ICOS", "KLRG1", "PD-1", "CD45RA", "CD62L", "CCR7")) %>% 
  mutate(fcs_name = str_replace_all(antibody, '[- ]', '_')) %>% 
  identity() ->
  predmarkers

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
  filter(biosource == 'PBMC' & disease == 'Ced' & disease_status != 'Challenge') %>%
  filter(sample_category != 'full' & (is.na(note) | note != 'exclude_ac'))

full_train_df <- ced_tet_samples %>%
  mutate(h5path = glue('samples/{unique_name}')) %>%
  .$h5path %>%
  map(partial(read_dataset, 'out/cytof.h5')) %>%
  map(~ as_tibble(.x[, predmarkers$fcs_name])) %>% 
  set_names(ced_tet_samples$unique_name) %>% 
  bind_rows(.id = 'unique_name')

ced_tet_samples %>% 
  distinct(unique_name, sample_category, disease_status) %>% 
  mutate(tetramer = if_else(sample_category == 'tetramer+', 'tetramer_pos', 'tetramer_neg')) %>% 
  dplyr::select(-sample_category) %>% 
  inner_join(full_train_df) %>% 
  mutate(tetramer = fct_relevel(tetramer, 'tetramer_neg')) %>% 
  identity() ->
  full_train_df

##
## prepare training set
##
registerDoMC(cores = 6)
SEED <- 1234
set.seed(SEED)
control <- trainControl(method="repeatedcv", number=10, repeats=3, classProbs= TRUE, summaryFunction = multiClassSummary, search = 'grid')
metric <- "logLoss"

##
## Train UCD without CCR5 and CD39
##
train_df <- full_train_df %>% 
  filter(disease_status == 'UCD') %>% 
  select(-unique_name, -disease_status, -CD39)
nparam <- ncol(train_df) - 1
mtry <- sqrt(nparam)
tunegrid <- expand.grid(.mtry = c(1:floor(mtry), mtry, ceiling(mtry):(nparam / 2)))
tic('Training random forest')
set.seed(SEED)

m_rf <- train(tetramer ~ ., data=train_df, method="rf", metric=metric, 
              tuneGrid = tunegrid, trControl=control) #, preProcess = c("center", "scale") )
# no scaling for random forest - https://stackoverflow.com/questions/8961586/do-i-need-to-normalize-or-scale-data-for-randomforest-r-package
toc()
print(m_rf)
saveRDS(m_rf, file = 'out/predict_tetp/model_rf_tetp_ucd_no_CD39.robj')

m_rf$results %>% 
  mutate(mtry = factor(mtry)) %>% 
  ggplot(aes(logLoss, AUC, color = mtry, label = mtry)) +
  geom_point() + 
  geom_text_repel()
