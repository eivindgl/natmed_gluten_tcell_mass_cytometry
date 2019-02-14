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
                    "ICOS", "KLRG1", "PD-1", "CD45RA", "CCR7", 'CCR5', 'CD137', 'CD244', 'CXCR4', 'OX40')) %>% 
  mutate(fcs_name = str_replace_all(antibody, '[- ]', '_')) %>% 
  identity() ->
  predmarkers
tibble(fcs_name = colnames(full_train_df)) %>%
  full_join(predmarkers) %>%
  arrange(antibody) %>%
  print(n = Inf)

full_train_df %>% 
  map(colnames) %>% 
  unlist() %>% 
  tibble(fcs_name = .) %>% 
  dplyr::count(fcs_name) %>% 
  arrange(n) 
full_train_df %>% 
  set_names(ced_tet_samples$unique_name) %>% 
  map_dfr(function(x) length(colnames(x)))
  unlist() %>% 
  

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
  filter(biosource == 'SCS' & disease_status == 'UCD') %>% 
  filter(disease != 'PotentialCed' & donor != 'CD1414')
ced_tet_samples %>% 
  arrange(donor, sample_category)

full_train_df <- ced_tet_samples %>%
  mutate(h5path = glue('samples/{unique_name}')) %>%
  .$h5path %>%
  map(partial(read_dataset, 'out/cytof.h5')) %>% 
  map(~ as_tibble(.x[, predmarkers$fcs_name])) %>%
  # map(as_tibble) %>% 
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
  select(-unique_name, -disease_status)
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
saveRDS(m_rf, file = 'out/predict_tetp/model_rf_gut_tetp_ucd.robj')

m_rf$results %>% 
  mutate(mtry = factor(mtry)) %>% 
  ggplot(aes(logLoss, AUC, color = mtry, label = mtry)) +
  geom_point() + 
  geom_text_repel()

m_rf$results

m_rf$finalModel

randomForest::importance(m_rf$finalModel) %>% 
  as.data.frame() %>% 
  rownames_to_column('antibody') %>% 
  arrange(desc(MeanDecreaseGini)) %>% 
  mutate(antibody = fct_reorder(antibody, MeanDecreaseGini)) %>% 
  ggplot(aes(antibody, MeanDecreaseGini)) +
  geom_point(shape = 5, size = 4) +
  coord_flip() +
  theme_minimal()
