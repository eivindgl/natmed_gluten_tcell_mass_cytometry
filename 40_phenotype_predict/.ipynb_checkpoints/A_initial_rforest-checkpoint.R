pacman::p_load(
  randomForest, 
  caret, 
  tidyverse,
  glue
)
# set.seed(1234)
# ind = sample(2, nrow(iris), replace=TRUE, prob=c(0.7,0.3))
# trainData = iris[ind==1,]
# testData = iris[ind==2,]
# 
# iris_rf = randomForest(Species~., data=trainData, ntree=100, proximity=T)
# table(predict(iris_rf), trainData$Species)


tetdf <- read_csv('out/predict_tetp/pbmc_ucd_tetp_and_tetn.csv')
tetdf %>% 
  mutate(sample_category = fct_relevel(sample_category, 'tetramer-')) %>% 
  identity() ->
  tetdf
set.seed(1234)
ind = sample(2, nrow(tetdf), replace=TRUE, prob=c(0.7,0.3))
trainData = select(tetdf, -unique_name)[ind==1,]
testData = select(tetdf, -unique_name)[ind==2,]

tet_rf = randomForest(sample_category~., data=trainData, proximity=T)
table(predict(iris_rf), trainData$Species)
tetPred = predict(tet_rf, newdata=testData)
table(tetPred, testData$sample_category)

CM <-  table(tetPred, testData$sample_category)
accuracy <- sum(diag(CM)) / sum(CM)

alldf <- read_csv('out/predict_tetp/all.csv.gz')

allPred <- predict(tet_rf, newdata = alldf)

rdf <- alldf %>% 
  select(unique_name, sample_category) %>% 
  mutate(predgroup = allPred)

rdf %>% 
  group_by(unique_name, sample_category) %>% 
  summarize(nAT = sum(predgroup == 'tetramer+'),
            nOther = sum(predgroup == 'tetramer-'),
            nTot = nAT + nOther,
            per_million = nAT * 1e6 / n()) %>% 
  mutate(prct = 100.0 * nAT / (nAT + nOther),
         prct = glue('{format(prct, digits = 3)}')) %>% 
  identity() ->
  sdf

all_samples <- read_csv('input_data/renamed/sample_meta_full.csv') %>% 
  dplyr::select(-path, -note, -date, -instrument)

p <- all_samples %>%
  filter(sample_category %in% c('full', 'tetramer-')) %>% 
  filter(biosource == 'PBMC') %>% 
  filter(!(biosource == 'PBMC' & sample_category == 'tetramer-')) %>% 
  left_join(sdf) %>% 
  # filter(nOther > 1000) %>% 
  # mutate(per_million = per_million + 1) %>% 
  mutate(sample_group = if_else(disease == 'Flu', str_c("Flu_", disease_status), sample_group)) %>% 
  arrange(nOther) %>% 
  ggplot(aes(sample_group, per_million, fill  = sample_category)) +
  geom_jitter(aes(color = sample_group, shape = sample_category, label = donor), width = 0.3, size = 3) +
  theme(axis.text = element_text(angle = 45, vjust = 3)) +
  scale_y_continuous(limits = c( 50, 4000), oob = scales::squish)
p

pacman::p_load(
  plotly,
  ggrepel
)
ggplotly(p)

alldf %>% 
  filter(sample)

marker_sigdf <- importance(tet_rf) %>% 
  as.data.frame() %>% 
  rownames_to_column('marker') %>% 
  as_data_frame() %>% 
  arrange(desc(MeanDecreaseGini))


sdf %>% 
  # filter(sample_category == 'full') %>%
  inner_join(all_samples) %>% 
  filter(str_detect(donor, 'CD1414'))
  ggplot(aes(nTot, per_million, label = donor)) +
  geom_point(aes(color = disease)) +
  scale_y_log10() +
  scale_x_log10() +
  geom_text_repel( )
