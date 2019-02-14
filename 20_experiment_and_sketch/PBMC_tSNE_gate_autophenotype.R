pacman::p_load(
  tidyverse
)
devtools::install_github("jkrijthe/Rtsne", ref = 'openmp')
library(Rtsne)
select <- dplyr::select
filter <- dplyr::filter
count <- dplyr::count

source('20_experiment_and_sketch/load_fcs_and_rename_according_to_panel.R')

samples_meta <- read_csv('input_data/renamed/sample_meta_full.csv') %>% 
  filter(biosource == 'PBMC' & !(sample_category %in% c('AutoPhenotype', 'tetramer-'))) %>% 
  filter(donor != 'BC11' | str_detect(path, '170419 UCDpat 1 PBM')) %>% 
  arrange(disease, donor, sample_category)

fcs <- read_all_fcs_files(sample_meta = samples_meta)
fs <- flowset_of_subgroup(fcs$df)

sample_ids <- rep(fcs$df$unique_name, fsApply(fs, nrow))
expr <- fsApply(fs, exprs)

MAX_CELLS_PER_SAMPLE <- 40000
## Find and skip duplicates
# dups <- which(!duplicated(expr))
## Data subsampling: create indices by sample
inds <- split(1:length(sample_ids), sample_ids)
## How many cells to downsample per-sample
tsne_ncells <- pmin(table(sample_ids), MAX_CELLS_PER_SAMPLE)
## Get subsampled indices
set.seed(1234)
tsne_inds <- lapply(names(inds), function(i){
  sample(inds[[i]], tsne_ncells[i], replace = FALSE)
})
tsne_inds <- unlist(tsne_inds)
tsne_expr <- expr[tsne_inds, ]

set.seed(12345)
tsne_out <- Rtsne(tsne_expr, check_duplicates = FALSE, pca = FALSE, num_threads = 6)
# Plot t-SNE colored by CD4 expression
dr <- data.frame(
  sample = sample_ids[tsne_inds],
  tSNE1 = tsne_out$Y[, 1], tSNE2 = tsne_out$Y[, 2],
  expr[tsne_inds, ], stringsAsFactors = FALSE
  ) %>% 
  as_data_frame() %>% 
  inner_join(
    select(fcs$df, sample=unique_name, disease, disease_status, donor, sample_category)
  )
dr %>% 
  write_csv('out/tmp/PBMC_tsne_data_40K.csv')
stop('all done')
p <- dr %>% 
  filter(sample_category == 'tetramer+') %>% 
  # filter(Disease != 'HC' & Disease != 'Ced') %>% 
  # filter(Disease != 'Ced') %>% 
  ggplot(aes(x = tSNE1, y = tSNE2, color = interaction(disease, disease_status))) +
  geom_point(aes(alpha = 0.3))
  stat_summary_2d(aes(z = frac_of_max), fun = median, bins = 50)+
  facet_wrap(~ variable, ncol = 4) +
  theme_minimal() +
  theme_minimal() +
  facet_wrap(~ interaction(Disease, `Disease Status`), ncol = 1) +
  theme(legend.position = 'none')

dr %>% 
  filter(sample_category == 'tetramer+' & disease_status == 'UCD') %>% 
  select(starts_with('tSNE')) %>% 
  summarize(x = median(tSNE1), y = median(tSNE2), x_sd = sd(tSNE1), y_sd = sd(tSNE2)) %>% 
  mutate(sd_f = 0.4,
         xmin = x - sd_f * x_sd,
         xmax = x + sd_f * x_sd,
         ymin = y - sd_f * y_sd,
         ymax = y + sd_f * y_sd) %>% 
  identity() ->
  mid_point
p + annotate('rect', xmin = mid_point$xmin, ymin = mid_point$ymin, xmax = mid_point$xmax, ymax = mid_point$ymax, size = 5, alpha = 0.7)

within_rectangle <- function(xs, ys, rect) {
  within_x <- rect$xmin < xs & xs < rect$xmax
  within_y <- rect$ymin < ys & ys < rect$ymax
  within_x & within_y
}

dr %>% 
  filter(sample_category == 'full') %>% 
  select(donor, disease, disease_status, starts_with('tSNE')) %>% 
  mutate(AutoPheno = within_rectangle(tSNE1, tSNE2, mid_point),
         phenotype = if_else(AutoPheno, 'AutoPheno', 'conventional')) %>% 
  count(donor, disease, disease_status, phenotype) %>% 
  spread(phenotype, n) %>% 
  mutate(log2_ratio = log2(AutoPheno) - log2(AutoPheno + conventional),
         cat = str_c(disease, disease_status, sep = '_')) %>% 
  ggplot() +
  geom_point(aes(cat, log2_ratio, color = donor))
