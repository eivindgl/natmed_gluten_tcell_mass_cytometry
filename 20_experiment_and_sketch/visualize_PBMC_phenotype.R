pacman::p_load(
  tidyverse,
  devtools,
  viridis
)
devtools::install_github("jkrijthe/Rtsne", ref = 'openmp')
library(Rtsne)


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

all_fcs <- read_all_fcs_files(only_functional_AB = TRUE)
all_fcs$df <- all_fcs$df %>% 
  mutate(unique_name = str_c(sample, sample_category, sep = '_'))
  
common_ab <- all_fcs$df$fcs %>% 
  map(colnames) %>% 
  unlist() %>% 
  tibble(Antibody = .) %>% 
  count(Antibody) %>% 
  filter(n == nrow(all_fcs$df)) %>% 
  magrittr::extract2('Antibody')

hc_and_ap <- all_fcs$df %>%
  filter(biosource == 'PBMC' &
           donor != 'BC11') %>% 
  filter(Disease == 'HC' | !(sample_category == 'full' | sample_category == 'tetramer-'))
  
flowset_of_subgroup <- function(fcs_df) {
  fcs <- fcs_df %>%
    mutate(fcs = select_fcs_subset(fcs, common_ab))
  fcs <- setNames((fcs$fcs), fcs$unique_name)
  as(fcs, 'flowSet')
}
fsha <- flowset_of_subgroup(hc_and_ap)
sample_ids <- rep(hc_and_ap$unique_name, fsApply(fsha, nrow))
expr <- fsApply(fsha, exprs)

MAX_CELLS_PER_SAMPLE <- 3000
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
    select(hc_and_ap, sample=unique_name, Disease, `Disease Status`, donor, sample_category)
  )
dr %>% 
  write_csv('out/tmp/tsne_data.csv')


dr %>% 
  # filter(Disease != 'HC' & Disease != 'Ced') %>% 
  # filter(Disease != 'Ced') %>% 
ggplot(aes(x = tSNE1, y = tSNE2, color = interaction(Disease, `Disease Status`))) +
  geom_point(aes(alpha = 0.3)) +
  theme_minimal() +
  facet_wrap(~ interaction(Disease, `Disease Status`), ncol = 1) +
  theme(legend.position = 'none')
ggsave('out/tmp/HCbg_tetpos_and_phenotype.png', width =4, height = 12)

dr %>% 
ggplot(aes(x = tSNE1, y = tSNE2, color = CD39, shape = Disease)) +
  geom_point() 


dr %>% 
  rownames_to_column('cell_id') %>% 
  select(-sample, -donor) %>% 
  gather(variable, value, -cell_id, -tSNE1, -tSNE2, -Disease, -`Disease Status`, -sample_category) %>% 
  identity() ->
  tdr

ptdr <- tdr %>% 
  group_by(variable) %>% 
  mutate(frac_of_max = value / max(value)) %>% 
  ungroup()

p_fracsummary <- ptdr %>% 
  # filter(variable %in% c('CD39', 'CD57')) %>%
  ggplot(aes(tSNE1, tSNE2)) +
  stat_summary_2d(aes(z = frac_of_max), fun = median, bins = 50)+
  facet_wrap(~ variable, ncol = 4) +
  theme_minimal() +
  theme(legend.position = 'bottom') +
  labs(fill = 'AB scaled intensity',
       title = 'Marker intensity distribution over t-SNE plot',
       subtitle = 'Composed of all CD4 T-cells from healthy controls and tetramer+ or phenotype+ cells from patients') +
  scale_fill_viridis(option = 'B')
# p_fracsummary
ggsave('out/tmp/cd4_by_scaled_intensity.png', width = 12, height = 14, plot = p_fracsummary) 

p_exprsummary <- ptdr %>% 
  # filter(variable %in% c('CD39', 'CD57')) %>%
  ggplot(aes(tSNE1, tSNE2)) +
  stat_summary_2d(aes(z = value), fun = median, bins = 50)+
  facet_wrap(~ variable, ncol = 4) +
  theme_minimal() +
  theme(legend.position = 'bottom') +
  labs(fill = 'AB absolute intensity',
       title = 'Marker intensity distribution over t-SNE plot',
       subtitle = 'Composed of all CD4 T-cells from healthy controls and tetramer+ or phenotype+ cells from patients') +
  scale_fill_viridis(option = 'B')
p_exprsummary
ggsave('out/tmp/cd4_by_absolute_intensity.png', width = 12, height = 14, plot = p_exprsummary) 
