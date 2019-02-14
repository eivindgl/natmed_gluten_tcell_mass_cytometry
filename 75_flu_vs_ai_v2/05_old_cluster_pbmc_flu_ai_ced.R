pacman::p_load(
  flowCore,
  magrittr,
  tidyverse,
  stringr,
  matrixStats,
  glue,
  FlowSOM,
  ConsensusClusterPlus,
  RColorBrewer,
  pheatmap,
  rhdf5,
  assertthat,
  lme4
)
##
## Ensure we use the right dplyr functions
##
select <- dplyr::select
filter <- dplyr::filter
count <- dplyr::count
##
## Initial params
##
SEED <- 123456
cnum <- 30# number of clusters
cluster_dir <- glue('out/ai_vs_flu/cluster/seed_{SEED}')
cluster_path <- file.path(cluster_dir, glue('cluster_c{cnum}_seed{SEED}.csv.gz'))
dir.create(cluster_dir, showWarnings = FALSE, recursive = TRUE)
samples <- read_csv('out/ai_vs_flu/samples_used.csv')
panel <- tibble(antibody = c("CCR4", "CCR6", "CD127", "CD161", "CD25", "CD27", "CD28", "CD38", 
                                   "CD39", "CD57", "CD69", "CD73", "CTLA-4", "CXCR3", "CXCR5", "HLA-DR", 
                                   "ICOS", "KLRG1", "PD-1", "CD45RA", "CD62L", "CCR7")) %>% 
  mutate(fcs_name = str_replace_all(antibody, '[- ]', '_')) %>% 
  filter(!(antibody %in% c('CD45RA', 'CD73')))
panel_order <- read_csv('input_data/original/asbjorn_markers_specific.csv')
panel %>% 
  inner_join(panel_order, ., by = c('common_name' = 'antibody' )) %>% 
  identity() ->
  panel_ordered
adf <- read_csv('out/ai_vs_flu/activated_tcells.csv.gz')

##
## Initialize
##
x <- adf %>% 
  select(-sample_group, -sample_category, -unique_name)
fcs <- base::split(x, adf$unique_name) %>% 
  map(function(x) flowFrame(as.matrix(x))) %>% 
  flowSet()

if(!file.exists(cluster_path)) {
  fsom <- ReadInput(fcs, transform = FALSE, scale = FALSE)
  set.seed(SEED)
  som <- BuildSOM(fsom, colsToUse = panel$fcs_name, xdim = 15, ydim = 15, rlen = 15) 
                  # colsToUse = discard(panel$fcs_name, ~ . %in% c('CTLA_4'CD73')))
  
  codes <- som$map$codes
  plot_outdir <- file.path(cluster_dir, "consensus")
  
  mc <- ConsensusClusterPlus(t(codes), maxK = cnum, reps = 10000,
                             pItem = 0.9, pFeature = 0.9, title = plot_outdir, plot = "png",
                             clusterAlg = "hc", innerLinkage = "average", finalLinkage = "average",
                             distance = "euclidean", seed = SEED)
  
  ## Get cluster ids for each cell
  
  code_clustering1 <- mc[[cnum]]$consensusClass
  clusterdf <- tibble(cell_n = fsApply(fcs, flowCore::exprs)[, 'cell_n'],
                      cell_clustering = code_clustering1[som$map$mapping[,1]])
  
  clusterdf %>% 
    write_csv(cluster_path)
} else {
  clusterdf <- read_csv(cluster_path)
  cell_clustering1 <- clusterdf$cell_clustering
}
clusterdf %>% 
  mutate(cell_clustering = as.integer(cell_clustering),
         cluster = case_when(cell_clustering %in% c(1, 6,18,19,20) ~ 'CM',
                             # cell_clustering == 1 ~ 'CM CD27- CCR4-',
                             cell_clustering == 14 ~ 'AT',
                             cell_clustering == 5 ~ 'CD28+CXCR3+HLA-DR-',
                             # cell_clustering == 3 ~ 'Exhausted',
                             cell_clustering == 21 ~ 'Flu',
                             cell_clustering == 25 ~ 'Flu',
                             TRUE ~ 'undef'),
         cell_clustering = factor(cell_clustering)) %>% 
  identity() ->
  clusterdf

clusterdf %>% 
  inner_join(adf) %>% 
  identity() ->
  cadf

tsnedf <- read_csv('out/ai_vs_flu/tsne_tweak_df.csv')
clusterdf %>% 
  inner_join(tsnedf, by = 'cell_n') %>% 
  mutate(cell_clustering = factor(cell_clustering)) %>% 
  inner_join(samples) %>% 
  mutate(sample_group = if_else(sample_category == 'tetramer+', 'PBMC_UCD_tetp', sample_group)) %>% 
  mutate(disease = if_else(sample_category == 'tetramer+', 'CeD_tet+', disease)) %>% 
  identity() ->
  sdf

sdf %>% 
  ggplot(aes(tSNE1, tSNE2)) +
  geom_point(mapping = aes(color = sample_category), size = 0.6) +
  facet_wrap(~ cell_clustering, ncol = 3) +
  theme_minimal()

clusterdf

sdf %>% 
  ggplot(aes(tSNE1, tSNE2)) +
  geom_point(data = filter(sdf, cluster == 'undef'), color = 'gray', alpha = 0.4, size = 0.5) +
  geom_point(data = filter(sdf, cluster != 'undef'), mapping = aes(color = cluster), size = 0.6, alpha = 0.9) +
  # geom_density_2d() +
  facet_wrap(~ disease, ncol = 2) +
  theme_minimal() +
  scale_color_brewer(palette = 'Set1', type = 'qual')
  # scale_color_viridis_d(option = 'D')
ggsave('out/ai_vs_flu/tSNE_by_group_with_clusters.png')
sdf %>% 
  filter(cluster == 'undef') %>% 
  ggplot(aes(tSNE1, tSNE2)) +
  geom_point(aes(color = cell_clustering), alpha = 0.5, size = 1) +
  # geom_density_2d() +
  facet_wrap(~ disease, ncol = 2) +
  theme_minimal()

sdf %>% 
  # filter(cell_clustering == 14) %>% 
  filter(disease_status == 'infected') %>% 
  ggplot(aes(tSNE1, tSNE2)) +
  geom_point(alpha = 0.5, size = 1) +
  # geom_bin2d() +
  facet_wrap(~ cell_clustering, ncol = 4) +
  theme_minimal()
  
sdf %>% 
  filter(sample_category == 'tetramer+') %>%
  janitor::tabyl(cell_clustering) %>% 
  mutate(percent = percent * 100) %>% 
  arrange(n)

sdf %>% 
  arrange(desc(cell_n)) %>% 
  glimpse()

expr <- fsApply(fcs, flowCore::exprs)


tsnedf %>% 
  arrange(cell_n) %>% 
  tail() %>% 
  glimpse()

expr[2406330, ]
dim(expr)


clusterdf %>% 
  inner_join(adf) %>% 
  select(-cell_n, -unique_name, -sample_group, -sample_category) %>% 
  gather(marker, expression, -cell_clustering, -cluster) %>% 
  group_by(cell_clustering, cluster, marker) %>% 
  summarise(expression = median(expression)) %>% 
  ungroup() %>% 
  identity() ->
  cluster_expr
cluster_expr %>% 
  group_by(marker) %>% 
  summarize(top_expr = max(expression)) %>% 
  inner_join(cluster_expr) %>% 
  mutate(prct = expression / top_expr * 100) %>% 
  identity() ->
  cluster_expr_prct
cluster_expr_prct %>% 
  mutate(cluster_name = str_c(cell_clustering, cluster, sep = ' ')) %>% 
  ggplot(aes(marker, cluster_name , fill = prct)) +
  geom_tile() +
  viridis::scale_fill_viridis() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

cluster_expr_prct %>% 
  filter(marker == 'CD25') %>% 
  arrange(desc(prct))
  
sdf %>% 
  filter(as.integer(cell_clustering) %in% c(9)) %>%
  # filter(as.integer(cell_clustering) %in% c(7, 16, 11, 10)) %>%
  # filter(sample_group == 'PBMC_UCD') %>% 
  ggplot(aes(tSNE1, tSNE2, color = cell_clustering)) +
  geom_point(alpha = 0.5, size = 1) +
  # geom_density_2d() +
  facet_wrap(~ disease, ncol = 2) +
  theme_minimal()

##
## Differential expression testing
##
cadf %>% 
  inner_join(samples) %>% 
  filter(sample_category == 'full') %>% 
  count(disease, cluster) %>% 
  add_count(disease, wt = n) %>% 
  mutate(frac = n / nn) %>% 
  rename(cells_in_cluster = n, cells_in_sample = nn) %>% 
  identity() ->
  counts_table

cluster_colors <- c(CM = '#1b9e77', 'CD28+CXCR3+HLA-DR-' = '#1f78b4', 'CD161+CCR6+' = '#d95f02', 'AT' = '#e78ac3')
counts_table %>% 
  # filter(cluster != 'undef') %>%
  # filter(disease != 'Flu_recovering') %>% 
  mutate(cluster = if_else(cluster %in% c('CM', 'CD28+CXCR3+HLA-DR-', 'CD161+CCR6+'), cluster, 'Other'),
  # mutate(cluster = if_else(cluster %in% c('CM', 'CD28+CXCR3+HLA-DR-'), cluster, 'Other'),
  # mutate(cluster = if_else(cluster %in% c('CM', 'CD28+CXCR3+HLA-DR-', 'AT'), cluster, 'Other'),
         cluster = fct_relevel(cluster, 'Other', 'CD161+CCR6+', 'CD28+CXCR3+HLA-DR-','CM'),
         disease = fct_relevel(disease, 'Flu_infected', 'Flu_recovering', 'SSc', 'SLE', 'Ced')) %>% 
  filter(cluster != 'Other') %>%
  ggplot(aes(cluster, frac, fill = disease)) +
  geom_bar(stat = 'identity', position = 'dodge2') +
  # scale_y_continuous(limits = c(0, 0.5), oob = scales::squish) +
  # scale_fill_brewer(palette = 'Set2') +
  theme_minimal() +
  # scale_fill_brewer(type = 'qual')
  scale_fill_manual(values = cluster_colors)
ggsave('out/ai_vs_flu/cluster_enrichment_per_disease.png', height = 6)
ggsave('out/ai_vs_flu/cluster_enrichment_per_disease.pdf', height = 6)

cadf %>% 
  inner_join(samples) %>% 
  filter(sample_category == 'full') %>% 
  count(unique_name, cluster) %>% 
  add_count(unique_name, wt = n) %>% 
  mutate(frac = n / nn) %>% 
  rename(cells_in_cluster = n, cells_in_sample = nn) %>% 
  identity() ->
  counts_table_sample
counts_table_sample %>% 
  inner_join(samples) %>% 
  # filter(cluster != 'undef' & cluster != 'AT') %>%
  filter(cluster != 'undef') %>%
  mutate(prct = frac * 100,
         cluster = fct_relevel(cluster, 'CD28+CXCR3+HLA-DR-', 'CD161+CCR6+', 'CM'),
         disease = fct_relevel(disease, 'Flu_infected', 'Flu_recovering', 'Ced', 'SLE', 'SSc')) %>% 
  # filter(cluster != 'undef')  %>% 
  ggplot(aes(disease, frac, fill = cluster)) +
  geom_boxplot(outlier.color = NA) +#, coef = 0) +
  # geom_boxplot(outlier.color = NA) +
  # scale_y_log10() +
  # scale_fill_brewer(palette = 'Set1', type = 'qual') +
  scale_fill_manual(values = cluster_colors) +
  scale_y_continuous(labels = scales::percent, oob = scales::squish, limits = c(0, 0.50)) +
  theme_minimal() +
  # theme(#legend.position = 'bottom',
  #       axis.text.x = element_blank()) +
  labs(x = 'Disease and cell cluster',
       y = 'Percentage of activated CD4+ T-cells in cluster',
       title = 'Percentage of actived T-cells within selected cell clusters per disease group ')
ggsave('out/ai_vs_flu/cell_cluster_enrichment_boxplot.png', height = 5, width = 7)
ggsave('out/ai_vs_flu/cell_cluster_enrichment_boxplot.pdf', height = 5, width = 7)

counts_table_sample %>% 
  filter(cluster %in% c('CM', 'CD28+CXCR3+HLA-DR-', 'CD161+CCR6+')) %>% 
  inner_join(samples) %>% 
  mutate(disease = factor(disease),
         cells_outside_cluster = cells_in_sample - cells_in_cluster) %>% 
  identity() ->
  x
by_clust <- base::split(x, x$cluster)
by_clust <- by_clust %>% 
  map(function(x) base::split(x, x$disease))
diseases <- as.character(unique(x$disease))
clusters <- as.character(unique(x$cluster))

pvals <- list()
expand.grid(diseases, diseases)
combdf <- combn(diseases, 2) %>% 
  t() %>% 
  as_tibble() %>% 
  dplyr::rename(disease1 = V1, disease2 = V2) %>% 
  merge(tibble(cluster = clusters)) %>% 
  as_tibble()

pacman::p_load(purrrlyr)

f <- function(row) {
  df <- bind_rows(by_clust[[row$cluster]][[row$disease1]],
                  by_clust[[row$cluster]][[row$disease2]])
  x <- glmer(cbind(cells_in_cluster, cells_outside_cluster) ~ disease + (1|unique_name), family = binomial, data = df) %>%
    summary()
  pval <- x$coefficients[2,4]
  pval
}
clustdf <- combdf %>% 
  by_row(f, .collate = 'row', .to = 'pval')
  
clustdf %>% 
  filter(disease1 == 'Flu_infected') %>% 
  arrange(cluster, disease2) %>% 
  mutate(padj = p.adjust(pval)) %>% 
  mutate(sig_05 = padj < 0.05) 


all_count_df <- cadf %>% 
  inner_join(samples) %>% 
  filter(sample_category == 'full') %>% 
  count(unique_name, cell_clustering) %>% 
  rename(count = n)

all_count_df %>% 
  add_count(unique_name, wt = count) %>% 
  mutate(frac = count / n) %>% 
  inner_join(samples) %>% 
  ggplot(aes(disease, frac, fill = disease)) +
  geom_boxplot() +
  facet_wrap(~ cell_clustering) +
  scale_y_log10()

## Sketch
set.seed(1234)
sdf %>% 
  filter(sample_category == 'full') %>% 
  group_by(disease) %>% 
  sample_n(4000) %>% 
  ungroup() %>% 
  bind_rows(filter(sdf, sample_category == 'tetramer+')) %>% 
  # mutate(cluster = case_when(sample_category == 'tetramer+' ~ 'undef', 
  #                            # cluster == 'AT' ~ 'undef',
  #                            TRUE ~ cluster)) %>% 
  mutate(cluster = fct_relevel(cluster, 'undef', 'CD28+CXCR3+HLA-DR-', 'CD161+CCR6+', 'CM'),
         disease = fct_relevel(disease, 'Flu_infected', 'Flu_recovering', 'Ced', 'CeD_tet+', 'SLE', 'SSc')) %>% 
  identity() ->
  subdf
# Accent, Dark2, Paired, Pastel1, Pastel2, Set1, Set2, Set3 
subdf %>%
  ggplot(aes(tSNE1, tSNE2)) +
  geom_point(data = filter(subdf, cluster == 'undef'), color = 'gray', alpha = 0.4, size = 0.5) +
  geom_point(data = filter(subdf, cluster != 'undef'), mapping = aes(color = cluster), size = 1, alpha = 0.9) +
  # geom_density_2d(alpha = 0.5, n = 20) +
  facet_wrap(~ disease, ncol = 2) +
  theme_minimal() +
  scale_color_manual(values = cluster_colors)
  # scale_color_brewer(palette = 'Dark2', type = 'qual')
ggsave('out/ai_vs_flu/tSNE_by_group_with_clusters_equal_cell_count.png')
ggsave('out/ai_vs_flu/tSNE_by_group_with_clusters_equal_cell_count.pdf')

cadf %>% 
  mutate(cluster = if_else(sample_category == 'tetramer+', 'tetramer+', cluster)) %>% 
  select(-cell_n, -cell_clustering, -sample_group, -sample_category) %>% 
  gather(mark, value, -cluster, -unique_name) %>% 
  group_by(cluster, unique_name, mark) %>% 
  summarise(value = median(value)) %>% 
  ungroup() %>% 
  filter(cluster %in% c('CD28+CXCR3+HLA-DR-', 'CD161+CCR6+', 'CM', 'tetramer+')) %>%
  identity() ->
  expr_per_sample_cluster

expr_per_sample_cluster %>% 
  filter(!(mark %in% c('CD57', 'CD69', 'KLRG1'))) %>% 
  mutate(cluster = fct_relevel(cluster, 'CM', 'CD161+CCR6+', 'CD28+CXCR3+HLA-DR-', 'tetramer+'),
         mark = factor(mark, panel_ordered$fcs_name)) %>% 
  group_by(mark) %>% 
  mutate(zscore = as.vector(scale(value))) %>% 
  ungroup() %>% 
  group_by(cluster, mark) %>% 
  summarize(zscore = median(zscore)) %>% 
  ungroup() %>% 
  ggplot(aes(mark, cluster)) +
  geom_tile(aes(fill = zscore)) +
  # scale_fill_viridis_c(limit = c(-1.5, 1.5), oob = scales::squish) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_gradient2(low = 'dodgerblue4', mid = 'grey20', high = 'yellow2', limits = c(-1.5, 1.5), oob=scales::squish)

ggsave('out/ai_vs_flu/cluster_zscore_hmap.png', height = 4)
ggsave('out/ai_vs_flu/cluster_zscore_hmap.pdf', height = 4)
