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
SEED <- 1234567
cnum <- 30# number of clusters
outdir <- '75_flu_vs_ai_v2/out/clustering'
cluster_dir <- glue('{outdir}/seed_{SEED}')
cluster_path <- file.path(cluster_dir, glue('cluster_c{cnum}_seed{SEED}.csv.gz'))
dir.create(cluster_dir, showWarnings = FALSE, recursive = TRUE)
samples <- read_csv('75_flu_vs_ai_v2/out/samples.csv')
panel <- read_csv('75_flu_vs_ai_v2/out/predmarkers.csv') %>% 
  filter(fcs_name != 'CD45RA')
panel_order <- read_csv('input_data/original/asbjorn_markers_specific.csv')
panel %>% 
  inner_join(panel_order, ., by = c('common_name' = 'antibody' )) %>% 
  identity() ->
  panel_ordered
adf <- read_csv('75_flu_vs_ai_v2/out/clust_subsampled_activated_tcells.csv.gz')
tsne_idx <- read_csv('75_flu_vs_ai_v2/out/tsne_idx.csv')

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
  som <- BuildSOM(fsom, colsToUse = panel$fcs_name, xdim = 20, ydim = 20, rlen = 20) 
  # colsToUse = discard(panel$fcs_name, ~ . %in% c('CTLA_4'CD73')))
  
  codes <- som$map$codes
  plot_outdir <- file.path(cluster_dir, "consensus")
  
  mc <- ConsensusClusterPlus(t(codes), maxK = cnum, reps = 10000,
                             pItem = 0.9, pFeature = 1, title = plot_outdir, plot = "png",
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
