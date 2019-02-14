pacman::p_load(
  tidyverse,
  glue,
  rhdf5,
  viridis
)
SEED <- 12345
all_samples <- read_csv('input_data/renamed/sample_meta_full.csv')

ch_donors <- c('CD1300')
ucd_donor <- c('CD1467')

all_samples %>% 
  filter(biosource == 'PBMC') %>% 
  filter(donor %in% ch_donors | sample_category == 'tetramer+' & donor == ucd_donor) %>% 
  filter(sample_category != 'tetramer-') %>% 
  arrange(sample_group, sample_category) %>% 
  identity() ->
  samples

read_dataset <- function(path, h5_path) {
  M <- h5read(path, h5_path)
  colnames(M) <- h5readAttributes(path, h5_path)$colnames
  M
}

read_sample <- function(unique_names, hdf_path = 'out/cytof.h5') {
  unique_names %>% 
    map(function(u) {
      h5_path <- glue("samples/{u}")
      read_dataset(hdf_path, h5_path)
    }) %>% 
    set_names(unique_names)
}

M <- read_sample(samples$unique_name)
rdf <- M %>% 
  map_dfr(as_tibble, .id = 'unique_name')

tetp_cells <- samples %>% 
  filter(sample_category == 'tetramer+') %>% 
  semi_join(rdf, .)
bg_cells <- samples %>% 
  filter(sample_category == 'tetramer+') %>% 
  anti_join(rdf, .)
bg_cells %>% gather(marker, expr, -unique_name) %>% filter(marker == 'CD127') %>% ggplot(aes(marker, expr)) + geom_boxplot(aes(fill = unique_name), outlier.color = NA) +
  ggtitle('CD127 in selected samples')

set.seed(SEED)
tsnedf <- bg_cells %>% 
  group_by(unique_name) %>% 
  sample_n(1e4) %>% 
  ungroup() %>% 
  bind_rows(tetp_cells)
nrow(tsnedf)

##
## Compute t-SNE
##
M_class <- tsnedf$unique_name
M <- as.matrix(select(tsnedf, -unique_name))

devtools::install_github("jkrijthe/Rtsne", ref = 'openmp')
library(Rtsne)

get_tsne <- function(M, SEED = 12345, perplexity = 30) {
  set.seed(SEED)
  Rtsne(M, check_duplicates = FALSE, pca = FALSE, num_threads = 6, perplexity = perplexity)
}
path_tsne_df <- 'out/rechallenge/tsne_tetp_ucd_chal_gfd.csv.gz'
if (! file.exists(path_tsne_df)) {
  tsne_raw <- get_tsne(M, perplexity = 40, SEED = 54321)
  tsne_df <- as_tibble(tsne_raw$Y) %>% 
    set_names(c('x', 'y')) %>% 
    mutate(unique_name = M_class) %>% 
    bind_cols(as_tibble(M)) %>% 
    inner_join(samples %>% 
                 dplyr::select(unique_name, donor, disease_status, sample_category))
  tsne_df %>% 
    write_csv(path_tsne_df)
} else {
  tsne_df <- read_csv(path_tsne_df)
}

tsne_df %>% 
  filter(sample_category == 'full') %>% 
  group_by(unique_name) %>% 
  sample_n(5000) %>%
  identity() ->
  sub_pre
tsne_df %>% 
  filter(sample_category != 'full') %>% 
  bind_rows(sub_pre) %>% 
  identity() ->
  df

df %>% 
  mutate(category = str_c(donor, disease_status, sep = '_')) %>% 
  ggplot(aes(x, y)) +
  geom_point(data = filter(df, sample_category != 'tetramer+'), color = 'gray', alpha = 1) +
  geom_point(data = filter(df, sample_category == 'tetramer+'), mapping = aes(color = disease_status), size = 2) +
  scale_color_manual(values = c('Challenge' = 'firebrick', 'GFD' = 'black', 'UCD' = '#f36d16')) +
  theme_minimal() +
  theme(legend.position = '',
        axis.title = element_blank(),
        axis.text = element_blank())
ggsave('out/rechallenge/tsne_tetp_ucd_chal_gfd.png', bg = 'transparent')
ggsave('out/rechallenge/tsne_tetp_ucd_chal_gfd.pdf', bg = 'transparent')

tsne_df %>% 
  gather(marker, expression, -unique_name, -donor, -x, -y, -disease_status, -sample_category) %>% 
  group_by(marker) %>% 
  mutate(frac_of_max = expression / max(expression)) %>% 
  ungroup() %>% 
  # filter(str_detect(marker, 'PD')) %>%
  # head(20000) %>%
  ggplot(aes(x, y)) +
  # geom_point(aes(color = scaled_intensity)) +
  # geom_point(aes(color = scaled_intensity), size = 1) +
  stat_summary_2d(aes(z = frac_of_max), fun = median, bins = 70) +
  # annotate("rect", xmin = tlim$PBMC$xmin, xmax = tlim$PBMC$xmax, 
  #            ymin = tlim$PBMC$ymin, ymax = tlim$PBMC$ymax,
  #            alpha = .1,
  #          color = 'black') +
  facet_wrap(~ marker, ncol = 6) +
  theme_minimal() +
  # scale_fill_viridis() +
  scale_fill_continuous(type = 'viridis', limits = c(0, 1)) +
  labs(fill = 'Scaled Intensity') +
  coord_fixed()
  
##
##
##
all_samples %>% 
  filter(sample_category == 'tetramer+' & disease == 'Ced' & biosource == 'PBMC') %>%
  identity() ->
  tetp_samples


tpMs <- read_sample(tetp_samples$unique_name) %>% 
  set_names(tetp_samples$unique_name)
prob_samples <- tpMs %>% 
  keep(function(M) 'TIGIT' %in% colnames(M)) %>% 
  names()
tetp_samples %>% 
  filter(unique_name %in% prob_samples) %>% 
  glimpse()
tpdf <- read_sample(tetp_samples$unique_name) %>% 
  discard(function(M) 'TIGIT' %in% colnames(M)) %>% 
  map_dfr(as_tibble, .id = 'unique_name') %>% 
  gather(marker, expr, -unique_name) %>% 
  inner_join(tetp_samples)
library(visdat)
exprdf <- tpdf %>% 
  mutate(expr = sinh(expr) * 5) %>% 
  group_by(disease_status, marker, unique_name) %>% 
  summarize(expr = mean(expr)) %>% 
  summarize(expr = mean(expr)) 

exprdf %>% 
  ggplot(aes(disease_status, marker, fill = expr)) +
  geom_tile() +
  scale_fill_viridis(limits = c(0, 75), oob = scales::squish)
ggsave('out/rechallenge/tetp_absolute_expr.png')
  
mark_order <- read_csv('input_data/original/asbjorn_markers_specific.csv') %>% 
  mutate(marker = str_replace_all(common_name, '[- ]', '_')) %>% 
  filter(category == 'CyTOF') %>% 
  dplyr::select(-gene_name, -category)

exprdf %>% 
  spread(disease_status, expr) %>% 
  mutate(ch_vs_ucd = log2(Challenge / UCD),
         ch_vs_gfd = log2(Challenge / GFD)) %>% 
  select(-Challenge, -GFD, -UCD) %>% 
  gather(comparison, logfc, -marker) %>% 
  inner_join(mark_order) %>% 
  mutate(AB = fct_relevel(common_name, rev(mark_order$common_name))) %>% 
  ggplot(aes(comparison, AB, fill = logfc)) +
  geom_tile() +
  scale_fill_gradient2(low = 'dodgerblue4', mid = 'grey20', high = 'yellow2', limits = c(-4, 4), oob = scales::squish) +
  # scale_fill_viridis(limits = c(-4, 4), oob = scales::squish) +
  theme_minimal()

ggsave('out/rechallenge/tetp_logfc_challenge_vs_gfd_or_ucd.png')
ggsave('out/rechallenge/tetp_logfc_challenge_vs_gfd_or_ucd.pdf')
