# not done.. 
# needs to read up on saving attributes (colnames)
# and easily read them back out again.
# I think this is as easy as using file handles...
pacman::p_load(
  tidyverse,
  rhdf5
)
source('20_experiment_and_sketch/load_fcs_and_rename_according_to_panel.R')

save_dataset <- function(M, path, h5_path) {
  h5write(M, path, h5_path)
  file <- H5Fopen(path)
  did <- H5Dopen(file, h5_path)
  h5writeAttribute(did, attr = colnames(M), name = 'colnames')
  H5Dclose(did)
  H5Fclose(file)
}

read_dataset <- function(path, h5_path) {
  M <- h5read(path, h5_path)
  colnames(M) <- h5readAttributes(path, h5_path)$colnames
  M
}

all_samples <- read_csv('input_data/renamed/sample_meta_full.csv')
outp <- 'out/cytof.h5'

if (file.exists(outp)) {
  file.remove(outp)    
}

h5createFile(outp)

samples_meta <- all_samples#head(all_samples, 1)
all_fcs <- read_all_fcs_files(sample_meta = samples_meta, fcs_trans = list(func = asinh, cofactor = 5))

h5createGroup(outp, 'meta')
h5write(samples_meta, outp, 'meta/samples')
h5write(all_fcs$cytof_panels, outp, 'meta/cytof_panels')
h5write(all_fcs$metals_that_varies, outp, 'meta/varying_metals')

h5createGroup(outp, 'samples')
all_fcs$df %>% 
  transmute(unique_name, expr = map(fcs, exprs)) %>% 
  identity() ->
  x
walk2(x$unique_name, x$expr, function(x, y) {
  h5_path <- glue('samples/{x}')
  print(h5_path)
  save_dataset(y, outp, h5_path)
  })
# 
# bc10f <- h5read(outp, 'samples/BC10_full')
# head(bc10f)
# h5ls(outp, all = TRUE)

h5ls(outp)
read_dataset(outp, 'samples/BC10_full') %>% 
  head()
h5readAttributes(outp, 'samples/BC10_full')$colnames
