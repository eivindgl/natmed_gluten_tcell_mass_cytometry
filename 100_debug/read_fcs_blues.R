# debugging mysterious bug that turned out to be ggplot::exprs shadowing flowCore::exprs...
# TODO for later. any way to isolate dependencies? or just always use explicit names?

st_fcs <- function(fcs, all_names, sub_markers, 
                                             transform_func = asinh, cofactor = 5) {
  # print(colnames(fcs))
  # print()
  # print(all_names)
  assert_that(length(colnames(fcs)) == length(all_names))
  assert_that(length(sub_markers) <= length(all_names))
  assert_that(length(intersect(all_names, sub_markers)) == length(sub_markers))
  
  colnames(fcs) <- all_names
  expr <- exprs(fcs)
  expr <- transform_func(expr[, sub_markers, drop = F] / cofactor)
  exprs(fcs) <- expr
  return(fcs)
}

read_fcs<- function(path) {
    
    with_R_friendly_names <- function(x) {
      colnames(x) <- str_replace_all(colnames(x), "[- ]", "_")
      return(x)
    }
    
    print(glue('reading FCS file {path}'))
    fcs_raw <- read.FCS(path, transformation = FALSE,
             truncate_max_range = FALSE)
    extract_metals_and_ab_df(fcs_raw, fcsAB_to_AB) %>% 
      select(metal, Antibody) %>% 
      identity() ->
      metals_with_desc_names
    
    markers <- ab_meta %>% 
      inner_join(metals_with_desc_names, by = 'Antibody')
    if (only_functional_AB) {
      markers <- markers %>%  
        filter(functional)
    } else {
      markers <- markers %>% 
        filter(!(Antibody %in% c('CD3', 'CD4', 'DNA1', 'DNA2')))
    }
    markers <- unique(markers$Antibody)
    
    f <- partial(st_fcs,
                 all_names = metals_with_desc_names$Antibody, sub_markers = markers, 
                 transform_func = fcs_trans$func, cofactor = fcs_trans$cofactor)
    fcs_raw %>% 
      f() %>% 
      with_R_friendly_names()
  }

read_fcs_files <- function(only_functional_AB = TRUE, sample_meta = 'input_data/renamed/sample_meta_full.csv',
                               fcs_trans = list(func = asinh, cofactor = 5)) {
  cytof_panels <- read_csv('input_data/original/cytof_panels.csv') %>% 
    mutate(metal_number = str_extract(Metal, '\\d+'))
  number_of_panels <- nrow(distinct(cytof_panels, Panel))
  assert_that(are_equal(number_of_panels, 2)) # , msg = 'This code assumes two cytof panels, abort if not the case'
  
  #
  ## Use subset that varies between panels to detect correct panel per FCS files
  #
  cytof_panels %>% 
    spread(Panel, Antibody) %>% 
    filter(CeD != Autoimmune) %>% 
    identity() ->
    metals_that_varies
  fcsAB_to_AB <- read_csv('input_data/original/fcs_AB_to_Antibody.csv')
  
  ab_meta <- read_csv('input_data/original/antibody_meta.csv')
  
  ## Test calls
  # x <- read_fcs('input_data/original/SCS/170519 SCS 2nd gating for R/UCD1462D.A Tetramer POS.fcs')
  # read_fcs('input_data/original/Autoimmune control and Flu CD4/Flu samples/22-002-01-12 Flu POS CD4 Phenotype.fcs')
  # read_fcs("input_data/original/SCS/170602 SCS Multiple pat for R/TCD1375D.D Tetramer POS.fcs") # was problematic - 1 cell
  # read_fcs('input_data/original/PBMC/170425 UCDpat 3 PBM/170425 BC11 CD4.fcs') # was problematic - 1 varying marker not matching (TODO check why)
  
  if (is_bare_character(sample_meta)) {
    sample_meta <- read_csv(sample_meta)
  }
  assert_that(is.data.frame(sample_meta))
  assert_that(nrow(sample_meta) > 0)
  assert_that(all(file.exists(sample_meta$path)))
  fcs = read_fcs(sample_meta$path[1])
  return(123)
  sample_meta %>% 
    mutate(fcs = map(path, read_fcs)) %>% 
    identity() ->
    sdf
  retval <- list(df = sdf,
                 cytof_panels = cytof_panels,
                 metals_that_varies = metals_that_varies)
  return (retval)
}
read_fcs(ced_samples_meta$path[1])

         
x <- read.FCS("input_data/original/PBMC/180302 PBMC Bergen/UCD1349 CD4 PRE Enriched.fcs", transformation = F, truncate_max_range = F)
flowCore::exprs(x)
