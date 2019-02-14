# read FCS file, asign antibody panel (possible check with expected panel)
# and rename entries
pacman::p_load(
  tidyverse,
  flowCore,
  assertthat,
  glue,
  stringr,
  assertthat
)

filter <- dplyr::filter
select <- dplyr::select

# new_names = metals_with_desc_names$Antibody
#
#
subset_and_transform_fcs <- function(fcs, all_names, sub_markers, 
                                             transform_func = asinh, cofactor = 5) {
  # print(colnames(fcs))
  # print()
  # print(all_names)
  assert_that(length(colnames(fcs)) == length(all_names))
  assert_that(length(sub_markers) <= length(all_names))
  assert_that(length(intersect(all_names, sub_markers)) == length(sub_markers))
  
  colnames(fcs) <- all_names
  expr <- flowCore::exprs(fcs)
  expr <- transform_func(expr[, sub_markers, drop = F] / cofactor)
  flowCore::exprs(fcs) <- expr
  return(fcs)
}

#
# Get a full list of metals and antibodies per FCS file
#
extract_metals_and_ab_df <- function(fcs, fcsAB_to_AB) {
  desc <- fcs %>% 
    parameters() %>% 
    magrittr::extract2('desc') %>% 
    as.character()
  metals <- colnames(fcs)
  abmatch <- str_match(desc, '^[[:digit:]]+[[:alpha:]]+_(.+)')
  fcs_ab = abmatch[, 2]
  tibble(metal = metals,
         ab_desc = desc,
         fcs_AB = fcs_ab) %>% 
    left_join(fcsAB_to_AB, by = 'fcs_AB')
}

detect_cytof_panel <- function(varying_metals, fcs_df) {
  df <- fcs_df %>% 
    inner_join(varying_metals, by = 'metal_number') 
  num_auto <- str_detect(df$met_ab, df$Autoimmune) %>% 
    sum()
  num_ced <- str_detect(df$met_ab, df$CeD) %>% 
    sum()
  acceptable_cutoff <- nrow(varying_metals) - 2
  if (num_ced >= acceptable_cutoff) {
    return('CeD') # TODO -> don't return just rename according to panel!
  } else if (num_auto > acceptable_cutoff) {
    return('Autoimmune') # TODO -> don't return just rename according to panel!
  } else {
    print(df, n = Inf)
    stop(glue('FCS file matches neither CeD ({num_ced}/{acceptable_cutoff})  nor Autimmune ({num_auto}/{acceptable_cutoff}) panel'))
  }
}

# uses the following values from the global environment.
# `cytof_panel`: list of antibodies and metals used in various panels
# `metals_that_varies` metals tagging different antibodies in different panels
#
# Could add a parameter that specifies the cytofpanel used (still check if sensible)
read_all_fcs_files <- function(only_functional_AB = TRUE, sample_meta = 'input_data/renamed/sample_meta_full.csv',
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
    
    f <- partial(subset_and_transform_fcs,
                 all_names = metals_with_desc_names$Antibody, sub_markers = markers, 
                 transform_func = fcs_trans$func, cofactor = fcs_trans$cofactor)
    fcs_raw %>% 
      f() %>% 
      with_R_friendly_names()
  }
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
  sample_meta %>% 
    mutate(fcs = map(path, read_fcs)) %>% 
    identity() ->
    sdf
  retval <- list(df = sdf,
                 cytof_panels = cytof_panels,
                 metals_that_varies = metals_that_varies)
  return (retval)
}

select_fcs_subset <- function(xs, markers) {
  f <- function(x) {
    expr <- flowCore::exprs(x)
    flowCore::exprs(x) <- expr[, markers, drop = F]
    x
  }
  xs %>% 
    map(f)
}

flowset_of_subgroup <- function(fcs_df, markers = NULL) {
  assert_that(has_name(fcs_df, 'unique_name'))
  if (is.null(markers)) {
    markers <- get_common_markers(fcs_df$fcs)
  }
  assert_that(is_bare_character(markers))
  fcs_df %>% 
    mutate(fcs = select_fcs_subset(fcs, markers),
           fcs = setNames(fcs, unique_name)) %>% 
    magrittr::extract2('fcs') %>% 
    as('flowSet')
}

get_common_markers <- function(xs) {
  nsamples <- length(xs)
  xs %>% 
    map(colnames) %>% 
    unlist() %>% 
    tibble(antibody = .) %>% 
    count(antibody) %>% 
    filter(n == nsamples) %>% 
    .$antibody
}