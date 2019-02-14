pacman::p_load(
  glue
)
df %>% 
  filter(disease_status == 'UCD') %>% 
  filter(biosource == 'SCS') %>% 
  anti_join(tibble(variable = skip_in_gut)) %>% 
  # mutate(value = asinh(value) / 5.0) %>%
  # identity() ->
  # xdf
  spread(sample_category, value) %>% 
  group_by(variable) %>% 
  summarize(pval = t.test(`tetramer+`, `tetramer-`, paired = T)$p.value) %>%
  # summarize(pval = wilcox.test(`tetramer+`, `tetramer-`, paired = T)$p.value) %>% 
  arrange(pval) %>% 
  mutate(FDR = p.adjust(pval, method = 'BH')) %>% 
  identity() ->
  abs_pval_scs

df %>% 
  filter(disease_status == 'UCD') %>% 
  filter(biosource == 'SCS') %>% 
  anti_join(tibble(variable = skip_in_gut)) %>% 
  mutate(value = asinh(value) / 5.0) %>%
  # identity() ->
  # xdf
  spread(sample_category, value) %>% 
  group_by(variable) %>% 
  summarize(pval = t.test(`tetramer+`, `tetramer-`, paired = T)$p.value) %>%
  # summarize(pval = wilcox.test(`tetramer+`, `tetramer-`, paired = T)$p.value) %>% 
  arrange(pval) %>% 
  mutate(FDR = p.adjust(pval, method = 'BH')) %>% 
  identity() ->
  asinh_pval_scs

xdf %>% 
  inner_join(asinh_pval_scs) %>% 
  mutate(variable = glue("{variable} {formatC(FDR, digits = 2, format='f')}")) %>% 
  mutate(value = asinh(value) / 5.0) %>%
  ggplot(aes(sample_category, value, group = donor, color = donor)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ variable) %>% 
  identity() ->
  p_asinh

xdf %>% 
  inner_join(abs_pval_scs) %>% 
  mutate(variable = glue("{variable} {formatC(FDR, digits = 2, format='f')}")) %>% 
  # mutate(value = asinh(value) / 5.0) %>%
  ggplot(aes(sample_category, value, group = donor, color = donor)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ variable) %>% 
  identity() ->
  p_abs

ggsave('out/tmp/FC_fdr_asinh.png', p_asinh)
ggsave('out/tmp/FC_fdr_abs.png', p_abs)

xdf %>% 
  filter(variable == 'CD39') %>% 
  mutate(value = asinh(value) / 5.0) %>%
  ggplot(aes(sample_category, value, group = donor, color = donor)) +
  geom_point() +
  geom_line()
