ceddf <- read_csv('out/paper_figs/main/prediction_figure/model_with_CD62L_CCR7_parameter_importance.csv')
p_ced <- ceddf %>% 
  mutate(antibody = fct_reorder(antibody, MeanDecreaseGini)) %>% 
  ggplot(aes(MeanDecreaseGini, antibody)) +
  geom_point(shape = 5) +
  theme_minimal()
ggsave('out/paper_figs/main/prediction_figure/model_with_CD62L_CCR7_parameter_importance.png', p_ced, width = 5, height = 6)
ggsave('out/paper_figs/main/prediction_figure/model_with_CD62L_CCR7_parameter_importance.pdf', p_ced, width = 5, height = 6)

aidf <- read_csv('out/paper_figs/main/prediction_figure/model_without_CD62L_CCR7_parameter_importance.csv')
p_ai <- aidf %>% 
  mutate(antibody = fct_reorder(antibody, MeanDecreaseGini)) %>% 
  ggplot(aes(MeanDecreaseGini, antibody)) +
  geom_point(shape = 5) +
  theme_minimal()
ggsave('out/paper_figs/main/prediction_figure/model_without_CD62L_CCR7_parameter_importance.png', p_ai, width = 5, height = 6)
ggsave('out/paper_figs/main/prediction_figure/model_without_CD62L_CCR7_parameter_importance.pdf', p_ai, width = 5, height = 6)
