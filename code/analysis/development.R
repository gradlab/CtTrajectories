t_test <- seq(-3, 0, 1/24) 

maxct <- 40
inf_ct <- 30

params_df_slice <- params_df %>% 
	filter(dp>(maxct-inf_ct)) %>%
	sample_n(1000)

eff_se <- unlist(lapply(t_test, get_effective_sensitivity, 
    lod=35, se=0.95, inf_ct=inf_ct, maxct=maxct, params_df_slice, event_duration=3/24))

ggplot(data=tibble(x=-t_test,y=eff_se), aes(x=x, y=y)) + 
  geom_point(size=0.1, alpha=0.2) + 
  geom_line(stat="smooth", method="loess", span=0.6) + 
  scale_y_continuous(limits=c(0,1)) + 
  scale_x_reverse() + 
  labs(title="Effective sensitivity", x="Days prior to event", y="Effective sensitivity") + 
  theme_minimal() + 
  theme(text=element_text(size=18))


