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


pop_pars <- c(
        n_attendees=500,
        prev=0.05
        )


ninf <- reduce(lapply(t_test, get_n_infectious, 
	lod=35, se=0.95, inf_ct=inf_ct, maxct=maxct, params_df_slice, pop_pars=pop_pars, event_duration=3/24), bind_rows) %>% 
	pivot_wider(names_from="statistic", values_from=c("value")) %>%
	mutate(mean_smooth=predict(loess(mean~t, data=., span=0.6))) %>%
	mutate(lwr_smooth=predict(loess(lwr~t, data=., span=0.6))) %>%
	mutate(upr_smooth=predict(loess(upr~t, data=., span=0.6))) %>%
	mutate(t=-t)

ggplot() + 
	geom_ribbon(
	  data=ninf, 
	  aes(x=t, ymin=lwr_smooth, ymax=upr_smooth), alpha=0.2, fill="grey") + 
	geom_point(data=ninf, aes(x=t, y=lwr), size=0.1, alpha=0) +
	geom_point(data=ninf, aes(x=t, y=upr), size=0.1, alpha=0) +
	geom_point(data=ninf, aes(x=t, y=mean), size=0.1, alpha=0) + 
	geom_line(data=ninf, aes(x=t, y=mean), stat="smooth", method="loess", span=0.6) + 
	coord_cartesian(ylim=c(0,max(ninf$upr)), expand=FALSE) + 
	scale_x_reverse() + 
	theme_minimal() + 
	theme(text=element_text(size=18)) + 
	labs(title="Number infectious at event", subtitle="(90% pred. interval)", x="Days prior to event", y="Number infectious at event")