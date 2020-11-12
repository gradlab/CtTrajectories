# =============================================================================
# Initialize effective sensitivity parameters
# =============================================================================

t_test <- seq(-3, 0, 1/24) 

params_df_slice <- params_df %>% 
	filter(dp>(maxct-inf_ct)) %>%
	sample_n(2500)

# =============================================================================
# Plot effective sensitivity
# =============================================================================

eff_se_40 <- unlist(lapply(t_test, get_effective_sensitivity, 
    lod=40, se=0.99, inf_ct=inf_ct, maxct=maxct, params_df_slice, event_duration=3/24))

eff_se_35 <- unlist(lapply(t_test, get_effective_sensitivity, 
    lod=35, se=0.95, inf_ct=inf_ct, maxct=maxct, params_df_slice, event_duration=3/24))

eff_se_df <- tibble(t=-t_test, eff_se_40=eff_se_40, eff_se_35=eff_se_35)

fig_eff_se <- eff_se_df %>% 
	pivot_longer(-t) %>%
	ggplot(aes(x=t, y=value, col=name)) + 
		geom_point(size=0.1, alpha=0.0) + 
		geom_line(stat="smooth", method="loess", span=0.6) + 
		scale_y_continuous(limits=c(0,1)) + 
		scale_x_reverse() + 
		labs(x="Days prior to event", y="Effective sensitivity") + 
		scale_color_manual(values=c(eff_se_40="red",eff_se_35="blue")) + 
		theme_minimal() + 
		theme(text=element_text(size=18), legend.position="none")

# =============================================================================
# Plot number infected
# =============================================================================

pop_pars <- c(
        n_attendees=1000,
        prev=0.02
        )

ninf_40 <- reduce(lapply(t_test, get_n_infectious, 
	lod=40, se=0.99, inf_ct=inf_ct, maxct=maxct, params_df_slice, pop_pars=pop_pars, event_duration=3/24), bind_rows) %>% 
	pivot_wider(names_from="statistic", values_from=c("value")) %>%
	mutate(mean_smooth=predict(loess(mean~t, data=., span=0.6))) %>%
	mutate(lwr_smooth=predict(loess(lwr~t, data=., span=0.6))) %>%
	mutate(upr_smooth=predict(loess(upr~t, data=., span=0.6))) %>%
	mutate(t=-t)

ninf_35 <- reduce(lapply(t_test, get_n_infectious, 
	lod=35, se=0.95, inf_ct=inf_ct, maxct=maxct, params_df_slice, pop_pars=pop_pars, event_duration=3/24), bind_rows) %>% 
	pivot_wider(names_from="statistic", values_from=c("value")) %>%
	mutate(mean_smooth=predict(loess(mean~t, data=., span=0.6))) %>%
	mutate(lwr_smooth=predict(loess(lwr~t, data=., span=0.6))) %>%
	mutate(upr_smooth=predict(loess(upr~t, data=., span=0.6))) %>%
	mutate(t=-t)

ninf <- rbind(mutate(ninf_40, test="pcr"), mutate(ninf_35, test="rapid"))

fig_ninf <- ggplot() + 
	geom_ribbon(
	  data=ninf, 
	  aes(x=t, ymin=lwr_smooth, ymax=upr_smooth, fill=test), alpha=0.2) + 
	geom_point(data=ninf, aes(x=t, y=lwr, col=test), size=0.1, alpha=0) +
	geom_point(data=ninf, aes(x=t, y=upr, col=test), size=0.1, alpha=0) +
	geom_point(data=ninf, aes(x=t, y=mean, col=test), size=0.1, alpha=0) + 
	geom_line(data=ninf, aes(x=t, y=mean, col=test), stat="smooth", method="loess", span=0.6) + 
	coord_cartesian(ylim=c(0,max(ninf$upr)), expand=FALSE) + 
	scale_x_reverse() + 
	scale_color_manual(values=c(pcr="red",rapid="blue")) + 
	theme_minimal() + 
	theme(text=element_text(size=18), legend.position="none") + 
	labs(x="Days prior to event", y="Number infectious at event")