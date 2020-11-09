# =============================================================================
# Import: 
# =============================================================================

library(tidyverse) 
source('code/utilities/utils_analysis.R')

# =============================================================================
# Generate best-fit distributions: 
# =============================================================================

library(fitdistrplus)

# Find distributions that align with the wp and dp: 
wp_gamma_fit <- fitdist(params_df$wp, "gamma") 
wr_gamma_fit <- fitdist(params_df$wr, "gamma") 
infection_duration_gamma_fit <- fitdist(params_df$wp+params_df$wr, "gamma") 
peak_ct_normal_fit <- fitdist(global_pars[["lod"]]-params_df$dp, "norm") 

detach(package:fitdistrplus, unload=TRUE)
detach(package:MASS, unload=TRUE)

wp_gamma_fit_curve <- geom_line(data=make_gamma_prior_df(shape=wp_gamma_fit$estimate[["shape"]], rate=wp_gamma_fit$estimate[["rate"]], pmin=0, pmax=0.999, step=0.1), aes(x=x, y=density), col="black")
wr_gamma_fit_curve <- geom_line(data=make_gamma_prior_df(shape=wr_gamma_fit$estimate[["shape"]], rate=wr_gamma_fit$estimate[["rate"]], pmin=0, pmax=0.999, step=0.1), aes(x=x, y=density), col="black")
infection_duration_gamma_fit_curve <- geom_line(data=make_gamma_prior_df(shape=infection_duration_gamma_fit$estimate[["shape"]], rate=infection_duration_gamma_fit$estimate[["rate"]], pmin=0, pmax=0.999, step=0.1), aes(x=x, y=density), col="black")
peak_ct_fit_curve <- geom_line(data=make_normal_prior_df(mean=peak_ct_normal_fit$estimate[["mean"]], sd=peak_ct_normal_fit$estimate[["sd"]], pmin=0.001, pmax=0.999, step=0.1), aes(x=x, y=density), col="black")

# =============================================================================
# Generate figures
# =============================================================================

figlist_ct_dat_refined <- with(as.list(global_pars),{
	ct_dat_clean %>% 
	clean_person_id() %>%
	left_join(kept_ids, by="Person.ID") %>%
	replace_na(list(CT.Mean=lod, kept=0)) %>% 
	mutate(kept=case_when(kept==1~"Acute", TRUE~"Persistent")) %>%
	mutate(plotgroup = floor((Person.ID.Clean-1)/9)+1) %>%
	split(.$plotgroup) %>% 
	map(~ ggplot(., aes(x=Date.Index, y=CT.Mean, col=kept)) + 
		geom_point(size=0.5) + 
		geom_line() + 
		scale_y_reverse(limits=c(lod, min(ct_dat_clean$CT.Mean, na.rm=TRUE))) + 
		scale_color_manual(values=c("Acute"="red", "Persistent"="black")) + 
		theme_minimal() + 
		theme(text=element_text(size=14), legend.title=element_blank(), legend.position="none") + 
		labs(x="Days since first sample", y="Ct") + 
		facet_wrap(~Person.ID))
	})

fig_ct_fit <- plot_ct_fit(params_df, global_pars, indiv_data, ctalpha=0.02)
fig_ct_fit_symp <- plot_ct_fit_symp(params_df, global_pars, indiv_data, ctalpha=0.02)
fig_ct_fit_list <- plot_ct_fit_panes(params_df, global_pars, indiv_data, ctalpha=0.02)

fig_peak_ct <- with(as.list(c(global_pars, prior_pars)),{
	params_df %>% 
	ggplot(aes(x=lod-dp)) + 
		# Plot posterior density: 
		kernelplot + 
		scale_x_continuous(limits=c(0,lod)) + 
		labs(x="Peak Ct", y="Density") + 
		y_ticks_off + 
		grid_off + 
		facet_wrap(~id)
	})

# Plot the peak time posterior distribution (relaative to min Ct time):
fig_t_peak <- with(as.list(prior_pars),{params_df %>% 
	left_join(data.frame(
		id_clean=as.character(1:n_indiv), 
		tp_prior=0, stringsAsFactors=FALSE), by="id_clean") %>% 
	mutate(tp_normed = tp-tp_prior) %>% 
	ggplot(aes(x=tp_normed)) + 
		# Plot posterior density: 
		kernelplot + 
		labs(x="Peak time relative to min Ct (days)", y="Density") + 
		y_ticks_off + 
		grid_off + 
		facet_wrap(~id)
	})

# Plot the onset-to-peak time posterior distribution:
fig_t_onset_to_peak <- with(as.list(prior_pars),{params_df %>% 
	ggplot(aes(x=wp)) + 
		# Plot posterior density: 
		kernelplot + 
		scale_x_continuous(limits=c(0,wpmax)) + 
		labs(x="Proliferation stage duration (days)", y="Density") + 
		y_ticks_off + 
		grid_off + 
		# xlim(0,30) + 
		facet_wrap(~id)
	})

# Plot the peak-to-recovery time posterior distribution:
fig_t_peak_to_recovery <- with(as.list(prior_pars),{params_df %>% 
	ggplot(aes(x=wr)) + 
		# Plot posterior density: 
		kernelplot + 
		scale_x_continuous(limits=c(0,wrmax)) + 
		labs(x="Clearance stage duration (days)", y="Density") + 
		y_ticks_off + 
		grid_off + 
		facet_wrap(~id)
	})

# Plot the duration of infection distribution:
fig_duration_of_infection <- with(as.list(prior_pars),{params_df %>% 
	ggplot(aes(x=wr+wp)) + 
		# Plot posterior density: 
		kernelplot + 
		scale_x_continuous(limits=c(0,wpmax+wrmax)) + 
		labs(x="Duration of infection (days)", y="Density") + 
		y_ticks_off + 
		grid_off + 
		facet_wrap(~id)
	})

# Plot the process noise: 
fig_sigma <- with(as.list(prior_pars), {data.frame(sigma=params$sigma) %>% 
	ggplot(aes(x=sigma)) + 
	# Plot posterior density: 
	kernelplot + 
	# Plot prior: 
	geom_line(
			data=make_cauchy_prior_df(location=0,scale=sigma_prior_scale,xmin=0,xmax=sigma_max,step=0.01),
			aes(x=x, y=density), col="cornflowerblue", alpha=0.5) + 
	y_ticks_off + 
	grid_off + 
	labs(x="Process noise SD (Ct)", y="Density")
	})

if(current_pars[["symptom_treatment"]]=="split"){
	fig_dpmean_withprior <- shared_params_df %>% 
		mutate(dpmean_symp=dpmeanS) %>% 
		mutate(dpmean_asymp=dpmeanA) %>% 
		select(dpmean_symp, dpmean_asymp) %>% 
		pivot_longer(everything()) %>% 
		ggplot(aes(x=global_pars[["lod"]]-value)) + 
			geom_histogram(aes(y=..density.., fill=name), alpha=0.2, position="identity", bins=50) + 
			geom_density(aes(col=name), adjust=2) + 
			geom_line(data=make_normal_prior_df(mean=prior_pars$dpmean_prior, sd=prior_pars$dpsd_prior, pmin=0.001, pmax=0.999, step=0.1), aes(x=x, y=density), col="black", linetype="dashed") + 
			scale_x_continuous(limits=c(0,NA)) + 
			scale_color_manual(values=c("dpmean_symp"="red","dpmean_asymp"="blue")) + 
			scale_fill_manual(values=c("dpmean_symp"="red","dpmean_asymp"="blue")) + 
			labs(x="Mean peak Ct", y="Density") + 
			theme_minimal() + 
			theme(legend.position="none", text=element_text(size=18)) + 
			y_ticks_off + 
			grid_off
} else {
	fig_dpmean_withprior <- shared_params_df %>% 
		ggplot(aes(x=global_pars[["lod"]]-dpmean)) + 
			geom_histogram(aes(y=..density..), alpha=0.2, position="identity", bins=50) + 
			geom_density(adjust=2) + 
			geom_line(data=make_normal_prior_df(mean=prior_pars$dpmean_prior, sd=prior_pars$dpsd_prior, pmin=0.001, pmax=0.999, step=0.1), aes(x=x, y=density), col="black", linetype="dashed") + 
			scale_x_continuous(limits=c(0,NA)) + 
			labs(x="Mean peak Ct", y="Density") + 
			theme_minimal() + 
			theme(legend.position="none", text=element_text(size=18)) + 
			y_ticks_off + 
			grid_off
}


if(current_pars[["symptom_treatment"]]=="split"){
	fig_wpmean_withprior <- shared_params_df %>% 
		mutate(wpmean_symp=wpmeanS) %>% 
		mutate(wpmean_asymp=wpmeanA) %>% 
		select(wpmean_symp, wpmean_asymp) %>% 
		pivot_longer(everything()) %>% 
		ggplot(aes(x=value)) + 
			geom_histogram(aes(y=..density.., fill=name), alpha=0.2, position="identity", bins=50) + 
			geom_density(aes(col=name), adjust=2) + 
			geom_line(data=make_normal_prior_df(mean=prior_pars$wpmean_prior, sd=prior_pars$wpsd_prior, pmin=0.001, pmax=0.999, step=0.1), aes(x=x, y=density), col="black", linetype="dashed") + 
			scale_x_continuous(limits=c(0,NA)) + 
			scale_color_manual(values=c("wpmean_symp"="red","wpmean_asymp"="blue")) + 
			scale_fill_manual(values=c("wpmean_symp"="red","wpmean_asymp"="blue")) + 
			labs(x="Mean proliferation stage duration (days)", y="Density") + 
			theme_minimal() + 
			theme(legend.position="none", text=element_text(size=18)) + 
			y_ticks_off + 
			grid_off
	} else {
	fig_wpmean_withprior <- shared_params_df %>% 
		ggplot(aes(x=wpmean)) + 
			geom_histogram(aes(y=..density..), alpha=0.2, position="identity", bins=50) + 
			geom_density(adjust=2) + 
			geom_line(data=make_normal_prior_df(mean=prior_pars$wpmean_prior, sd=prior_pars$wpsd_prior, pmin=0.001, pmax=0.999, step=0.1), aes(x=x, y=density), col="black", linetype="dashed") + 
			scale_x_continuous(limits=c(0,NA)) + 
			labs(x="Mean proliferation stage duration (days)", y="Density") + 
			theme_minimal() + 
			theme(legend.position="none", text=element_text(size=18)) + 
			y_ticks_off + 
			grid_off
	}


if(current_pars[["symptom_treatment"]]=="split"){
	fig_wrmean_withprior <- shared_params_df %>% 
		mutate(wrmean_symp=wrmeanS) %>% 
		mutate(wrmean_asymp=wrmeanA) %>% 
		select(wrmean_symp, wrmean_asymp) %>% 
		pivot_longer(everything()) %>% 
		ggplot(aes(x=value)) + 
			geom_histogram(aes(y=..density.., fill=name), alpha=0.2, position="identity", bins=50) + 
			geom_density(aes(col=name), adjust=2) + 
			geom_line(data=make_normal_prior_df(mean=prior_pars$wrmean_prior, sd=prior_pars$wrsd_prior, pmin=0.001, pmax=0.999, step=0.1), aes(x=x, y=density), col="black", linetype="dashed") + 
			scale_x_continuous(limits=c(0,NA)) + 
			scale_color_manual(values=c("wrmean_symp"="red","wrmean_asymp"="blue")) + 
			scale_fill_manual(values=c("wrmean_symp"="red","wrmean_asymp"="blue")) + 
			labs(x="Mean clearance stage duration (days)", y="Density") + 
			theme_minimal() + 
			theme(legend.position="none", text=element_text(size=18)) + 
			y_ticks_off + 
			grid_off
	} else {
	fig_wrmean_withprior <- shared_params_df %>% 
		ggplot(aes(x=wrmean)) + 
			geom_histogram(aes(y=..density..), alpha=0.2, position="identity", bins=50) + 
			geom_density(adjust=2) + 
			geom_line(data=make_normal_prior_df(mean=prior_pars$wrmean_prior, sd=prior_pars$wrsd_prior, pmin=0.001, pmax=0.999, step=0.1), aes(x=x, y=density), col="black", linetype="dashed") + 
			scale_x_continuous(limits=c(0,NA)) + 
			labs(x="Mean clearance stage duration (days)", y="Density") + 
			theme_minimal() + 
			theme(legend.position="none", text=element_text(size=18)) + 
			y_ticks_off + 
			grid_off
		}


if(current_pars[["symptom_treatment"]]=="split"){
	fig_dpmean <- shared_params_df %>% 
		mutate(dpmean_symp=dpmeanS) %>% 
		mutate(dpmean_asymp=dpmeanA) %>% 
		select(dpmean_symp, dpmean_asymp) %>% 
		pivot_longer(everything()) %>% 
		ggplot(aes(x=global_pars[["lod"]]-value)) + 
			geom_histogram(aes(y=..density.., fill=name), alpha=0.2, position="identity", bins=50) + 
			geom_density(aes(col=name), adjust=2) + 
			# scale_x_continuous(limits=c(0,NA)) + 
			scale_color_manual(values=c("dpmean_symp"="red","dpmean_asymp"="blue")) + 
			scale_fill_manual(values=c("dpmean_symp"="red","dpmean_asymp"="blue")) + 
			labs(x="Mean peak Ct", y="Density") + 
			theme_minimal() + 
			theme(legend.position="none", text=element_text(size=18)) + 
			y_ticks_off + 
			grid_off
	} else{
	fig_dpmean <- shared_params_df %>% 
		ggplot(aes(x=global_pars[["lod"]]-dpmean)) + 
			geom_histogram(aes(y=..density..), alpha=0.2, position="identity", bins=50) + 
			geom_density(adjust=2) + 
			# scale_x_continuous(limits=c(0,NA)) + 
			labs(x="Mean peak Ct", y="Density") + 
			theme_minimal() + 
			theme(legend.position="none", text=element_text(size=18)) + 
			y_ticks_off + 
			grid_off
		}

if(current_pars[["symptom_treatment"]]=="split"){
	fig_wpmean <- shared_params_df %>% 
		mutate(wpmean_symp=wpmeanS) %>% 
		mutate(wpmean_asymp=wpmeanA) %>% 
		select(wpmean_symp, wpmean_asymp) %>% 
		pivot_longer(everything()) %>% 
		ggplot(aes(x=value)) + 
			geom_histogram(aes(y=..density.., fill=name), alpha=0.2, position="identity", bins=50) + 
			geom_density(aes(col=name), adjust=2) + 
			# scale_x_continuous(limits=c(0,NA)) + 
			scale_color_manual(values=c("wpmean_symp"="red","wpmean_asymp"="blue")) + 
			scale_fill_manual(values=c("wpmean_symp"="red","wpmean_asymp"="blue")) + 
			labs(x="Mean proliferation stage duration (days)", y="Density") + 
			theme_minimal() + 
			theme(legend.position="none", text=element_text(size=18)) + 
			y_ticks_off + 
			grid_off
	} else {
	fig_wpmean <- shared_params_df %>% 
		ggplot(aes(x=wpmean)) + 
			geom_histogram(aes(y=..density..), alpha=0.2, position="identity", bins=50) + 
			geom_density(adjust=2) + 
			# scale_x_continuous(limits=c(0,NA)) + 
			labs(x="Mean proliferation stage duration (days)", y="Density") + 
			theme_minimal() + 
			theme(legend.position="none", text=element_text(size=18)) + 
			y_ticks_off + 
			grid_off
		}

if(current_pars[["symptom_treatment"]]=="split"){
fig_wrmean <- shared_params_df %>% 
	mutate(wrmean_symp=wrmeanS) %>% 
	mutate(wrmean_asymp=wrmeanA) %>% 
	select(wrmean_symp, wrmean_asymp) %>% 
	pivot_longer(everything()) %>% 
	ggplot(aes(x=value)) + 
		geom_histogram(aes(y=..density.., fill=name), alpha=0.2, position="identity", bins=50) + 
		geom_density(aes(col=name), adjust=2) + 
		# scale_x_continuous(limits=c(0,NA)) + 
		scale_color_manual(values=c("wrmean_symp"="red","wrmean_asymp"="blue")) + 
		scale_fill_manual(values=c("wrmean_symp"="red","wrmean_asymp"="blue")) + 
		labs(x="Mean clearance stage duration (days)", y="Density") + 
		theme_minimal() + 
		theme(legend.position="none", text=element_text(size=18)) + 
		y_ticks_off + 
		grid_off
	} else {
	fig_wrmean <- shared_params_df %>% 
		ggplot(aes(x=wrmean)) + 
			geom_histogram(aes(y=..density..), alpha=0.2, position="identity", bins=50) + 
			geom_density(adjust=2) + 
			# scale_x_continuous(limits=c(0,NA)) + 
			labs(x="Mean clearance stage duration (days)", y="Density") + 
			theme_minimal() + 
			theme(legend.position="none", text=element_text(size=18)) + 
			y_ticks_off + 
			grid_off
	}

if(current_pars[["symptom_treatment"]]=="split"){
	fig_infdurmean <- shared_params_df %>% 
		mutate(infdurmean_symp=wpmeanS+wrmeanS) %>% 
		mutate(infdurmean_asymp=wpmeanA+wrmeanA) %>% 
		select(infdurmean_symp, infdurmean_asymp) %>% 
		pivot_longer(everything()) %>% 
		ggplot(aes(x=value)) + 
			geom_histogram(aes(y=..density.., fill=name), alpha=0.2, position="identity", bins=50) + 
			geom_density(aes(col=name), adjust=2) + 
			# scale_x_continuous(limits=c(0,NA)) + 
			scale_color_manual(values=c("infdurmean_symp"="red","infdurmean_asymp"="blue")) + 
			scale_fill_manual(values=c("infdurmean_symp"="red","infdurmean_asymp"="blue")) + 
			labs(x="Mean acute infection duration (days)", y="Density") + 
			theme_minimal() + 
			theme(legend.position="none", text=element_text(size=18)) + 
			y_ticks_off + 
			grid_off
		} else {
		fig_infdurmean <- shared_params_df %>% 
			mutate(infdurmean=wpmean+wrmean) %>% 
			select(infdurmean) %>% 
			pivot_longer(everything()) %>% 
			ggplot(aes(x=value)) + 
				geom_histogram(aes(y=..density..), alpha=0.2, position="identity", bins=50) + 
				geom_density( adjust=2) + 
				# scale_x_continuous(limits=c(0,NA)) + 
				labs(x="Mean acute infection duration (days)", y="Density") + 
				theme_minimal() + 
				theme(legend.position="none", text=element_text(size=18)) + 
				y_ticks_off + 
				grid_off
		}


# Overall posterior peak Ct density:
fig_peak_ct_overall <- with(as.list(c(global_pars, prior_pars)),{
	params_df %>% 
	ggplot(aes(x=lod-dp)) + 
		# Plot posterior density: 
		kernelplotfun(3) + 
		# Plot best estimate:  
		peak_ct_fit_curve + 
		scale_x_continuous(limits=c(0,40)) + 
		labs(x="Peak Ct", y="Density") + 
		y_ticks_off + 
		grid_off + 
		theme(text=element_text(size=18))
	})
# ggsave(figspot("peak_ct_overall.pdf"), fig_peak_ct_overall, width=8, height=3)

# Overall posterior onset-to-peak density:
fig_t_onset_to_peak_overall <- with(as.list(prior_pars),{params_df %>% 
	ggplot(aes(x=wp)) + 
		# Plot posterior density: 
		kernelplotfun(3) + 
		# Plot best estimate:  
		wp_gamma_fit_curve + 
		scale_x_continuous(limits=c(0,wpmax)) + 
		labs(x="Time from onset to peak (days)", y="Density") + 
		y_ticks_off + 
		grid_off + 
		theme(text=element_text(size=18))
	})
# ggsave(figspot("t_onset_to_peak_overall.pdf"), fig_t_onset_to_peak_overall, width=8, height=3)

# Overall posterior peak-to-recovery density:
fig_t_peak_to_recovery_overall <- with(as.list(prior_pars),{params_df %>% 
	ggplot(aes(x=wr)) + 
		# Plot posterior density: 
		kernelplotfun(3) + 
		# Plot prior: 
		# wr_prior_curve +
		# Plot best estimate:  
		wr_gamma_fit_curve + 
		scale_x_continuous(limits=c(0,wrmax)) + 
		labs(x="Time from peak to recovery (days)", y="Density") + 
		y_ticks_off + 
		grid_off + 
		theme(text=element_text(size=18))
	})
# ggsave(figspot("t_peak_to_recovery_overall.pdf"), fig_t_peak_to_recovery_overall, width=8, height=3)

fig_infection_duration_overall <- with(as.list(prior_pars),{params_df %>% 
	ggplot(aes(x=wr+wp)) + 
		# Plot posterior density: 
		kernelplotfun(5) + 
		# Plot best estimate:
		infection_duration_gamma_fit_curve + 
		scale_x_continuous(limits=c(0,wrmax+wpmax)) + 
		labs(x="Duration of infection (days)", y="Density") + 
		y_ticks_off + 
		grid_off + 
		theme(text=element_text(size=18))
	})
# ggsave(figspot("infection_duration_overall.pdf"), fig_infection_duration_overall, width=8, height=3)

# See if anything separates for the people who do/don't have symptoms: 

fig_peak_symptomatic <- with(as.list(global_pars),{params_df %>% 
	mutate(peak=lod-dp) %>% 
	plot_symptomatic_hists("peak", adj=5) + 
		labs(x="Peak Ct", y="Density") + 
		theme(text=element_text(size=18)) + 
		y_ticks_off + 
		grid_off 
		})
# ggsave(figspot("peak_symptomatic.pdf"), fig_peak_symptomatic, width=8, height=3)


fig_wp_symptomatic <- params_df %>% 
	plot_symptomatic_hists("wp", adj=5) + 
		labs(x="Time from onset to peak (days)", y="Density") + 
		theme(text=element_text(size=18)) + 
		y_ticks_off + 
		grid_off 
# ggsave(figspot("wp_symptomatic.pdf"), fig_wp_symptomatic, width=8, height=3)

fig_wr_symptomatic <- params_df %>% 
	plot_symptomatic_hists("wr", adj=5)  + 
		labs(x="Time from peak to recovery (days)", y="Density") + 
		theme(text=element_text(size=18)) + 
		y_ticks_off + 
		grid_off 
# ggsave(figspot("wr_symptomatic.pdf"), fig_wr_symptomatic, width=8, height=3)

fig_infection_duration_symptomatic <- with(as.list(global_pars),{params_df %>% 
	mutate(duration=wp+wr) %>% 
	plot_symptomatic_hists("duration", adj=5) + 
		labs(x="Duration of Infection (days)", y="Density") + 
		theme(text=element_text(size=18)) + 
		y_ticks_off + 
		grid_off 
		})
# ggsave(figspot("infection_duration_symptomatic.pdf"), fig_infection_duration_symptomatic, width=8, height=3)

# make the PPV figures: 
ct_probs_ma <- make_ct_probs_ma(ct_dat_clean, params_df, global_pars)
paired_ct_probs_ma <- make_paired_ct_probs_ma(ct_dat_clean, params_df, global_pars)

fig_infection_ppv_ma <- ct_probs_ma %>% 
	mutate(Ct.Change="Raw") %>%
	bind_rows(paired_ct_probs_ma) %>% 
	mutate(Infection=Onset+Resolution) %>%
	select(-Onset,-Resolution,-Chatter) %>%
	mutate(Infection_Upr=make_wald_v(Overall_n, Infection, 0.1, "upr")) %>%
	mutate(Infection_Lwr=make_wald_v(Overall_n, Infection, 0.1, "lwr")) %>%
	ggplot(aes(x=CT.Bin, y=Infection, col=Ct.Change)) + 
		geom_point() + 
		# geom_line() + 
		geom_line(stat="smooth", method="loess", span=0.5) + 
		geom_errorbar(aes(ymin=Infection_Lwr, ymax=Infection_Upr), width=0.2, alpha=0.4) + 
		scale_x_reverse() + 
		scale_y_continuous(breaks=seq(0,1,0.2), limits=c(-0.1,1.1)) + 
		scale_color_manual(values=c("Raw"="black","Increase"="blue","Decrease"="red")) + 
		theme_minimal() + 
		theme(text=element_text(size=18)) + 
		labs(x="Mean Ct for 5-unit sliding window", y="Probability of acute infection")

fig_infection_ppv_ma_raw <- ct_probs_ma %>% 
	mutate(Ct.Change="Raw") %>%
	mutate(Infection=Onset+Resolution) %>%
	select(-Onset,-Resolution,-Chatter) %>%
	mutate(Infection_Upr=make_wald_v(Overall_n, Infection, 0.1, "upr")) %>%
	mutate(Infection_Lwr=make_wald_v(Overall_n, Infection, 0.1, "lwr")) %>%
	ggplot(aes(x=CT.Bin, y=Infection)) + 
		geom_point() + 
		# geom_line() + 
		geom_line(stat="smooth", method="loess", span=0.5) + 
		geom_errorbar(aes(ymin=Infection_Lwr, ymax=Infection_Upr), width=0.2, alpha=0.4) + 
		scale_x_reverse() + 
		scale_y_continuous(breaks=seq(0,1,0.2), limits=c(-0.1,1.1)) + 
		scale_color_manual(values=c("Raw"="black","Increase"="blue","Decrease"="red")) + 
		theme_minimal() + 
		theme(text=element_text(size=18)) + 
		labs(x="Mean Ct for 5-unit sliding window", y="Probability of acute infection")

fig_infection_ppv_ma_paired <- paired_ct_probs_ma %>% 
	mutate(Infection=Onset+Resolution) %>%
	select(-Onset,-Resolution,-Chatter) %>%
	mutate(Infection_Upr=make_wald_v(Overall_n, Infection, 0.1, "upr")) %>%
	mutate(Infection_Lwr=make_wald_v(Overall_n, Infection, 0.1, "lwr")) %>%
	ggplot(aes(x=CT.Bin, y=Infection, col=Ct.Change)) + 
		geom_point() + 
		# geom_line() + 
		geom_line(stat="smooth", method="loess", span=0.5) + 
		geom_errorbar(aes(ymin=Infection_Lwr, ymax=Infection_Upr), width=0.2, alpha=0.4) + 
		scale_x_reverse() + 
		scale_y_continuous(breaks=seq(0,1,0.2), limits=c(-0.1,1.1)) + 
		scale_color_manual(values=c("Raw"="black","Increase"="blue","Decrease"="red")) + 
		theme_minimal() + 
		theme(text=element_text(size=18), legend.position="none") + 
		labs(x="Mean Ct for 5-unit sliding window", y="Probability of acute infection")

fig_onset_given_infection_ppv_ma <- ct_probs_ma %>% 
	select(-Chatter) %>%
	mutate(Denom = Onset+Resolution) %>% 
	mutate(Onset = Onset / Denom) %>% 
	mutate(Resolution = Resolution / Denom) %>% 
	select(-Denom) %>%
	mutate(Ct.Change="Raw") %>%
	bind_rows(paired_ct_probs_ma %>% 
		select(-Chatter) %>%
		mutate(Denom=Onset+Resolution) %>% 
		mutate(Onset = Onset / Denom) %>% 
		mutate(Resolution = Resolution / Denom) %>% 
		select(-Denom)) %>% 
	mutate(Onset_Upr=make_wald_v(Overall_n, Onset, 0.1, "upr")) %>%
	mutate(Onset_Lwr=make_wald_v(Overall_n, Onset, 0.1, "lwr")) %>%
	ggplot(aes(x=CT.Bin, y=Onset, col=Ct.Change)) + 
		geom_point() + 
		# geom_line() + 
		geom_line(stat="smooth", method="loess", span=1) + 
		geom_errorbar(aes(ymin=Onset_Lwr, ymax=Onset_Upr), width=0.2, alpha=0.4) + 
		scale_x_reverse() + 
		scale_y_continuous(breaks=seq(0,1,0.2), limits=c(-0.1,1.1)) + 
		scale_color_manual(values=c("Raw"="black","Increase"="blue","Decrease"="red")) + 
		theme_minimal() + 
		theme(text=element_text(size=18)) + 
		labs(x="Mean Ct for 5-unit sliding window", y="Probability of proliferation stage\ngiven acute infection")

fig_onset_given_infection_ppv_ma_raw <- ct_probs_ma %>% 
	select(-Chatter) %>%
	mutate(Denom = Onset+Resolution) %>% 
	mutate(Onset = Onset / Denom) %>% 
	mutate(Resolution = Resolution / Denom) %>% 
	select(-Denom) %>%
	mutate(Ct.Change="Raw") %>%
	mutate(Onset_Upr=make_wald_v(Overall_n, Onset, 0.1, "upr")) %>%
	mutate(Onset_Lwr=make_wald_v(Overall_n, Onset, 0.1, "lwr")) %>%
	ggplot(aes(x=CT.Bin, y=Onset)) + 
		geom_point() + 
		# geom_line() + 
		geom_line(stat="smooth", method="loess", span=1) + 
		geom_errorbar(aes(ymin=Onset_Lwr, ymax=Onset_Upr), width=0.2, alpha=0.4) + 
		scale_x_reverse() + 
		scale_y_continuous(breaks=seq(0,1,0.2), limits=c(-0.1,1.1)) + 
		scale_color_manual(values=c("Raw"="black","Increase"="blue","Decrease"="red")) + 
		theme_minimal() + 
		theme(text=element_text(size=18)) + 
		labs(x="Mean Ct for 5-unit sliding window", y="Probability of proliferation stage\ngiven acute infection")

fig_onset_given_infection_ppv_ma_paired <- paired_ct_probs_ma %>% 
	select(-Chatter) %>%
	mutate(Denom = Onset+Resolution) %>% 
	mutate(Onset = Onset / Denom) %>% 
	mutate(Resolution = Resolution / Denom) %>% 
	select(-Denom) %>%
	mutate(Onset_Upr=make_wald_v(Overall_n, Onset, 0.1, "upr")) %>%
	mutate(Onset_Lwr=make_wald_v(Overall_n, Onset, 0.1, "lwr")) %>%
	ggplot(aes(x=CT.Bin, y=Onset, col=Ct.Change)) + 
		geom_point() + 
		# geom_line() + 
		geom_line(stat="smooth", method="loess", span=1) + 
		geom_errorbar(aes(ymin=Onset_Lwr, ymax=Onset_Upr), width=0.2, alpha=0.4) + 
		scale_x_reverse() + 
		scale_y_continuous(breaks=seq(0,1,0.2), limits=c(-0.1,1.1)) + 
		scale_color_manual(values=c("Raw"="black","Increase"="blue","Decrease"="red")) + 
		theme_minimal() + 
		theme(text=element_text(size=18), legend.position="none") + 
		labs(x="Mean Ct for 5-unit sliding window", y="Probability of proliferation stage\ngiven acute infection")

srise <- function(x, dp, wp){
	out <- dp*(1+x/wp)
	return(out)
}
sfall <- function(x, dp, wr){
	out <- dp*(1-x/wr)
	return(out)
}

make_sample_trajectory_symp <- function(shared_params_df, global_pars, siglevel=0.9){
	# For asymptomatic:
	with(as.list(global_pars),{

	wp_mean_A <- mean(shared_params_df$wpmeanA)
	wp_lwr_A <- quantile(shared_params_df$wpmeanA,(1-siglevel)/2)
	wp_upr_A <- quantile(shared_params_df$wpmeanA,1-(1-siglevel)/2)

	wr_mean_A <- mean(shared_params_df$wrmeanA)
	wr_lwr_A <- quantile(shared_params_df$wrmeanA,(1-siglevel)/2)
	wr_upr_A <- quantile(shared_params_df$wrmeanA,1-(1-siglevel)/2)

	dp_mean_A <- mean(shared_params_df$dpmeanA)
	dp_lwr_A <- quantile(shared_params_df$dpmeanA,(1-siglevel)/2)
	dp_upr_A <- quantile(shared_params_df$dpmeanA,1-(1-siglevel)/2)

	xvals_proliferation_A <- seq(from=-wp_upr_A, 0, length.out=500)
	xvals_clearance_A <- seq(from=0, wr_upr_A, length.out=500)

	yvals_upr_proliferation_A <- unlist(lapply(xvals_proliferation_A, 
	function(x) quantile(srise(x, shared_params_df$dpmeanA, shared_params_df$wpmeanA),0.9)))
	yvals_lwr_proliferation_A <- unlist(lapply(xvals_proliferation_A, 
	function(x) quantile(srise(x, shared_params_df$dpmeanA, shared_params_df$wpmeanA),0.1)))
	yvals_upr_clearance_A <- unlist(lapply(xvals_clearance_A, 
	function(x) quantile(sfall(x, shared_params_df$dpmeanA, shared_params_df$wrmeanA),0.9)))
	yvals_lwr_clearance_A <- unlist(lapply(xvals_clearance_A, 
	function(x) quantile(sfall(x, shared_params_df$dpmeanA, shared_params_df$wrmeanA),0.1)))

	# For symptomatic:
	wp_mean_S <- mean(shared_params_df$wpmeanS)
	wp_lwr_S <- quantile(shared_params_df$wpmeanS,(1-siglevel)/2)
	wp_upr_S <- quantile(shared_params_df$wpmeanS,1-(1-siglevel)/2)

	wr_mean_S <- mean(shared_params_df$wrmeanS)
	wr_lwr_S <- quantile(shared_params_df$wrmeanS,(1-siglevel)/2)
	wr_upr_S <- quantile(shared_params_df$wrmeanS,1-(1-siglevel)/2)

	dp_mean_S <- mean(shared_params_df$dpmeanS)
	dp_lwr_S <- quantile(shared_params_df$dpmeanS,(1-siglevel)/2)
	dp_upr_S <- quantile(shared_params_df$dpmeanS,1-(1-siglevel)/2)

	xvals_proliferation_S <- seq(from=-wp_upr_S, 0, length.out=500)
	xvals_clearance_S <- seq(from=0, wr_upr_S, length.out=500)

	yvals_upr_proliferation_S <- unlist(lapply(xvals_proliferation_S, 
	function(x) quantile(srise(x, shared_params_df$dpmeanS, shared_params_df$wpmeanS),0.9)))
	yvals_lwr_proliferation_S <- unlist(lapply(xvals_proliferation_S, 
	function(x) quantile(srise(x, shared_params_df$dpmeanS, shared_params_df$wpmeanS),0.1)))
	yvals_upr_clearance_S <- unlist(lapply(xvals_clearance_S, 
	function(x) quantile(sfall(x, shared_params_df$dpmeanS, shared_params_df$wrmeanS),0.9)))
	yvals_lwr_clearance_S <- unlist(lapply(xvals_clearance_S, 
	function(x) quantile(sfall(x, shared_params_df$dpmeanS, shared_params_df$wrmeanS),0.1)))


	out <- ggplot() + 
		geom_ribbon(
			data=tibble(
				xvals=xvals_proliferation_A, 
				yvals_lwr=yvals_lwr_proliferation_A, 
				yvals_upr=yvals_upr_proliferation_A),
			aes(x=xvals, ymin=lod-yvals_lwr, ymax=lod-yvals_upr), alpha=0.2, fill="blue") + 
		geom_segment(aes(x=-wp_mean_A,xend=0,y=lod,yend=lod-dp_mean_A),col="blue") + 
		geom_ribbon(
			data=tibble(
				xvals=xvals_proliferation_S, 
				yvals_lwr=yvals_lwr_proliferation_S, 
				yvals_upr=yvals_upr_proliferation_S),
			aes(x=xvals, ymin=lod-yvals_lwr, ymax=lod-yvals_upr), alpha=0.2, fill="red") + 
		geom_segment(aes(x=-wp_mean_S,xend=0,y=lod,yend=lod-dp_mean_S),col="red") + 
		geom_ribbon(
			data=tibble(
				xvals=xvals_clearance_A, 
				yvals_lwr=yvals_lwr_clearance_A, 
				yvals_upr=yvals_upr_clearance_A),
			aes(x=xvals, ymin=lod-yvals_lwr, ymax=lod-yvals_upr), alpha=0.2, fill="blue") + 
		geom_segment(aes(x=0,xend=wr_mean_A,y=lod-dp_mean_A,yend=lod),col="blue") + 
		geom_ribbon(
			data=tibble(
				xvals=xvals_clearance_S, 
				yvals_lwr=yvals_lwr_clearance_S, 
				yvals_upr=yvals_upr_clearance_S),
			aes(x=xvals, ymin=lod-yvals_lwr, ymax=lod-yvals_upr), alpha=0.2, fill="red") + 
		geom_segment(aes(x=0,xend=wr_mean_S,y=lod-dp_mean_S,yend=lod),col="red") + 
		coord_cartesian(ylim=c(40,20), expand=FALSE) + 
		theme_minimal() + 
		labs(x="Days from peak", y="Ct") + 
		scale_y_reverse() + 
		theme(text=element_text(size=18))

	return(out)

	})

}


make_sample_trajectory <- function(shared_params_df, global_pars, siglevel=0.9){
	# For asymptomatic:
	with(as.list(global_pars),{

	wp_mean <- mean(shared_params_df$wpmean)
	wp_lwr <- quantile(shared_params_df$wpmean,(1-siglevel)/2)
	wp_upr <- quantile(shared_params_df$wpmean,1-(1-siglevel)/2)

	wr_mean <- mean(shared_params_df$wrmean)
	wr_lwr <- quantile(shared_params_df$wrmean,(1-siglevel)/2)
	wr_upr <- quantile(shared_params_df$wrmean,1-(1-siglevel)/2)

	dp_mean <- mean(shared_params_df$dpmean)
	dp_lwr <- quantile(shared_params_df$dpmean,(1-siglevel)/2)
	dp_upr <- quantile(shared_params_df$dpmean,1-(1-siglevel)/2)

	xvals_proliferation <- seq(from=-wp_upr, 0, length.out=500)
	xvals_clearance <- seq(from=0, wr_upr, length.out=500)

	yvals_upr_proliferation <- unlist(lapply(xvals_proliferation, 
	function(x) quantile(srise(x, shared_params_df$dpmean, shared_params_df$wpmean),0.9)))
	yvals_lwr_proliferation <- unlist(lapply(xvals_proliferation, 
	function(x) quantile(srise(x, shared_params_df$dpmean, shared_params_df$wpmean),0.1)))
	yvals_upr_clearance <- unlist(lapply(xvals_clearance, 
	function(x) quantile(sfall(x, shared_params_df$dpmean, shared_params_df$wrmean),0.9)))
	yvals_lwr_clearance <- unlist(lapply(xvals_clearance, 
	function(x) quantile(sfall(x, shared_params_df$dpmean, shared_params_df$wrmean),0.1)))


	out <- ggplot() + 
		geom_ribbon(
			data=tibble(
				xvals=xvals_proliferation, 
				yvals_lwr=yvals_lwr_proliferation, 
				yvals_upr=yvals_upr_proliferation),
			aes(x=xvals, ymin=lod-yvals_lwr, ymax=lod-yvals_upr), alpha=0.2, fill="grey") + 
		geom_segment(aes(x=-wp_mean,xend=0,y=lod,yend=lod-dp_mean),col="black") + 
		geom_ribbon(
			data=tibble(
				xvals=xvals_clearance, 
				yvals_lwr=yvals_lwr_clearance, 
				yvals_upr=yvals_upr_clearance),
			aes(x=xvals, ymin=lod-yvals_lwr, ymax=lod-yvals_upr), alpha=0.2, fill="grey") + 
		geom_segment(aes(x=0,xend=wr_mean,y=lod-dp_mean,yend=lod),col="black") + 
		coord_cartesian(ylim=c(40,20), expand=FALSE) + 
		theme_minimal() + 
		labs(x="Days from peak", y="Ct") + 
		scale_y_reverse() + 
		theme(text=element_text(size=18))

	return(out)

	})

}


if(current_pars[["symptom_treatment"]]=="split"){
	fig_ct_trajectory_inference <- make_sample_trajectory_symp(shared_params_df, global_pars)
	} else {
	fig_ct_trajectory_inference <- make_sample_trajectory(shared_params_df, global_pars)
	}

