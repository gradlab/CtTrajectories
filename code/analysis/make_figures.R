# =============================================================================
# Import: 
# =============================================================================

library(tidyverse) 
library(scales)
source('code/utilities/utils_analysis.R')
source('code/utilities/utils_effse.R')

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

jitterfactor <- 0.01
meanvalsindiv <- params_df %>% 
	mutate(infdur=wp+wr) %>%
	group_by(id) %>% 
	summarise(tp=mean(tp), dp=mean(dp), wp=mean(wp), wr=mean(wr), infdur=mean(infdur), symptomatic=first(symptomatic))

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
		labs(x="Days since min. observed Ct", y="Ct") + 
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
		mutate(dpmean_symp=dpmeanS_trans) %>% 
		mutate(dpmean_asymp=dpmeanA_trans) %>% 
		select(dpmean_symp, dpmean_asymp) %>% 
		pivot_longer(everything()) %>% 
		ggplot(aes(x=global_pars[["lod"]]-value)) + 
			geom_histogram(aes(y=..density.., fill=name), alpha=0.2, position="identity", bins=50) + 
			geom_density(aes(col=name), adjust=2) + 
			geom_line(data=make_normal_prior_df(mean=prior_pars$dpmean_prior, sd=prior_pars$dpsd_prior, pmin=0.001, pmax=0.999, step=0.1), aes(x=x, y=density), col="black", linetype="dashed") + 
			scale_x_continuous(limits=c(0,NA)) + 
			scale_color_manual(values=c("dpmean_symp"="red","dpmean_asymp"="blue","1"="red","0"="blue")) + 
			scale_fill_manual(values=c("dpmean_symp"="red","dpmean_asymp"="blue")) + 
			geom_jitter(data=meanvalsindiv, aes(x=global_pars[["lod"]]-dp, y=jitterfactor*symptomatic, col=factor(symptomatic)), alpha=0.4, width=0, height=0.003) + 
			labs(x="Mean peak Ct", y="Density") + 
			theme_minimal() + 
			theme(legend.position="none", text=element_text(size=18)) + 
			y_ticks_off + 
			grid_off
} else {
	fig_dpmean_withprior <- shared_params_df %>% 
		ggplot(aes(x=global_pars[["lod"]]-dpmean_trans)) + 
			geom_histogram(aes(y=..density..), alpha=0.2, position="identity", bins=50) + 
			geom_density(adjust=2) + 
			geom_line(data=make_normal_prior_df(mean=prior_pars$dpmean_prior, sd=prior_pars$dpsd_prior, pmin=0.001, pmax=0.999, step=0.1), aes(x=x, y=density), col="black", linetype="dashed") + 
			scale_x_continuous(limits=c(0,NA)) + 
			geom_jitter(data=meanvalsindiv, aes(x=global_pars[["lod"]]-dp, y=0), alpha=0.4, width=0, height=0.003)  + 
			labs(x="Mean peak Ct", y="Density") + 
			theme_minimal() + 
			theme(legend.position="none", text=element_text(size=18)) + 
			y_ticks_off + 
			grid_off
}


if(current_pars[["symptom_treatment"]]=="split"){
	fig_wpmean_withprior <- shared_params_df %>% 
		mutate(wpmean_symp=wpmeanS_trans) %>% 
		mutate(wpmean_asymp=wpmeanA_trans) %>% 
		select(wpmean_symp, wpmean_asymp) %>% 
		pivot_longer(everything()) %>% 
		ggplot(aes(x=value)) + 
			geom_histogram(aes(y=..density.., fill=name), alpha=0.2, position="identity", bins=50) + 
			geom_density(aes(col=name), adjust=2) + 
			geom_line(data=make_normal_prior_df(mean=prior_pars$wpmean_prior, sd=prior_pars$wpsd_prior, pmin=0.001, pmax=0.999, step=0.1), aes(x=x, y=density), col="black", linetype="dashed") + 
			scale_x_continuous(limits=c(0,NA)) + 
			scale_color_manual(values=c("wpmean_symp"="red","wpmean_asymp"="blue","1"="red","0"="blue")) + 
			scale_fill_manual(values=c("wpmean_symp"="red","wpmean_asymp"="blue")) + 
			geom_jitter(data=meanvalsindiv, aes(x=wp, y=jitterfactor*symptomatic, col=factor(symptomatic)), alpha=0.4, width=0, height=0.003)  + 
			labs(x="Mean proliferation stage duration (days)", y="Density") + 
			theme_minimal() + 
			theme(legend.position="none", text=element_text(size=18)) + 
			y_ticks_off + 
			grid_off
	} else {
	fig_wpmean_withprior <- shared_params_df %>% 
		ggplot(aes(x=wpmean_trans)) + 
			geom_histogram(aes(y=..density..), alpha=0.2, position="identity", bins=50) + 
			geom_density(adjust=2) + 
			geom_line(data=make_normal_prior_df(mean=prior_pars$wpmean_prior, sd=prior_pars$wpsd_prior, pmin=0.001, pmax=0.999, step=0.1), aes(x=x, y=density), col="black", linetype="dashed") + 
			scale_x_continuous(limits=c(0,NA)) + 
			geom_jitter(data=meanvalsindiv, aes(x=wp, y=jitterfactor*symptomatic), alpha=0.4, width=0, height=0.003)  +
			labs(x="Mean proliferation stage duration (days)", y="Density") + 
			theme_minimal() + 
			theme(legend.position="none", text=element_text(size=18)) + 
			y_ticks_off + 
			grid_off
	}


if(current_pars[["symptom_treatment"]]=="split"){
	fig_wrmean_withprior <- shared_params_df %>% 
		mutate(wrmean_symp=wrmeanS_trans) %>% 
		mutate(wrmean_asymp=wrmeanA_trans) %>% 
		select(wrmean_symp, wrmean_asymp) %>% 
		pivot_longer(everything()) %>% 
		ggplot(aes(x=value)) + 
			geom_histogram(aes(y=..density.., fill=name), alpha=0.2, position="identity", bins=50) + 
			geom_density(aes(col=name), adjust=2) + 
			geom_line(data=make_normal_prior_df(mean=prior_pars$wrmean_prior, sd=prior_pars$wrsd_prior, pmin=0.001, pmax=0.999, step=0.1), aes(x=x, y=density), col="black", linetype="dashed") + 
			scale_x_continuous(limits=c(0,NA)) + 
			scale_color_manual(values=c("wrmean_symp"="red","wrmean_asymp"="blue","1"="red","0"="blue")) + 
			scale_fill_manual(values=c("wrmean_symp"="red","wrmean_asymp"="blue")) + 
			geom_jitter(data=meanvalsindiv, aes(x=wr, y=jitterfactor*symptomatic, col=factor(symptomatic)), alpha=0.4, width=0, height=0.003)  + 
			labs(x="Mean clearance stage duration (days)", y="Density") + 
			theme_minimal() + 
			theme(legend.position="none", text=element_text(size=18)) + 
			y_ticks_off + 
			grid_off
	} else {
	fig_wrmean_withprior <- shared_params_df %>% 
		ggplot(aes(x=wrmean_trans)) + 
			geom_histogram(aes(y=..density..), alpha=0.2, position="identity", bins=50) + 
			geom_density(adjust=2) + 
			geom_line(data=make_normal_prior_df(mean=prior_pars$wrmean_prior, sd=prior_pars$wrsd_prior, pmin=0.001, pmax=0.999, step=0.1), aes(x=x, y=density), col="black", linetype="dashed") + 
			scale_x_continuous(limits=c(0,NA)) + 
			geom_jitter(data=meanvalsindiv, aes(x=wr, y=jitterfactor*symptomatic), alpha=0.4, width=0, height=0.003)  + 
			labs(x="Mean clearance stage duration (days)", y="Density") + 
			theme_minimal() + 
			theme(legend.position="none", text=element_text(size=18)) + 
			y_ticks_off + 
			grid_off
		}


if(current_pars[["symptom_treatment"]]=="split"){
	fig_dpmean <- shared_params_df %>% 
		mutate(dpmean_symp=dpmeanS_trans) %>% 
		mutate(dpmean_asymp=dpmeanA_trans) %>% 
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
		ggplot(aes(x=global_pars[["lod"]]-dpmean_trans)) + 
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
		mutate(wpmean_symp=wpmeanS_trans) %>% 
		mutate(wpmean_asymp=wpmeanA_trans) %>% 
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
		ggplot(aes(x=wpmean_trans)) + 
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
	mutate(wrmean_symp=wrmeanS_trans) %>% 
	mutate(wrmean_asymp=wrmeanA_trans) %>% 
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
		ggplot(aes(x=wrmean_trans)) + 
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
		mutate(infdurmean_symp=wpmeanS_trans+wrmeanS_trans) %>% 
		mutate(infdurmean_asymp=wpmeanA_trans+wrmeanA_trans) %>% 
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
			mutate(infdurmean=wpmean_trans+wrmean_trans) %>% 
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



if(current_pars[["symptom_treatment"]]=="split"){
	fig_gemlmean <- shared_params_df %>% 
		mutate(dpmean_symp=dpmeanS_trans) %>% 
		mutate(dpmean_asymp=dpmeanA_trans) %>% 
		select(dpmean_symp, dpmean_asymp) %>% 
		pivot_longer(everything()) %>% 
		ggplot(aes(x=10^convert_Ct_logGEML(global_pars[["lod"]]-value))) + 
			geom_histogram(aes(y=..density.., fill=name), alpha=0.2, position="identity", bins=50) + 
			geom_density(aes(col=name), adjust=2) + 
			scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + 
			# scale_x_continuous(limits=c(0,NA)) + 
			scale_color_manual(values=c("dpmean_symp"="red","dpmean_asymp"="blue")) + 
			scale_fill_manual(values=c("dpmean_symp"="red","dpmean_asymp"="blue")) + 
			labs(x="Genome equivalents/ml", y="Density") + 
			theme_minimal() + 
			theme(legend.position="none", text=element_text(size=18)) + 
			y_ticks_off + 
			grid_off
	} else{
	fig_gemlmean <- shared_params_df %>% 
		ggplot(aes(x=10^convert_Ct_logGEML(global_pars[["lod"]]-dpmean_trans))) + 
			geom_histogram(aes(y=..density..), alpha=0.2, position="identity", bins=50) + 
			geom_density(adjust=2) + 
			scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) +
			# scale_x_continuous(limits=c(0,NA)) + 
			labs(x="Genome equivalents/ml", y="Density") + 
			theme_minimal() + 
			theme(legend.position="none", text=element_text(size=18)) + 
			y_ticks_off + 
			grid_off
		}



if(current_pars[["symptom_treatment"]]=="split"){
	fig_apmean <- shared_params_df %>% 
		mutate(apmean_symp=dpmeanS_trans/wpmeanS_trans) %>% 
		mutate(apmean_asymp=dpmeanA_trans/wpmeanA_trans) %>% 
		select(apmean_symp, apmean_asymp) %>% 
		pivot_longer(everything()) %>% 
		ggplot(aes(x=value)) + 
			geom_histogram(aes(y=..density.., fill=name), alpha=0.2, position="identity", bins=50) + 
			geom_density(aes(col=name), adjust=2) + 
			# scale_x_continuous(limits=c(0,NA)) + 
			scale_color_manual(values=c("apmean_symp"="red","apmean_asymp"="blue")) + 
			scale_fill_manual(values=c("apmean_symp"="red","apmean_asymp"="blue")) + 
			labs(x="Mean proliferation stage slope (Ct/day)", y="Density") + 
			theme_minimal() + 
			theme(legend.position="none", text=element_text(size=18)) + 
			y_ticks_off + 
			grid_off
	} else {
	fig_apmean <- shared_params_df %>% 
		ggplot(aes(x=dpmean_trans/wpmean_trans)) + 
			geom_histogram(aes(y=..density..), alpha=0.2, position="identity", bins=50) + 
			geom_density(adjust=2) + 
			# scale_x_continuous(limits=c(0,NA)) + 
			labs(x="Mean proliferation stage slope (Ct/day)", y="Density") + 
			theme_minimal() + 
			theme(legend.position="none", text=element_text(size=18)) + 
			y_ticks_off + 
			grid_off
		}


if(current_pars[["symptom_treatment"]]=="split"){
	fig_armean <- shared_params_df %>% 
		mutate(armean_symp=-dpmeanS_trans/wrmeanS_trans) %>% 
		mutate(armean_asymp=-dpmeanA_trans/wrmeanA_trans) %>% 
		select(armean_symp, armean_asymp) %>% 
		pivot_longer(everything()) %>% 
		ggplot(aes(x=value)) + 
			geom_histogram(aes(y=..density.., fill=name), alpha=0.2, position="identity", bins=50) + 
			geom_density(aes(col=name), adjust=2) + 
			# scale_x_continuous(limits=c(0,NA)) + 
			scale_color_manual(values=c("armean_symp"="red","armean_asymp"="blue")) + 
			scale_fill_manual(values=c("armean_symp"="red","armean_asymp"="blue")) + 
			labs(x="Mean clearance stage slope (Ct/day)", y="Density") + 
			theme_minimal() + 
			theme(legend.position="none", text=element_text(size=18)) + 
			y_ticks_off + 
			grid_off
	} else {
	fig_armean <- shared_params_df %>% 
		ggplot(aes(x=-dpmean_trans/wrmean_trans)) + 
			geom_histogram(aes(y=..density..), alpha=0.2, position="identity", bins=50) + 
			geom_density(adjust=2) + 
			# scale_x_continuous(limits=c(0,NA)) + 
			labs(x="Mean clearance stage slope (Ct/day)", y="Density") + 
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


if(current_pars[["symptom_treatment"]]=="split"){
	fig_ct_trajectory_inference <- make_sample_trajectory_symp(shared_params_df, global_pars)
	} else {
	fig_ct_trajectory_inference <- make_sample_trajectory(shared_params_df, global_pars)
	}


if(current_pars[["symptom_treatment"]]=="split"){
	fig_ge_trajectory_inference <- make_sample_trajectory_symp(shared_params_df, global_pars, ge=TRUE)
	} else {
	fig_ge_trajectory_inference <- make_sample_trajectory(shared_params_df, global_pars, ge=TRUE)
	}



# =============================================================================
# Initialize effective sensitivity parameters
# =============================================================================

t_test <- seq(-3, 0, 1/24) 
maxct <- 40
inf_ct <- 30

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