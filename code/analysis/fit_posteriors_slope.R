
if(current_pars[["symptom_treatment"]]=="split"){
	ct_model <- stan_model("code/analysis/fit_posteriors_symptomatic.stan") 
} else {
	ct_model <- stan_model("code/analysis/fit_posteriors_combined.stan") 
}
fit_startq <- Sys.time()
ct_fit <- sampling(ct_model, 
	data=list(
		N=nrow(indiv_data), 
		n_id=length(unique(indiv_data$id_clean)),
		lod=as.list(global_pars)$lod, 
		id=indiv_data$id_clean,
		symp=as.list(prior_pars)$symp,
		t=indiv_data$t, 
		y=indiv_data$y, 
		tpsd=as.list(prior_pars)$tpsd,
		dpmean_prior=as.list(prior_pars)$dpmean_prior,
		dpsd_prior=as.list(prior_pars)$dpsd_prior,
		wpmax=as.list(prior_pars)$wpmax,
		wpmean_prior=as.list(prior_pars)$wpmean_prior,
		wpsd_prior=as.list(prior_pars)$wpsd_prior,
		wrmax=as.list(prior_pars)$wrmax,
		wrmean_prior=as.list(prior_pars)$wrmean_prior,
		wrsd_prior=as.list(prior_pars)$wrsd_prior,
		sigma_max=as.list(prior_pars)$sigma_max,
		sigma_prior_scale=as.list(prior_pars)$sigma_prior_scale,
		lambda=as.list(prior_pars)$lambda,
		fpmean=as.list(prior_pars)$fpmean,
		epsilon=(indiv_data$adjusted)*(as.list(global_pars)$adjusted_sd)), 
	iter=5000, chains=4, control = list(adapt_delta=0.99))

fit_endq <- Sys.time()
print(paste0("Fit time: ",difftime(fit_endq, fit_startq, units="min")," mins"))

params <- rstan::extract(ct_fit)
indiv_params_df <- make_indiv_params_df(params, c("tp","dp","wp","wr"), n_indiv) %>% 
	rename(id_clean=id) %>% 
	left_join(id_map, by="id_clean") %>%
	left_join(symptom_map, by="id")
if(current_pars[["symptom_treatment"]]=="split"){
	shared_params_df <- make_shared_params_df(params, c("dpmeanA","wpmeanA","wrmeanA","dpmeanS","wpmeanS","wrmeanS","dpsd","wpsd","wrsd")) 
	} else {
	shared_params_df <- make_shared_params_df(params, c("dpmean","wpmean","wrmean","dpsd","wpsd","wrsd")) 
	}

params_df <- indiv_params_df %>% 
	left_join(shared_params_df, by="iteration") %>% 
	select(-iteration) 

params_summary <- summary(ct_fit)$summary %>% 
	as.data.frame() %>% 
	setDT(keep.rownames=TRUE) %>% 
	as_tibble() %>% 
	rename(param=rn) %>%
	filter(
		grepl('dp',param) | 
		grepl('wp',param) | 
		grepl('wr',param) | 
		grepl('tp',param) | 
		grepl('sigma',param))

fit_summary <- get_sampler_params(ct_fit)
n_divergent <- 0
for(fsindex in 1:length(fit_summary)){
	fit_summary_chain <- as_tibble(fit_summary[[fsindex]])
	nits <- nrow(fit_summary_chain) 
	fit_summary_chain <- fit_summary_chain[floor(nits/2):nits,]
	n_divergent <- n_divergent + sum(fit_summary_chain$divergent__)
}

# launch_shinystan_nonblocking(ct_fit)
	