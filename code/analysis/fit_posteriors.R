
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
	iter=500, chains=4)
# , control = list(adapt_delta=0.99)

fit_endq <- Sys.time()
print(paste0("Fit time: ",difftime(fit_endq, fit_startq, units="min")," mins"))

params <- rstan::extract(ct_fit)
source('code/analysis/summarize_pars.R')

# launch_shinystan_nonblocking(ct_fit)
	