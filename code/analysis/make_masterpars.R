library(tidyverse) 
# Define the set of sensitivities to run: 

masterpars <- list()

masterpars[[1]] <- c(symptom_treatment="split", 
	parametrization="waittime",
	toexclude=c(""),
	dpmean_prior=global_pars[["lod"]]/2,
	dpsd_prior=global_pars[["lod"]]/6,
	wpmax=14,
	wpmean_prior=14/2,
	wpsd_prior=14/6,
	wrmax=30,
	wrmean_prior=30/2,
	wrsd_prior=30/6,
	sigma_max=10,
	sigma_prior_scale=5,
	lambda=0.01)

masterpars[[2]] <- c(symptom_treatment="combined", 
	parametrization="waittime",
	toexclude=c(""),
	dpmean_prior=global_pars[["lod"]]/2,
	dpsd_prior=global_pars[["lod"]]/6,
	wpmax=14,
	wpmean_prior=14/2,
	wpsd_prior=14/6,
	wrmax=30,
	wrmean_prior=30/2,
	wrsd_prior=30/6,
	sigma_max=10,
	sigma_prior_scale=5,
	lambda=0.01)

# Base case but excluding person 3047
masterpars[[3]] <- c(symptom_treatment="split", 
	parametrization="waittime",
	toexclude=c(3047),
	dpmean_prior=global_pars[["lod"]]/2,
	dpsd_prior=global_pars[["lod"]]/6,
	wpmax=14,
	wpmean_prior=14/2,
	wpsd_prior=14/6,
	wrmax=30,
	wrmean_prior=30/2,
	wrsd_prior=30/6,
	sigma_max=10,
	sigma_prior_scale=5,
	lambda=0.01)

# Base case but with 95% PCR sensitivity
masterpars[[4]] <- c(symptom_treatment="split", 
	parametrization="waittime",
	toexclude=c(""),
	dpmean_prior=global_pars[["lod"]]/2,
	dpsd_prior=global_pars[["lod"]]/6,
	wpmax=14,
	wpmean_prior=14/2,
	wpsd_prior=14/6,
	wrmax=30,
	wrmean_prior=30/2,
	wrsd_prior=30/6,
	sigma_max=10,
	sigma_prior_scale=5,
	lambda=0.05)

# Base case but without the constraints on wp and wr
masterpars[[5]] <- c(symptom_treatment="split", 
	parametrization="waittime",
	toexclude=c(""),
	dpmean_prior=global_pars[["lod"]]/2,
	dpsd_prior=global_pars[["lod"]]/6,
	wpmax=100,
	wpmean_prior=14/2,
	wpsd_prior=14/6,
	wrmax=100,
	wrmean_prior=30/2,
	wrsd_prior=30/6,
	sigma_max=10,
	sigma_prior_scale=5,
	lambda=0.01)

# Base case but with low priors on wp and wr
masterpars[[6]] <- c(symptom_treatment="split", 
	parametrization="waittime",
	toexclude=c(""),
	dpmean_prior=global_pars[["lod"]]/2,
	dpsd_prior=global_pars[["lod"]]/6,
	wpmax=14,
	wpmean_prior=14/4,
	wpsd_prior=14/6,
	wrmax=30,
	wrmean_prior=30/4,
	wrsd_prior=30/6,
	sigma_max=10,
	sigma_prior_scale=5,
	lambda=0.01)

# Base case but with high priors on wp and wr
masterpars[[7]] <- c(symptom_treatment="split", 
	parametrization="waittime",
	toexclude=c(""),
	dpmean_prior=global_pars[["lod"]]/2,
	dpsd_prior=global_pars[["lod"]]/6,
	wpmax=14,
	wpmean_prior=3*14/4,
	wpsd_prior=14/6,
	wrmax=30,
	wrmean_prior=3*30/4,
	wrsd_prior=30/6,
	sigma_max=10,
	sigma_prior_scale=5,
	lambda=0.01)