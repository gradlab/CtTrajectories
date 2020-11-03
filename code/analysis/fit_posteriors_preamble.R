library(tidyverse)
library(lazyeval)
library(rstan) 
library(shinystan) 
library(purrr)
library(data.table)
options(mc.cores=parallel::detectCores())
source('code/utilities/utils_analysis.R')

# Store the number of people we've kept: 
n_indiv <- length(unique(ct_dat_refined$Person.ID))

# Store the ids of the people we've kept in the refined set: 
kept_ids <- data.frame(Person.ID=sort(unique(ct_dat_refined$Person.ID)), kept=1)

# Define a pared-down dataset for passing to Stan: 
indiv_data <- with(as.list(global_pars),{
	ct_dat_refined %>% 
	select(Person.ID, Person.ID.Clean, Date.Index, CT.Mean, Adjusted) %>%
	rename(id=Person.ID) %>%
	rename(id_clean=Person.ID.Clean) %>% 
	rename(t=Date.Index) %>%
	rename(y=CT.Mean) %>%
	rename(adjusted=Adjusted) %>%
	replace_na(list(y=lod, adjusted=0)) %>%
	trim_negatives(global_pars)
	})

# Useful dataframe for mapping official ID to Stan ID:
id_map <- ct_dat_refined %>% 
	group_by(Person.ID) %>%
	summarise(Person.ID.Clean=first(Person.ID.Clean)) %>% 
	select(Person.ID, Person.ID.Clean) %>%
	rename(id=Person.ID) %>% 
	rename(id_clean=Person.ID.Clean) %>%
	mutate(id_clean=as.character(id_clean))

# Useful dataframe for mapping official ID to symptoms
symptom_map <- ct_dat_refined %>% 
	mutate(Symptomatic=case_when(Symptomatic=="Yes"~1, TRUE~0)) %>% 
	group_by(Person.ID) %>%
	summarise(Symptomatic=max(Symptomatic)) %>%
	rename(id=Person.ID) %>%
	rename(symptomatic=Symptomatic)

# Append symptoms onto indiv_data:
indiv_data <- left_join(indiv_data, symptom_map, by="id")

# Generatae list of individuals with symptoms for passing to Stan:
symp <- indiv_data %>% 
	group_by(id) %>% 
	slice(1) %>% 
	select(id, symptomatic) %>% 
	arrange(id) %>% 
	pull(symptomatic)

# Store key prior parameters for passing to Stan: 
prior_pars <- with(as.list(current_pars),{
	list(
	symp=symp,
	tpsd=2,
	dpmean_prior=global_pars[["lod"]]/2,
	dpsd_prior=global_pars[["lod"]]/6,
	apmean_prior=5,
	apsd_prior=3,
	armean_prior=-5,
	arsd_prior=3,
	sigma_max=10,
	sigma_prior_scale=5,
	lambda=0.01,
	fpmean=1/log(10)  # so that 90% of mass is <1 and 99% is <2
	) 
	})

sensitivity_pars <- with(as.list(current_pars),{
	list(
		splitsymptoms=if_else(current_pars[["symptom_treatment"]]=="split",1,0)
	)
	})