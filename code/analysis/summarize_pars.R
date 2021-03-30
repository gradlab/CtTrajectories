indiv_params_df <- make_indiv_params_df(params, c("tp","dp","wp","wr"), n_indiv) %>% 
	rename(id_clean=id) %>% 
	left_join(id_map, by="id_clean") %>%
	left_join(symptom_map, by="id")
if(current_pars[["symptom_treatment"]]=="split"){
	shared_params_df <- make_shared_params_df(params, c("dpmeanA","wpmeanA","wrmeanA","dpmeanS","wpmeanS","wrmeanS","dpsd","wpsd","wrsd")) %>% 
		mutate(dpmeanA_trans=truncnormmean(dpmeanA,dpsd,0,global_pars[["lod"]])) %>%
		mutate(dpmeanS_trans=truncnormmean(dpmeanS,dpsd,0,global_pars[["lod"]])) %>%
		mutate(wpmeanA_trans=truncnormmean(wpmeanA,wpsd,0,prior_pars$wpmax)) %>%
		mutate(wpmeanS_trans=truncnormmean(wpmeanS,wpsd,0,prior_pars$wpmax)) %>%
		mutate(wrmeanA_trans=truncnormmean(wrmeanA,wrsd,0,prior_pars$wrmax)) %>%
		mutate(wrmeanS_trans=truncnormmean(wrmeanS,wrsd,0,prior_pars$wrmax))
	} else {
	shared_params_df <- make_shared_params_df(params, c("dpmean","wpmean","wrmean","dpsd","wpsd","wrsd"))  %>% 
		mutate(dpmean_trans=truncnormmean(dpmean,dpsd,0,global_pars[["lod"]])) %>%
		mutate(wpmean_trans=truncnormmean(wpmean,wpsd,0,prior_pars$wpmax)) %>%
		mutate(wrmean_trans=truncnormmean(wrmean,wrsd,0,prior_pars$wrmax)) 
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