# source('code/utilities/utils_effse.R')

library(tidyverse) 
library(purrr)

# Vectorize the uniform draw function:
runif_V = Vectorize(runif)

drawonsets <- function(trajectory_df, inf_ct, maxct, event_duration){
	out <- trajectory_df %>% 
		mutate(to_lwr = (wr/dp)*(maxct-dp-inf_ct) - wp ) %>% 
		mutate(to_upr = -(wp/dp)*(maxct-inf_ct) + event_duration) %>%
		mutate(to=runif_V(1,to_lwr,to_upr)) %>%
		select(-to_lwr, -to_upr) %>%
		select(to, wp, wr, dp)
}

drawonsets_infected <- function(trajectory_df, event_duration){
	out <- trajectory_df %>% 
		mutate(to=runif_V(1,-(wp+wr),event_duration)) %>%
		select(to, wp, wr, dp)
}


getct <- function(trajectory_df, t, maxct){

	out <- trajectory_df %>% 
		mutate(t=t) %>% 
		mutate(ct=case_when(
			(t>to) & (t<=(to+wp)) ~ maxct - (dp/wp)*(t-to),
			(t>(to+wp)) & (t<=(to+wp+wr)) ~ maxct - (dp-(dp/wr)*(t-(to+wp))),
			TRUE ~ maxct
			))

	return(out) 

}

runtest <- function(ct_df, lod, se){
	out <- ct_df %>% 
		mutate(falseneg=rbinom(n(),1,1-se)) %>% 
		mutate(screened=case_when(
			ct<lod & falseneg==0 ~ 1,
			TRUE ~ 0
		)) %>%
		select(-falseneg)

	return(out)
}

# In general, will want to put in a slice of params_df rather than the full thing.
get_effective_sensitivity <- function(t_test, lod, se, inf_ct, maxct, params_df, event_duration=3/24){

	trajectory_df <- params_df %>% 
		# sample_n(1000) %>%
		select(wp, wr, dp) %>%
		drawonsets(inf_ct, maxct, event_duration)

	out <- trajectory_df %>% 
		getct(t_test, maxct) %>% 
		runtest(lod=lod, se=se) %>%
		summarise(sensitivity=sum(screened)/n()) %>% 
		pull(sensitivity)

	return(out)

}

get_effective_sensitivity_infected <- function(t_test, lod, se, inf_ct, maxct, params_df, event_duration=3/24){

	trajectory_df <- params_df %>% 
		# sample_n(1000) %>%
		select(wp, wr, dp) %>%
		drawonsets_infected(event_duration)

	out <- trajectory_df %>% 
		getct(t_test, maxct) %>% 
		runtest(lod=lod, se=se) %>%
		summarise(sensitivity=sum(screened)/n()) %>% 
		pull(sensitivity)

	return(out)

}

# In general, will want to put in a slice of params_df rather than the full thing.
# get_n_infectious <- function(t_test, lod, se, inf_ct, maxct, params_df, pop_pars, event_duration=3/24, ngames=1000, siglevel=0.9){

# 	with(as.list(c(pop_pars)),{

# 		n_infectious_baseline <- rbinom(ngames, n_attendees, prev)
# 		eff_se <- get_effective_sensitivity(t_test, lod, se, inf_ct, maxct, params_df, event_duration)
# 		n_infectious <- unlist(lapply(n_infectious_baseline, function(x) rbinom(1, x, 1-eff_se)))
# 		out <- tibble(
# 			mean=mean(n_infectious),
# 			lwr=quantile(n_infectious, (1-siglevel)/2), 
# 			upr=quantile(n_infectious, 1-(1-siglevel)/2)) %>%
# 			pivot_longer(everything(), names_to="statistic", values_to="value") %>%
# 			mutate(t=t_test)

# 		return(out)

# 		})

# }

get_minct <- function(trajectory_df, maxct, event_duration){
	out <- trajectory_df %>% 
		mutate(minct = case_when(
			to+wp < 0 ~ maxct - dp*((to+wp)/wr + 1),
			to+wp >= 0 & to+wp <= event_duration ~ maxct - dp,
			to+wp > event_duration ~ maxct - dp*((event_duration - to - wp)/wp + 1)
			))
	return(out)
}

get_n_infectious <- function(t_test, lod, se, inf_ct, maxct, params_df, pop_pars, event_duration=3/24, ngames=1000, siglevel=0.9){

	with(as.list(c(pop_pars)),{

		n_infected_baseline <- rbinom(ngames, n_attendees, prev)
		out <- params_df %>%
			sample_n(sum(n_infected_baseline), replace=TRUE) %>% 
			drawonsets_infected(event_duration=event_duration) %>% 
			mutate(game=unlist(map2(1:ngames, n_infected_baseline, rep))) %>%
			get_minct(maxct=maxct, event_duration=event_duration) %>% 
			mutate(infectious=case_when(minct < inf_ct~1, TRUE~0)) %>% 
			getct(t_test, maxct) %>% 
			runtest(lod=lod, se=se) %>% 
			mutate(escaped=case_when(infectious==1 & screened==0~1, TRUE~0)) %>%
			group_by(game) %>% 
			summarise(
				n_infectious_baseline=sum(infectious),
				n_infectious=sum(escaped)) %>%
			summarise(mean=mean(n_infectious),
				lwr=quantile(n_infectious, (1-siglevel)/2), 
				upr=quantile(n_infectious, 1-(1-siglevel)/2),
				mean_baseline=mean(n_infectious_baseline),
				lwr_baseline=quantile(n_infectious_baseline, (1-siglevel)/2), 
				upr_baseline=quantile(n_infectious_baseline, 1-(1-siglevel)/2)) %>%
			pivot_longer(everything(), names_to="statistic", values_to="value") %>%
			mutate(t=t_test)

		# lapply(n_infected_baseline, rtrajectory_infected, trajectory_pars=trajectory_pars, event_duration=event_duration) %>%
		# 	bind_rows(.id="game") %>% 
		# 	mutate(game=as.numeric(game)) %>% 
		# 	# Calculate each person's min ct at the event: 
		# 	get_minct(trajectory_pars=trajectory_pars, event_duration=event_duration) %>% 
		# 	# Tag if the person is infectious at the event: 
		# 	mutate(infectious=case_when(minct < inf_ct~1, TRUE~0)) %>% 
		# 	# Get ct at time of test
		# 	getct(t_test, trajectory_pars) %>% 
		# 	# run the test
		# 	runtest(lod=lod, se=se) %>% 
		# 	# See how many infectious escape detection
		# 	mutate(escaped=case_when(infectious==1 & screened==0~1, TRUE~0)) %>%
		# 	# summarize results
		# 	group_by(game) %>% 
		# 	summarise(
		# 		n_infectious_baseline=sum(infectious),
		# 		n_infectious=sum(escaped)) %>%
		# 	summarise(mean=mean(n_infectious),
		# 		lwr=quantile(n_infectious, (1-siglevel)/2), 
		# 		upr=quantile(n_infectious, 1-(1-siglevel)/2),
		# 		mean_baseline=mean(n_infectious_baseline),
		# 		lwr_baseline=quantile(n_infectious_baseline, (1-siglevel)/2), 
		# 		upr_baseline=quantile(n_infectious_baseline, 1-(1-siglevel)/2)) %>%
		# 	pivot_longer(everything(), names_to="statistic", values_to="value") %>%
		# 	mutate(t=t_test)

		# eff_se <- get_effective_sensitivity_infected(t_test, lod, se, inf_ct, maxct, params_df, event_duration)
		# n_infectious <- unlist(lapply(n_infectious_baseline, function(x) rbinom(1, x, 1-eff_se)))
		# out <- tibble(
		# 	mean=mean(n_infectious),
		# 	lwr=quantile(n_infectious, (1-siglevel)/2), 
		# 	upr=quantile(n_infectious, 1-(1-siglevel)/2)) %>%
		# 	pivot_longer(everything(), names_to="statistic", values_to="value") %>%
		# 	mutate(t=t_test)

		return(out)

		})

}