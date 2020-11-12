library(tidyverse) 
library(purrr)
library(truncnorm)

# Vectorize the uniform draw function:
runif_V = Vectorize(runif)

# Draw gamma-distributed values given a mean and sd; 
rgamma_musd <- function(n, mean, sd){
	shape <- (mean^2)/(sd^2)
	rate <- mean/(sd^2)
	out <- rgamma(n, shape=shape, rate=rate) 
	return(out) 
}

qgamma_musd <- function(p, mean, sd){
	shape <- (mean^2)/(sd^2)
	rate <- mean/(sd^2)
	out <- qgamma(p, shape=shape, rate=rate) 
	return(out) 
}

# Randomly sample trajectories (output is a dataframe with n rows)
rtrajectory <- function(n, trajectory_pars, event_duration){
	with(as.list(c(trajectory_pars)),{
 
		# Assign trajectory values to all individuals 
		wpvec <- rgamma_musd(n, wpmean, wpsd)
		wrvec <- rgamma_musd(n, wrmean, wrsd)
		dpvec <- rtruncnorm(n, a=maxct-inf_ct, b=maxct, dpmean, dpsd)

		trajectory_df <- tibble(
			wp=wpvec,
			wr=wrvec,
			dp=dpvec)

		trajectory_df <- trajectory_df %>% 
			mutate(to_lwr = (wr/dp)*(maxct-dp-inf_ct) - wp ) %>% 
			mutate(to_upr = -(wp/dp)*(maxct-inf_ct) + event_duration) %>%
			mutate(to=runif_V(1,to_lwr,to_upr)) %>%
			select(-to_lwr, -to_upr) %>%
			mutate(id=1:n()) %>%
			select(id, to, wp, wr, dp)
	
		return(trajectory_df)	
	})
}

getct <- function(trajectory_df, t, trajectory_pars){

	with(as.list(trajectory_pars),{

	out <- trajectory_df %>% 
		mutate(t=t) %>% 
		mutate(ct=case_when(
			(t>to) & (t<=(to+wp)) ~ maxct - (dp/wp)*(t-to),
			(t>(to+wp)) & (t<=(to+wp+wr)) ~ maxct - (dp-(dp/wr)*(t-(to+wp))),
			TRUE ~ maxct
			))

	return(out) 

	})
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

get_effective_sensitivity <- function(t_test, lod, se, trajectory_pars, event_duration=3/24, ndraws=1000){

	trajectory_df <- rtrajectory(ndraws, trajectory_pars, event_duration)

	out <- trajectory_df %>% 
		getct(t_test, trajectory_pars) %>% 
		runtest(lod=lod, se=se) %>%
		summarise(sensitivity=sum(screened)/n()) %>% 
		pull(sensitivity)

	return(out)

}

get_n_infectious <- function(t_test, lod, se, trajectory_pars, pop_pars, event_duration=3/24, ngames=1000, siglevel=0.9){

	with(as.list(c(trajectory_pars, pop_pars)),{

		n_infectious_baseline <- rbinom(ngames, n_attendees, prev)
		eff_se <- get_effective_sensitivity(t_test, lod, se, trajectory_pars, event_duration, ndraws=5000)
		n_infectious <- unlist(lapply(n_infectious_baseline, function(x) rbinom(1, x, 1-eff_se)))
		out <- tibble(
			mean=mean(n_infectious),
			lwr=quantile(n_infectious, (1-siglevel)/2), 
			upr=quantile(n_infectious, 1-(1-siglevel)/2)) %>%
			pivot_longer(everything(), names_to="statistic", values_to="value") %>%
			mutate(t=t_test)

		return(out)

		})

}