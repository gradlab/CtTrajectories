library(tidyverse) 

# Define a function to calculate the mean Ct given the parameters: 
mufun <- function(indiv_pars){
	with(indiv_pars,{
	# Initialize
	out <- rep(NA, length(t))
	# Ct = lod prior to infection: 
	out[t<=to] <- lod
	# Viral load rises between onset and peak: 
	riseinds <- (t>to & t<=(to+wp))
	out[riseinds] <- lod - dp/wp*(t[riseinds]-to)
	# Viral load falls between peak and recovery: 
	fallinds <- (t>(to+wp) & t<=(to+wp+wr))
	out[fallinds] <- (lod-dp) + dp/wr*(t[fallinds]-(to+wp))
	# Ct = lod after recovery: 
	out[t>(to+wp+wr)] <- lod
	# Return Ct values 
	return(out)
	})
} 

clean_test_result <- function(ct_dat){
	out <- ct_dat %>% 
		mutate(Test.Result=case_when(
			Test.Result%in%c("Positive for 2019-nCoV","Presumptive Positive for 2019-nCoV","DETECTED")~"Positive",
			Test.Result%in%c("Not Detected","NOT DETECTED")~"Not Detected",
			TRUE~"Inconclusive"
			))
	return(out)
}

clean_person_id <- function(ct_dat){
	id_list <- sort(unique(ct_dat$Person.ID))
	id_df <- data.frame(Person.ID=id_list, Person.ID.Clean=1:length(id_list), stringsAsFactors=FALSE)
	out <- left_join(ct_dat, id_df, by="Person.ID")
	return(out)
}

save_figlist <- function(obj,dir,name,driver,width=8,height=5){
	mapply(function(x,y) ggsave(x,
		file=paste0(dir,name,"_",y,".",driver), width=width, height=height),
	obj,1:length(obj))
}

clean_ct_dat <- function(ct_dat, global_pars){
	out <- with(as.list(global_pars),{
	ct_dat %>% 
		# ---------------------------------------------------------------------
		# Keep only people with >5 total observations: 
		group_by(Person.ID) %>%
		mutate(nobs=n()) %>% 
		filter(nobs>5) %>% 
		select(-nobs) %>%
		ungroup() %>% 
		# ---------------------------------------------------------------------
		# Keep only people whose Ct reaches below the LOD:
		group_by(Person.ID) %>%
		mutate(Min.CT = min(CT.Mean, na.rm=TRUE)) %>% 
		filter(Min.CT < lod) %>% 
		select(-Min.CT) %>%
		ungroup() %>%
		# ---------------------------------------------------------------------
		# Zero-out any Ct values above the limit of detection:
		mutate(CT.Mean=case_when(CT.Mean>lod~lod, TRUE~CT.Mean)) %>%
		# ---------------------------------------------------------------------
		# Replace dates with days since first sample: 
		make_date_index() %>% 
		# ---------------------------------------------------------------------
		# Replace NA CT.Means with the LOD:
		replace_na(list(CT.Mean=lod)) %>%
		# ---------------------------------------------------------------------
		# Keep just one record per accession number:
		group_by(Internal.Order.ID) %>% 
		filter(CT.Mean==min(CT.Mean)) %>%
		slice(1) %>%
		ungroup() %>% 
		# ---------------------------------------------------------------------
		# Reduce to a single observation per person per day, taking the min Ct: 
		# group_by(Person.ID, Date.Index) %>%
		# filter(CT.Mean==min(CT.Mean)) %>%
		# slice(1) %>%
		# ungroup() %>% 
		# ---------------------------------------------------------------------
		# Give people new ids: 
		clean_person_id()
		})
	return(out)
}

refine_ct_dat <- function(ct_dat, global_pars){
	out <- with(as.list(global_pars),{
	ct_dat %>% 
		# ---------------------------------------------------------------------
		# Extract useful columns:
		select(Person.ID, Symptomatic, Date.Index, CT.Mean, Adjusted) %>% 
		# ---------------------------------------------------------------------
		# Keep only people whose Ct reaches below 30 or have 2 consecutive <35:
		group_by(Person.ID) %>%
		mutate(Min.CT = min(CT.Mean, na.rm=TRUE)) %>% 
		mutate(OneUnder30 = case_when(Min.CT<30~1, TRUE~0)) %>% 
		select(-Min.CT) %>%
		group_by(Person.ID) %>% 
		mutate(CT.Mean.Lag1=lag(CT.Mean)) %>% 
		mutate(Consecutive=case_when(
			(CT.Mean<35) & (CT.Mean.Lag1<35) ~ 1,
			TRUE ~ 0
			)) %>% 
		group_by(Person.ID) %>%
		mutate(TwoUnder35 = max(Consecutive)) %>%
		filter((OneUnder30==1) | (TwoUnder35==1)) %>%
		select(-OneUnder30, -CT.Mean.Lag1, -Consecutive, -TwoUnder35) %>%
		ungroup() %>%
		# ---------------------------------------------------------------------
		# Give people new ids: 
		clean_person_id()
		})
	return(out)
}

make_date_index <- function(ct_dat){
	out <- ct_dat %>%
		group_by(Person.ID) %>% 
		arrange(Person.ID, Test.Date) %>% 
		mutate(First.Date=first(Test.Date)) %>% 
		mutate(Date.Index=as.integer(difftime(Test.Date,First.Date,units="days"))) %>%
		select(-Test.Date, -First.Date) %>% 
		ungroup()
	return(out)
}

center_date_index <- function(ct_dat){
	tminCt_df <- ct_dat %>% 
		group_by(Person.ID) %>% 
		arrange(Date.Index) %>% 
		filter(CT.Mean==min(CT.Mean, na.rm=TRUE)) %>%
		slice(1) %>%
		rename(tminCt=Date.Index) %>% 
		select(Person.ID, tminCt)
	out <- ct_dat %>%
		left_join(tminCt_df, by="Person.ID") %>%
		mutate(Date.Index=Date.Index-tminCt) %>% #print(n=50)
		select(-tminCt)
	return(out)
}

remove_one_off_positives <- function(ct_dat){
	out <- ct_dat %>% 
		group_by(Person.ID) %>% 
		mutate(CT.Mean.Lag1=lag(CT.Mean)) %>% 
		mutate(Consecutive=case_when(
			!is.na(CT.Mean) & !is.na(CT.Mean.Lag1) ~ 1,
			TRUE ~ 0
			)) %>% 
		mutate(Total.Consecutive=sum(Consecutive)) %>% 
		filter(Total.Consecutive>0) %>%
		ungroup() %>% 
		select(-CT.Mean.Lag1, -Consecutive, -Total.Consecutive)
	return(out) 
}

clean_ids <- function(indiv_data){
	old_id_list <- sort(unique(indiv_data$id)) 
	id_df <- data.frame(
		id=old_id_list, 
		new_id=1:length(old_id_list))
	out <- indiv_data %>% 
		left_join(id_df, by="id") %>% 
		select(-id) %>%
		rename(id=new_id) %>%
		select(names(indiv_data))
	return(out)
}

get_true_values <- function(indiv_pars, param){
	vals <- rep(NA, length(indiv_pars)) 
	for(id in 1:length(indiv_pars)){
		vals[id] <- indiv_pars[[id]][[param]]
	}
	out <- data.frame(id=1:length(indiv_pars), value=vals)
	names(out) <- c("id",param)
	return(out) 
}

# Function to jitter Ct values reasonably:
Ct_jitter_scalar <- function(y, global_pars){
	with(as.list(global_pars),{
		out <- min(round(if_else(y<lod,rnorm(1,y,sigma),lod)),lod)
		return(out)
		})}
Ct_jitter <- Vectorize(Ct_jitter_scalar, vectorize.args="y")

# For generating names for a matrix-turned-data frame: 
makenames <- function(parname, n_indiv){
	unlist(lapply(1:n_indiv, function(x) paste0(parname, "_", x)))
}

# For parsing the 'parameters' output from Stan: 
parseparam <- function(extracted_params, parname, n_indiv){
	as_tibble(setNames(as.data.frame(extracted_params[[parname]]), makenames(parname,n_indiv)))
}

make_params_df <- function(extracted_params, parnames, n_indiv){
	# Use "reduce" here
	out <- reduce(lapply(parnames, function(x) parseparam(extracted_params,x,n_indiv)), cbind) %>% 
		as_tibble %>% 
		mutate(iteration=1:n()) %>% 
		pivot_longer(-iteration) %>% 
		separate(name, c("param","id"), sep="_") %>%
		pivot_wider(c("id","iteration"), names_from="param", values_from="value") %>%
		select(-iteration) 
}

make_indiv_params_df <- function(extracted_params, parnames, n_indiv){
	# Use "reduce" here
	out <- reduce(lapply(parnames, function(x) parseparam(extracted_params,x,n_indiv)), cbind) %>% 
		as_tibble %>% 
		mutate(iteration=1:n()) %>% 
		pivot_longer(-iteration) %>% 
		separate(name, c("param","id"), sep="_") %>%
		pivot_wider(c("id","iteration"), names_from="param", values_from="value")
}


make_shared_params_df <- function(extracted_params, parnames){
	# Use "reduce" here
	out <- reduce(lapply(parnames, function(x) 
		as_tibble(setNames(as.data.frame(extracted_params[[x]]),x))
		), cbind) %>%
		as_tibble() %>%
		mutate(iteration=1:n())
}

# Plot the fitted Ct trajectories
plot_ct_fit <- function(params_df, global_pars, indiv_data, ctalpha=0.01, ntraces=500){
	with(as.list(global_pars),{
	params_df %>% 
		sample_n(ntraces) %>% 
		ggplot() + 
			# Plot traces:
			geom_segment(aes(x=-Inf, y=lod, xend=tp-wp, yend=lod), alpha=ctalpha) + 
			geom_segment(aes(x=tp-wp, y=lod, xend=tp, yend=lod-dp), alpha=ctalpha) + 
			geom_segment(aes(x=tp, y=lod-dp, xend=tp+wr, yend=lod), alpha=ctalpha) + 
			geom_segment(aes(x=tp+wr, y=lod, xend=Inf, yend=lod), alpha=ctalpha) + 
			# Plot data:
			geom_point(data=indiv_data, aes(x=t, y=y), size=0.5) + 
			theme_minimal() + 
			theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + 
			labs(x="Time (days)", y="Ct") + 
			scale_y_reverse() + 
			facet_wrap(~id)
			})
}

plot_ct_fit_symp <- function(params_df, global_pars, indiv_data, ctalpha=0.01, ntraces=500){
	with(as.list(global_pars),{
	params_df %>% 
		sample_n(ntraces) %>% 
		ggplot() + 
			# Plot traces:
			geom_segment(aes(x=-Inf, y=lod, xend=tp-wp, yend=lod), alpha=ctalpha) + 
			geom_segment(aes(x=tp-wp, y=lod, xend=tp, yend=lod-dp), alpha=ctalpha) + 
			geom_segment(aes(x=tp, y=lod-dp, xend=tp+wr, yend=lod), alpha=ctalpha) + 
			geom_segment(aes(x=tp+wr, y=lod, xend=Inf, yend=lod), alpha=ctalpha) + 
			# Plot data:
			geom_point(data=indiv_data, aes(x=t, y=y,col=factor(symptomatic)), size=0.5) + 
			scale_color_manual(values=c("1"="red","0"="blue")) + 
			theme_minimal() + 
			theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + 
			labs(x="Time (days)", y="Ct") + 
			scale_y_reverse() + 
			facet_wrap(~id)
			})
}

grid_off <- list(theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()))

y_ticks_off <- list(theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()))


plot_ct_fit_panes <- function(params_df, global_pars, indiv_data, ctalpha=0.01, ntraces=500, figsperpane=9){
	indiv_data <- indiv_data %>% mutate(plotgroup = floor((as.numeric(id_clean)-1)/figsperpane)+1)
	out <- with(as.list(global_pars),{
	params_df %>% 
		mutate(plotgroup = floor((as.numeric(id_clean)-1)/figsperpane)+1) %>%
		sample_n(ntraces) %>% 
		split(.$plotgroup) %>%
		imap(~ 
		ggplot(.x) + 
			# Plot traces:
			geom_segment(aes(x=-Inf, y=lod, xend=tp-wp, yend=lod), alpha=ctalpha) + 
			geom_segment(aes(x=tp-wp, y=lod, xend=tp, yend=lod-dp), alpha=ctalpha) + 
			geom_segment(aes(x=tp, y=lod-dp, xend=tp+wr, yend=lod), alpha=ctalpha) + 
			geom_segment(aes(x=tp+wr, y=lod, xend=Inf, yend=lod), alpha=ctalpha) + 
			# Plot data:
			geom_point(data=filter(indiv_data, plotgroup==as.numeric(.y)), aes(x=t, y=y), size=0.5) + 
			theme_minimal() + 
			theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + 
			labs(x="Time (days)", y="Ct") + 
			scale_y_reverse() + 
			facet_wrap(~id) )
			})
	return(out)
}

printnrows <- function(x, msg){
	print(paste(msg, nrow(x)))
	return(x)
}


printnpeople <- function(x, msg){
	npeople <- length(unique(x$Person.ID))
	print(paste(msg, npeople))
	return(x)
}

kernelplotfun <- function(adj){
	out <- list(
	geom_histogram(aes(y=..density..), bins=30, fill="white", col="gray"),
	geom_density(adjust=adj, col="gray", alpha=0.4),
	theme_minimal())
	return(out)
}

kernelplot <- list(
	geom_histogram(aes(y=..density..), bins=30, fill="white", col="gray"),
	geom_density(adjust=1.5),
	theme_minimal())

kernelplot_weighted <- list(
	geom_histogram(aes(y=..density.., weight=weight), bins=30, fill="white", col="gray"),
	geom_density(aes(weight=weight), adjust=1.5),
	theme_minimal())

plot_violins <- function(params_df, parname){
	params_df <- params_df %>% mutate(id=as.character(id))
	params_summary <- params_df %>% 
		group_by(id) %>% 
		summarise_(
			mean = interp(~mean(var), var=as.name(parname)),
			lwr = interp(~quantile(var,0.025), var=as.name(parname)),
			upr = interp(~quantile(var,0.975), var=as.name(parname))
			)
	fig_violin <- ggplot() + 
		geom_violin(data=params_df, aes_string(x="id", y=parname), col="gray", size=0.3) + 
		geom_segment(data=params_summary, aes(x=id, y=lwr, xend=id, yend=upr), col="black") + 
		geom_point(data=params_summary, aes(x=id, y=mean)) + 		
		theme_minimal() + 
		theme(panel.grid.major.x=element_blank(), panel.grid.major.y=element_line(size=0.2), panel.grid.minor=element_line(size=0.05))
	return(fig_violin)
}

plot_violins_overall <- function(params_df, parname, intwidth=0.9){
	params_df <- params_df %>% 
		mutate(id=factor(id)) %>%
		rbind(mutate(params_df, id="Overall")) 

	params_summary <- params_df %>% 
		group_by(id) %>% 
		summarise_(
			mean = interp(~mean(var), var=as.name(parname)),
			lwr = interp(~quantile(var,(1-intwidth)/2), var=as.name(parname)),
			upr = interp(~quantile(var,1-(1-intwidth)/2), var=as.name(parname))
			)
	fig_violin <- ggplot() + 
		geom_violin(data=params_df, aes_string(x="id", y=parname), col="gray", fill="lightgray", size=0.3) + 
		geom_segment(data=params_summary, aes(x=id, y=lwr, xend=id, yend=upr), col="black") + 
		geom_point(data=params_summary, aes(x=id, y=mean)) + 		
		theme_minimal() + 
		theme(panel.grid.major.x=element_blank(), panel.grid.major.y=element_line(size=0.2), panel.grid.minor=element_line(size=0.05))
	return(fig_violin)
}

plot_symptomatic_hists <- function(params_df, parname, adj=1.5){
	params_df_symp0 <- filter(params_df, symptomatic==0)
	params_df_symp1 <- filter(params_df, symptomatic==1)
	ggplot() + 
		geom_histogram(data=params_df_symp0, 
			aes_string(x=parname, y="..density.."), bins=30, fill="blue", alpha=0.1) + 
		geom_density(data=params_df_symp0, aes_string(x=parname), col="blue",adjust=adj) +
		geom_histogram(data=params_df_symp1, 
			aes_string(x=parname, y="..density.."), bins=30, fill="red", alpha=0.1) + 
		geom_density(data=params_df_symp1, aes_string(x=parname), col="red",adjust=adj) +
		scale_fill_identity(name = 'the fill', guide = 'legend', labels = c('m1')) + 
		theme_minimal() 
}

# Make data frame for dp prior distribution: 
make_unif_prior_df <- function(dmin, dmax, pmin, pmax, step){
	xmin <- qunif(p=pmin, min=dmin, max=dmax)
	xmax <- qunif(p=pmax, min=dmin, max=dmax)
	xvals <- seq(from=xmin, to=xmax, by=step)
	out <- as_tibble(data.frame(
		x=xvals, 
		density=dunif(xvals, min=dmin, max=dmax)))
	return(out)
}

make_gamma_prior_df <- function(shape, rate, pmin, pmax, step){
	xmin <- qgamma(pmin, shape=shape, rate=rate)
	xmax <- qgamma(pmax, shape=shape, rate=rate)
	xvals <- seq(from=xmin, to=xmax, by=step)
	out <- as_tibble(data.frame(
		x=xvals, 
		density=dgamma(xvals, shape=shape, rate=rate)))
	return(out)
}

make_normal_prior_df <- function(mean, sd, pmin, pmax, step){
	xmin <- qnorm(pmin, mean=mean, sd=sd)
	xmax <- qnorm(pmax, mean=mean, sd=sd)
	xvals <- seq(from=xmin, to=xmax, by=step)
	out <- as_tibble(data.frame(
		x=xvals, 
		density=dnorm(xvals, mean=mean, sd=sd)))
	return(out)
}

make_cauchy_prior_df <- function(location, scale, xmin, xmax, step){
	xvals <- seq(from=xmin, to=xmax, by=step)
	out <- as_tibble(data.frame(
		x=xvals, 
		density=dcauchy(xvals, location=location, scale=scale)))
	return(out)
}

# trim_negatives <- function(indiv_data, global_pars){
# 	out <- with(as.list(global_pars),{indiv_data %>% 
# 		split(.$id) %>% 
# 		map(~ arrange(., t)) %>% 
# 		map(~ mutate(., rowindex=1:n())) %>% 
# 		map(~ mutate(., ispositive=case_when(y<lod~1, TRUE~0))) %>% 
# 		map(~ list(data=., posindices=pull(filter(., ispositive==1),rowindex))) %>%
# 		map(~ list(data=.[[1]], slicemin=min(.[[2]]), slicemax=max(.[[2]]))) %>%
# 		map(~ slice(.[[1]], (.[[2]]-1):(.[[3]]+1))) %>% 
# 		map(~ select(., id, id_clean, t, y, adjusted)) %>%
# 		bind_rows()
# 	})
# 	return(out)
# }

# from https://www.martinmodrak.cz/2018/02/13/launch-shiny-app-without-blocking-the-session/

launch_shinystan_nonblocking <- function(fit) {
  library(future)
  plan(multisession)
  future(
    launch_shinystan(fit) #You can replace this with any other Shiny app
  )
}

trim_negatives <- function(indiv_data, global_pars){
	out <- with(as.list(global_pars),{indiv_data %>% 
		split(.$id) %>% 
		map(~ arrange(., t)) %>% 
		map(~ mutate(., rowindex=1:n())) %>% 
		map(~ mutate(., ispositive=case_when(y<lod~1, TRUE~0))) %>% 
		map(~ mutate(., ispositive_lag=lag(ispositive))) %>%
		map(~ mutate(., ispositive_lag2=lag(ispositive,2))) %>%
		map(~ mutate(., ispositive_lead=lead(ispositive))) %>%
		map(~ mutate(., ispositive_lead2=lead(ispositive,2))) %>%
		map(~ filter(., ispositive==1 | ispositive_lag==1 | ispositive_lag2==1 | ispositive_lead==1 | ispositive_lead2==1)) %>%
		map(~ select(., id, id_clean, t, y, adjusted)) %>%
		bind_rows()
	})
	return(out)
}


make_ct_probs <- function(ct_dat_clean, global_pars, ct_binwidth=2){
	out <- with(as.list(global_pars),{
		ct_dat_clean %>% 
		select(Person.ID, CT.Mean, Date.Index) %>% 
		filter(CT.Mean<lod) %>% 
		left_join(
			select(params_df, tp, wp, wr, id), by=c("Person.ID"="id")) %>%
		replace_na(list(tp=Inf, wp=Inf, wr=Inf)) %>%
		mutate(In.Onset=case_when(Date.Index>(tp-wp) & Date.Index<=tp ~ 1, TRUE~0)) %>%
		mutate(In.Resolution=case_when(Date.Index>tp & Date.Index<=(tp+wr) ~ 1, TRUE~0)) %>%
		group_by(Person.ID, Date.Index) %>%
		summarise(CT.Mean=first(CT.Mean), 
			Onset.Prob=sum(In.Onset)/n(), 
			Resolution.Prob=sum(In.Resolution)/n(), 
			Chatter.Prob=1-(Onset.Prob+Resolution.Prob)) %>% 
		mutate(CT.Bin=floor(CT.Mean/ct_binwidth)*ct_binwidth) %>%
		group_by(CT.Bin) %>% 
		summarise(
			Onset=sum(Onset.Prob)/n(), 
			Resolution=sum(Resolution.Prob)/n(),
			Chatter=sum(Chatter.Prob)/n()
			)			
		})
}

make_ct_probs_expanded <- function(ct_dat_clean, params_df, global_pars){
	with(as.list(global_pars),{
		ct_probs_expanded <- ct_dat_clean %>% 
			select(Person.ID, CT.Mean, Date.Index) %>% 
			filter(CT.Mean<lod) %>% 
			left_join(
				select(params_df, tp, wp, wr, id), by=c("Person.ID"="id")) %>%
			replace_na(list(tp=Inf, wp=Inf, wr=Inf)) %>%
			mutate(In.Onset=case_when(Date.Index>(tp-wp) & Date.Index<=tp ~ 1, TRUE~0)) %>%
			mutate(In.Resolution=case_when(Date.Index>tp & Date.Index<=(tp+wr) ~ 1, TRUE~0)) %>%
			group_by(Person.ID, Date.Index) %>%
			summarise(CT.Mean=first(CT.Mean), 
				Onset.Prob=sum(In.Onset)/n(), 
				Resolution.Prob=sum(In.Resolution)/n(), 
				Chatter.Prob=1-(Onset.Prob+Resolution.Prob)) %>%
			ungroup() 
		return(ct_probs_expanded)
		})
}

make_ct_probs_ma <- function(ct_dat_clean, params_df, global_pars, ct_binwidth=5, ct_increment=1){

	with(as.list(global_pars),{

		midpoints <- seq(from=lod-ct_binwidth/2,to=0+ct_binwidth/2,by=-ct_increment)

		ct_probs_expanded <- make_ct_probs_expanded(ct_dat_clean, params_df, global_pars)

		out <- data.frame() # Initialize output 

		for(mp in midpoints){
			new_out <- ct_probs_expanded %>% 
				filter(
					CT.Mean<=(mp+ct_binwidth/2) & 
					CT.Mean>(mp-ct_binwidth/2)) %>%
				summarise(
					Onset=sum(Onset.Prob, na.rm=TRUE)/n(),
					Resolution=sum(Resolution.Prob, na.rm=TRUE)/n(),
					Chatter=sum(Chatter.Prob, na.rm=TRUE)/n(),
					Overall_n=n()
					) %>%
				mutate(CT.Bin=mp) %>%
				select(CT.Bin, Onset, Resolution, Chatter, Overall_n) 
			out <- rbind(out, new_out)
			}

		return(filter(out, !is.na(Onset)))

		})
}

make_paired_ct_probs <- function(ct_dat_clean, params_df, global_pars, ct_binwidth=2){
	out <- with(as.list(global_pars),{
		ct_dat_clean %>% 
		select(Person.ID, CT.Mean, Date.Index) %>% 
		arrange(Person.ID, Date.Index) %>% 
		group_by(Person.ID) %>% 
		mutate(Date.Index.Diff=Date.Index-lag(Date.Index)) %>% 
		mutate(Ct.Change=case_when(lag(CT.Mean)<CT.Mean~"Increase", TRUE~"Decrease")) %>% 
		filter(Date.Index.Diff<=2) %>% 
		filter(CT.Mean<lod) %>% 
		left_join(
			select(params_df, tp, wp, wr, id), by=c("Person.ID"="id")) %>%
		replace_na(list(tp=Inf, wp=Inf, wr=Inf)) %>%
		mutate(In.Onset=case_when(Date.Index>(tp-wp) & Date.Index<=tp ~ 1, TRUE~0)) %>%
		mutate(In.Resolution=case_when(Date.Index>tp & Date.Index<=(tp+wr) ~ 1, TRUE~0)) %>%
		group_by(Person.ID, Date.Index) %>%
		summarise(CT.Mean=first(CT.Mean), 
			Date.Index.Diff=first(Date.Index.Diff), 
			Ct.Change=first(Ct.Change), 
			Onset.Prob=sum(In.Onset)/n(), 
			Resolution.Prob=sum(In.Resolution)/n(), 
			Chatter.Prob=1-(Onset.Prob+Resolution.Prob)) %>%
		mutate(CT.Bin=floor(CT.Mean/ct_binwidth)*ct_binwidth) %>%  
		group_by(CT.Bin, Ct.Change) %>%
		summarise(
			Onset=sum(Onset.Prob)/n(),
			Resolution=sum(Resolution.Prob)/n(),
			Chatter=sum(Chatter.Prob)/n(),
			)
		})
}

make_paired_ct_probs_expanded <- function(ct_dat_clean, params_df, global_pars){
	with(as.list(global_pars),{
		paired_ct_probs_expanded <- ct_dat_clean %>% 
			select(Person.ID, CT.Mean, Date.Index) %>% 
			arrange(Person.ID, Date.Index) %>% 
			group_by(Person.ID) %>% 
			mutate(Date.Index.Diff=lead(Date.Index)-Date.Index) %>% 
			mutate(Ct.Change=case_when(lead(CT.Mean)>=CT.Mean~"Increase", TRUE~"Decrease")) %>% 
			filter(Date.Index.Diff<=2) %>% 
			filter(CT.Mean<lod) %>% 
			left_join(
				select(params_df, tp, wp, wr, id), by=c("Person.ID"="id")) %>%
			replace_na(list(tp=Inf, wp=Inf, wr=Inf)) %>%
			mutate(In.Onset=case_when(Date.Index>(tp-wp) & Date.Index<=tp ~ 1, TRUE~0)) %>%
			mutate(In.Resolution=case_when(Date.Index>tp & Date.Index<=(tp+wr) ~ 1, TRUE~0)) %>%
			group_by(Person.ID, Date.Index) %>%
			summarise(CT.Mean=first(CT.Mean), 
				Date.Index.Diff=first(Date.Index.Diff), 
				Ct.Change=first(Ct.Change), 
				Onset.Prob=sum(In.Onset)/n(), 
				Resolution.Prob=sum(In.Resolution)/n(), 
				Chatter.Prob=1-(Onset.Prob+Resolution.Prob)) %>%
			ungroup() 
		return(paired_ct_probs_expanded)
		})
}

make_paired_ct_probs_ma <- function(ct_dat_clean, params_df, global_pars, ct_binwidth=5, ct_increment=1){

	with(as.list(global_pars),{

		midpoints <- seq(from=lod-ct_binwidth/2,to=0+ct_binwidth/2,by=-ct_increment)

		paired_ct_probs_expanded <- make_paired_ct_probs_expanded(ct_dat_clean, params_df, global_pars) 

		out <- data.frame() # Initialize output 

		for(mp in midpoints){
			new_out <- paired_ct_probs_expanded %>% 
				filter(
					CT.Mean<=(mp+ct_binwidth/2) & 
					CT.Mean>(mp-ct_binwidth/2)) %>%
				group_by(Ct.Change) %>% 
				summarise(
					Onset=sum(Onset.Prob, na.rm=TRUE)/n(),
					Resolution=sum(Resolution.Prob, na.rm=TRUE)/n(),
					Chatter=sum(Chatter.Prob, na.rm=TRUE)/n(),
					Overall_n=n()
					) %>%
				mutate(CT.Bin=mp) %>%
				select(CT.Bin, Ct.Change, Onset, Resolution, Chatter, Overall_n) 
			out <- rbind(out, new_out)
			}

		return(filter(out, !is.na(Onset)))

		})
}

plot_regime_probs <- function(ct_probs_df, lsspan=1){
	ct_probs_df %>% 
	pivot_longer(c("Onset","Resolution","Chatter"), names_to="Regime", values_to="Probability") %>% 
	ggplot(aes(x=CT.Bin, y=Probability, col=Regime)) + 
		geom_point() + 
		geom_line(stat="smooth", method="loess", span=lsspan) +
		scale_x_reverse() + 
		scale_y_continuous(breaks=seq(0,1,0.2), limits=c(-0.1,1.1)) + 
		scale_color_manual(values=c("Onset"="red","Resolution"="blue","Chatter"="black")) + 
		theme_minimal()
}

plot_infection_probs <- function(ct_probs_df, lsspan=1){
	ct_probs_df %>% 
	mutate(Infection=Onset+Resolution) %>% 
	select(-Onset, -Resolution) %>%
	pivot_longer(c("Infection","Chatter"), names_to="Regime", values_to="Probability") %>% 
	ggplot(aes(x=CT.Bin, y=Probability, col=Regime)) + 
		geom_point() + 
		geom_line(stat="smooth", method="loess", span=lsspan) +
		scale_x_reverse() + 
		scale_y_continuous(breaks=seq(0,1,0.2), limits=c(-0.1,1.1)) + 
		scale_color_manual(values=c("Infection"="red","Chatter"="black")) + 
		theme_minimal()
}

make_wald <- function(n, p, alpha, type){
	z <- qnorm(1-alpha/2)
	if(type=="upr"){
			out <- min(p + z*sqrt((p*(1-p))/n), 1)
			return(out)
		} else if(type=="lwr"){
			out <- max(p - z*sqrt((p*(1-p))/n), 0)
			return(out)
		} else {
			print("Invalid type")
			return(NaN)
		}
}

make_wald_v <- Vectorize(make_wald)

get_ks_vectors <- function(params_df, var){
	params_df %>% 
		split(.$symptomatic) %>% 
		map(~ .[[var]])	
}

reportsummary <- function(msg, data, interval, ndigits=1){
	print(paste0(msg," ",
	round(mean(data),digits=ndigits),
	" [",
	round(quantile(data,(1-interval)/2),digits=ndigits),", ",
	round(quantile(data,1-(1-interval)/2),digits=ndigits),
	"]"
	))
}

convert_Ct_logGEML <- function(Ct, m_conv=-3.609714286, b_conv=40.93733333){
	out <- (Ct-b_conv)/m_conv * log10(10) + log10(250)
	return(out) 
}


# the following four functions are for generating sample trajectory figures:
srise <- function(x, dp, wp){
	out <- dp*(1+x/wp)
	return(out)
}
sfall <- function(x, dp, wr){
	out <- dp*(1-x/wr)
	return(out)
}

make_sample_trajectory_symp <- function(shared_params_df, global_pars, siglevel=0.9, ge=FALSE){
	# For asymptomatic:
	with(as.list(global_pars),{

	shared_params_df <- shared_params_df %>% 
		select(-wpmeanA, -wpmeanS, -wrmeanA, -wrmeanS, -dpmeanA, -dpmeanS) %>%
		mutate(wpmeanA=wpmeanA_trans,
			wpmeanS=wpmeanS_trans,
			wrmeanA=wrmeanA_trans,
			wrmeanS=wrmeanS_trans,
			dpmeanA=dpmeanA_trans,
			dpmeanS=dpmeanS_trans)


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
	function(x) quantile(srise(x, shared_params_df$dpmeanA, shared_params_df$wpmeanA),1-(1-siglevel)/2)))
	yvals_lwr_proliferation_A <- unlist(lapply(xvals_proliferation_A, 
	function(x) quantile(srise(x, shared_params_df$dpmeanA, shared_params_df$wpmeanA),(1-siglevel)/2)))
	yvals_upr_clearance_A <- unlist(lapply(xvals_clearance_A, 
	function(x) quantile(sfall(x, shared_params_df$dpmeanA, shared_params_df$wrmeanA),1-(1-siglevel)/2)))
	yvals_lwr_clearance_A <- unlist(lapply(xvals_clearance_A, 
	function(x) quantile(sfall(x, shared_params_df$dpmeanA, shared_params_df$wrmeanA),(1-siglevel)/2)))

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
	function(x) quantile(srise(x, shared_params_df$dpmeanS, shared_params_df$wpmeanS),1-(1-siglevel)/2)))
	yvals_lwr_proliferation_S <- unlist(lapply(xvals_proliferation_S, 
	function(x) quantile(srise(x, shared_params_df$dpmeanS, shared_params_df$wpmeanS),(1-siglevel)/2)))
	yvals_upr_clearance_S <- unlist(lapply(xvals_clearance_S, 
	function(x) quantile(sfall(x, shared_params_df$dpmeanS, shared_params_df$wrmeanS),1-(1-siglevel)/2)))
	yvals_lwr_clearance_S <- unlist(lapply(xvals_clearance_S, 
	function(x) quantile(sfall(x, shared_params_df$dpmeanS, shared_params_df$wrmeanS),(1-siglevel)/2)))

	if(ge==FALSE){
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
			coord_cartesian(ylim=c(40,15), expand=FALSE) + 
			theme_minimal() + 
			labs(x="Days from peak", y="Ct") + 
			scale_y_reverse(sec.axis=sec_axis(~convert_Ct_logGEML(.), name=expression(log[10]~RNA~copies/ml))) + 
			theme(text=element_text(size=18))
	} else {
		out <- ggplot() + 
			geom_ribbon(
				data=tibble(
					xvals=xvals_proliferation_A, 
					yvals_lwr=(yvals_lwr_proliferation_A), 
					yvals_upr=(yvals_upr_proliferation_A)),
				aes(x=xvals, ymin=10^convert_Ct_logGEML(lod-yvals_lwr), ymax=10^convert_Ct_logGEML(lod-yvals_upr)), alpha=0.2, fill="blue") + 
			geom_segment(aes(x=-wp_mean_A,xend=0,y=10^convert_Ct_logGEML(lod),yend=10^convert_Ct_logGEML(lod-dp_mean_A)),col="blue") + 
			geom_ribbon(
				data=tibble(
					xvals=xvals_proliferation_S, 
					yvals_lwr=(yvals_lwr_proliferation_S), 
					yvals_upr=(yvals_upr_proliferation_S)),
				aes(x=xvals, ymin=10^convert_Ct_logGEML(lod-yvals_lwr), ymax=10^convert_Ct_logGEML(lod-yvals_upr)), alpha=0.2, fill="red") + 
			geom_segment(aes(x=-wp_mean_S,xend=0,y=10^convert_Ct_logGEML(lod),yend=10^convert_Ct_logGEML(lod-dp_mean_S)),col="red") + 
			geom_ribbon(
				data=tibble(
					xvals=xvals_clearance_A, 
					yvals_lwr=(yvals_lwr_clearance_A), 
					yvals_upr=(yvals_upr_clearance_A)),
				aes(x=xvals, ymin=10^convert_Ct_logGEML(lod-yvals_lwr), ymax=10^convert_Ct_logGEML(lod-yvals_upr)), alpha=0.2, fill="blue") + 
			geom_segment(aes(x=0,xend=wr_mean_A,y=10^convert_Ct_logGEML(lod-dp_mean_A),yend=10^convert_Ct_logGEML(lod)),col="blue") + 
			geom_ribbon(
				data=tibble(
					xvals=xvals_clearance_S, 
					yvals_lwr=(yvals_lwr_clearance_S), 
					yvals_upr=(yvals_upr_clearance_S)),
				aes(x=xvals, ymin=10^convert_Ct_logGEML(lod-yvals_lwr), ymax=10^convert_Ct_logGEML(lod-yvals_upr)), alpha=0.2, fill="red") + 
			geom_segment(aes(x=0,xend=wr_mean_S,y=10^convert_Ct_logGEML(lod-dp_mean_S),yend=10^convert_Ct_logGEML(lod)),col="red") + 
			coord_cartesian(ylim=c(10^convert_Ct_logGEML(40),10^convert_Ct_logGEML(15)), expand=FALSE) + 
			theme_minimal() + 
			labs(x="Days from peak", y="Genome equivalents per ml") + 
			scale_y_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + 
			theme(text=element_text(size=18))
	}

	return(out)

	})

}


make_sample_trajectory <- function(shared_params_df, global_pars, siglevel=0.9, ge=FALSE){
	# For asymptomatic:
	with(as.list(global_pars),{

	shared_params_df <- shared_params_df %>% 
		select(-wpmean, -wrmean, -dpmean) %>%
		mutate(wpmean=wpmean_trans,
			wrmean=wrmean_trans,
			dpmean=dpmean_trans)

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
	function(x) quantile(srise(x, shared_params_df$dpmean, shared_params_df$wpmean),1-(1-siglevel)/2)))
	yvals_lwr_proliferation <- unlist(lapply(xvals_proliferation, 
	function(x) quantile(srise(x, shared_params_df$dpmean, shared_params_df$wpmean),(1-siglevel)/2)))
	yvals_upr_clearance <- unlist(lapply(xvals_clearance, 
	function(x) quantile(sfall(x, shared_params_df$dpmean, shared_params_df$wrmean),1-(1-siglevel)/2)))
	yvals_lwr_clearance <- unlist(lapply(xvals_clearance, 
	function(x) quantile(sfall(x, shared_params_df$dpmean, shared_params_df$wrmean),(1-siglevel)/2)))

	if(ge==FALSE){
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
			coord_cartesian(ylim=c(40,15), expand=FALSE) + 
			theme_minimal() + 
			labs(x="Days from peak", y="Ct") + 
			scale_y_reverse(sec.axis=sec_axis(~convert_Ct_logGEML(.), name=expression(log[10]~RNA~copies/ml))) + 
			theme(text=element_text(size=18))
	} else {
		out <- ggplot() + 
			geom_ribbon(
				data=tibble(
					xvals=xvals_proliferation, 
					yvals_lwr=(yvals_lwr_proliferation), 
					yvals_upr=(yvals_upr_proliferation)),
				aes(x=xvals, ymin=10^convert_Ct_logGEML(lod-yvals_lwr), ymax=10^convert_Ct_logGEML(lod-yvals_upr)), alpha=0.2, fill="grey") + 
			geom_segment(aes(x=-wp_mean,xend=0,y=10^convert_Ct_logGEML(lod),yend=10^convert_Ct_logGEML(lod-dp_mean)),col="black") + 
			geom_ribbon(
				data=tibble(
					xvals=xvals_clearance, 
					yvals_lwr=(yvals_lwr_clearance), 
					yvals_upr=(yvals_upr_clearance)),
				aes(x=xvals, ymin=10^convert_Ct_logGEML(lod-yvals_lwr), ymax=10^convert_Ct_logGEML(lod-yvals_upr)), alpha=0.2, fill="grey") + 
			geom_segment(aes(x=0,xend=wr_mean,y=10^convert_Ct_logGEML(lod-dp_mean),yend=10^convert_Ct_logGEML(lod)),col="black") + 
			coord_cartesian(ylim=c(10^convert_Ct_logGEML(40),10^convert_Ct_logGEML(15)), expand=FALSE) + 
			theme_minimal() + 
			labs(x="Days from peak", y="Genome equivalents per ml") + 
			scale_y_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + 
			theme(text=element_text(size=18))
	}

	return(out)

	})

}

truncnormmean <- function(mu, sigma, T0, T1){
	alpha <- (T0 - mu)/sigma
	beta <- (T1 - mu)/sigma
	out <- mu + 
		(dnorm(alpha, 0, 1) - dnorm(beta, 0, 1))/
		(pnorm(beta, 0, 1) - pnorm(alpha, 0, 1))*sigma
	return(out)
}


# Given coefficients, this calculates the value of the breakpoint regression line with breakpoint at 't' for a given vector of values 'x': 
bpfunc <- function(params,t,x){
	b0 <- params[1]
	b1 <- params[2]
	b2 <- params[3]
	out <- b0 + b1*x + b2*(x-t)*as.numeric(x>t)
	return(out)
}

# This calculates the SSE for a breakpoint regression defined by the coefficients in 'params', the breakpoint at 't', and a data frame with independent variable x and dependent variable y: 
bpfunc_sse <- function(params,t,x,y){
	out <- sum((y - bpfunc(params,t,x))^2)
}

# This calculates the log likelihood of the best-fit breakpoint regression with breakpoint at time t and data given in 'df' with independent variable in 'xcol' and dependent variable in 'ycol': 
bpll <- function(df, t, xcol, ycol){
	x <- df[,xcol]
	y <- df[,ycol]
	z <- optim(c(1, 0, 0), bpfunc_sse, t=t, x=x, y=y)
	sseval <- bpfunc_sse(z$par,t,x,y)
	return(-sseval)
}

# This yields the best-fit breakpoint linear regression given data in 'df', a breakpoint at 't', dependentt variable in 'xcol', and dependent variable in 'ycol'. 
bpfit <- function(df, t, xcol, ycol){
	z <- optim(c(1, 0, 0), bpfunc_sse, t=t, x=df[,xcol], y=df[,ycol])
	xfitvals <- seq(min(df[,xcol]), max(df[,xcol]), length.out=100)
	yfitvals <- bpfunc(z$par, t, xfitvals)
	fitdf <- data.frame(xcoltemp=xfitvals, ycoltemp=yfitvals) %>% mutate(data.type="fit")
	names(fitdf)[1] <- xcol
	names(fitdf)[2] <- ycol
	df <- df %>% mutate(data.type="raw") %>% bind_rows(fitdf)
	return(df)
}


