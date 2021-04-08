# =============================================================================
# Import
# =============================================================================

library(tidyverse) 
source('code/utilities/utils_analysis.R')
source('code/analysis/make_masterpars.R')
source("code/data_parsing/parse_Ct_data.R")
source("code/analysis/set_global_pars.R")

current_pars <- masterpars[[1]]

source('code/analysis/refine_data.R')

# =============================================================================
# Basic statistics
# =============================================================================

# Number of positive  tests per person: ---------------------------------------
n_positive_df <- ct_dat_clean %>% 
	mutate(ispositive=case_when(CT.Mean<40~1, TRUE~0)) %>%
	group_by(Person.ID) %>%
	summarise(n_positive=sum(ispositive), Symptomatic=max(Symptomatic)) %>%
	ungroup() %>% 
	mutate(Symptomatic=case_when(Symptomatic=="Yes"~1, TRUE~0))

n_positive_df %>% 
	ggplot(aes(x=n_positive)) + 
		geom_histogram(binwidth=1, fill="white", col="black") + 
		scale_x_continuous(breaks=1:9) + 
		theme_minimal() 

n_positive_df %>% 
	ungroup() %>% 
	summarise(mean=mean(n_positive), min=min(n_positive), Q25=quantile(n_positive, 0.25), median=median(n_positive), Q75=quantile(n_positive, 0.75), max=max(n_positive), n=n())

n_positive_df %>% 
	group_by(Symptomatic) %>% 
	summarise(mean=mean(n_positive), min=min(n_positive), Q25=quantile(n_positive, 0.25), median=median(n_positive), Q75=quantile(n_positive, 0.75), max=max(n_positive), n=n())

# Min Ct per person: ----------------------------------------------------------
min_ct_df <- ct_dat_clean %>% 
	group_by(Person.ID) %>% 
	summarise(minCt = min(CT.Mean), Symptomatic=max(Symptomatic)) %>%
	ungroup() %>% 
	mutate(Symptomatic=case_when(Symptomatic=="Yes"~1, TRUE~0))

min_ct_df %>% 
	ggplot(aes(x=minCt)) + 
		geom_histogram(binwidth=2, fill="white", col="black") + 
		# scale_x_continuous(breaks=1:9) + 
		theme_minimal() 

min_ct_df %>% 
	ungroup() %>% 
	summarise(mean=mean(minCt), min=min(minCt), Q25=quantile(minCt, 0.25), median=median(minCt), Q75=quantile(minCt, 0.75), max=max(minCt), n=n())

min_ct_df %>% 
	group_by(Symptomatic) %>% 
	summarise(mean=mean(minCt), min=min(minCt), Q25=quantile(minCt, 0.25), median=median(minCt), Q75=quantile(minCt, 0.75), max=max(minCt), n=n())

# =============================================================================
# Plot the raw data with some rough trajectories 
# =============================================================================

symptom_map <- ct_dat_refined %>% 
	mutate(Symptomatic=case_when(Symptomatic=="Yes"~1, TRUE~0)) %>% 
	group_by(Person.ID) %>%
	summarise(Symptomatic=max(Symptomatic)) %>%
	rename(id=Person.ID) %>%
	rename(symptomatic=Symptomatic)

symptom_color_df <- symptom_map %>% mutate(color=case_when(symptomatic==1~"red", TRUE~"blue"))
symptom_color_vec <- symptom_color_df$color
names(symptom_color_vec) <- symptom_color_df$id

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
	}) %>% 
	left_join(symptom_map, by="id") %>%
	filter(t >= -14)

bpfits <- indiv_data %>% 
	split(.$id) %>% 
	map(~ bpfit(., 0, "t","y")) %>% 
	map(~ filter(., data.type=="fit")) %>%
	map(~ select(., t, y)) %>% 
	bind_rows(.id="id") %>%
	mutate(id=as.integer(id)) %>% 
	left_join(symptom_map, by="id")

fig_bpfits <- bpfits %>% 
	ggplot(aes(x=t, y=y, col=factor(id))) + 
		geom_line(size=0.2, alpha=0.5) + 
		scale_y_reverse(limits=c(40,15)) + 
		scale_color_manual(values=c(symptom_color_vec,"symptomatic"="red","asymptomatic"="blue")) + 
		scale_x_continuous(breaks=seq(from=-14,to=35,by=7)) + 
		theme_minimal() + 
		theme(legend.position="none") + 
		labs(x="Days from min Ct", y="Ct")

fig_bpfits_withpoints <- fig_bpfits + 
	geom_point(data=indiv_data, aes(col=factor(id)), alpha=0.5, size=0.5) 

bpfits_overall <- indiv_data %>% 
	split(.$symptomatic) %>% 
	map(~ bpfit(., 0, "t","y")) %>% 
	map(~ filter(., data.type=="fit")) %>%
	map(~ select(., t, y)) %>% 
	bind_rows(.id="symptomatic") %>%
	mutate(symptomstatus=case_when(symptomatic==1~"symptomatic",TRUE~"asymptomatic"))

fig_bpfits_withpoints_withoverall <- fig_bpfits_withpoints + 
	geom_line(data=bpfits_overall, aes(col=symptomstatus)) 

fig_bpfits_withpoints_facet <- bpfits %>% 
	ggplot(aes(x=t, y=y, col=factor(id))) + 
		geom_line(size=0.2, alpha=0.5) + 
		scale_y_reverse(limits=c(40,15)) + 
		scale_color_manual(values=c(symptom_color_vec,"symptomatic"="red","asymptomatic"="blue")) + 
		scale_x_continuous(breaks=seq(from=-14,to=35,by=7)) + 
		theme_minimal() + 
		theme(legend.position="none") + 
		labs(x="Days from min Ct", y="Ct") + 
		geom_point(data=indiv_data, aes(col=factor(id)), alpha=0.5, size=0.5) +
		facet_wrap(~ factor(id))
