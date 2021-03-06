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
n_positive_df <- ct_dat_refined %>% 
	mutate(ispositive=case_when(CT.Mean<40~1, TRUE~0)) %>%
	group_by(Person.ID) %>%
	summarise(n_positive=sum(ispositive), Symptomatic=max(Symptomatic)) %>%
	ungroup() %>% 
	mutate(Symptomatic=case_when(Symptomatic=="Yes"~1, TRUE~0))

fig_npositive <- n_positive_df %>% 
	ggplot(aes(x=n_positive)) + 
		geom_histogram(binwidth=1, fill="white", col="black") + 
		scale_x_continuous(breaks=1:9) + 
		theme_minimal() +
		theme(text=element_text(size=16))
ggsave(fig_npositive, file="figures/datviz/npositive.pdf", width=8, height=5)

n_positive_df %>% 
	ungroup() %>% 
	summarise(mean=mean(n_positive), min=min(n_positive), Q25=quantile(n_positive, 0.25), median=median(n_positive), Q75=quantile(n_positive, 0.75), max=max(n_positive), n=n())

n_positive_df %>% 
	group_by(Symptomatic) %>% 
	summarise(mean=mean(n_positive), min=min(n_positive), Q25=quantile(n_positive, 0.25), median=median(n_positive), Q75=quantile(n_positive, 0.75), max=max(n_positive), n=n())

# Min Ct per person: ----------------------------------------------------------
min_ct_df <- ct_dat_refined %>% 
	group_by(Person.ID) %>% 
	summarise(minCt = min(CT.Mean), Symptomatic=max(Symptomatic)) %>%
	ungroup() %>% 
	mutate(Symptomatic=case_when(Symptomatic=="Yes"~1, TRUE~0))

fig_minct <- min_ct_df %>% 
	ggplot(aes(x=minCt)) + 
		geom_histogram(binwidth=2, fill="white", col="black") + 
		# scale_x_continuous(breaks=1:9) + 
		theme_minimal() +
		theme(text=element_text(size=16))
ggsave(fig_minct, file="figures/datviz/minct.pdf", width=8, height=5)

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
# write_csv(bpfits, "figure_data/Fig2/bpfits.csv")


fig_bpfits <- bpfits %>% 
	ggplot(aes(x=t, y=y, col=factor(id))) + 
		geom_line(size=0.2, alpha=0.5) + 
		# scale_y_reverse(limits=c(40,15)) + 
		scale_y_reverse(limits=c(40,15),breaks=c(40,35,30,25,20,15), labels=c("(-)","35","30","25","20","15"), sec.axis=sec_axis(~convert_Ct_logGEML(.), name=expression(log[10]~RNA~copies/ml))) +
		scale_color_manual(values=c(symptom_color_vec,"symptomatic"="red","asymptomatic"="blue")) + 
		scale_x_continuous(breaks=seq(from=-14,to=35,by=7)) + 
		theme_minimal() + 
		theme(legend.position="none", text=element_text(size=16)) + 
		labs(x="Days from lowest Ct", y="Ct")

fig_bpfits_withpoints <- fig_bpfits + 
	geom_point(data=indiv_data, aes(col=factor(id)), alpha=0.5, size=0.5) 
ggsave(fig_bpfits_withpoints,file="figures/datviz/bpfits_withpoints.pdf", width=8, height=5)


bpfits_overall <- indiv_data %>% 
	split(.$symptomatic) %>% 
	map(~ bpfit(., 0, "t","y")) %>% 
	map(~ filter(., data.type=="fit")) %>%
	map(~ select(., t, y)) %>% 
	bind_rows(.id="symptomatic") %>%
	mutate(symptomstatus=case_when(symptomatic==1~"symptomatic",TRUE~"asymptomatic"))
# write_csv(indiv_data, "figure_data/Fig2/indiv_data.csv")

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
		labs(x="Days from lowest Ct", y="Ct") + 
		geom_point(data=indiv_data, aes(col=factor(id)), alpha=0.5, size=0.5) +
		facet_wrap(~ factor(id))


# =============================================================================
# Gap between observations 
# =============================================================================

fig_interval <- ct_dat_clean %>% 
	group_by(Person.ID) %>% 
	mutate(Date.Index.Lag1=lag(Date.Index)) %>%
	mutate(DateDiff=Date.Index - Date.Index.Lag1) %>%
	filter(!is.na(DateDiff)) %>% 
	ggplot(aes(x=DateDiff)) + 
		geom_histogram(aes(y=..density..), binwidth=1, col="black", fill="gray", size=0.2) + 
		scale_x_continuous(limits=c(0,12.5), breaks=0:12, minor_breaks=0:12) + 
		labs(x="Interval between consecutive tests (days)", y="Proportion of tests") + 
		theme_minimal() + 
		theme(text=element_text(size=14))
ggsave(fig_interval,file="figures/datviz/interval.pdf", width=8, height=5)

# write_csv(ct_dat_clean %>% 
# 	group_by(Person.ID) %>% 
# 	mutate(Date.Index.Lag1=lag(Date.Index)) %>%
# 	mutate(DateDiff=Date.Index - Date.Index.Lag1) %>%
# 	filter(!is.na(DateDiff)) %>% 
# 	select(DateDiff),
# 	file="figure_data/FigS1/datediff.csv")

interval_summary <- ct_dat_clean %>% 
	group_by(Person.ID) %>% 
	mutate(Date.Index.Lag1=lag(Date.Index)) %>%
	mutate(DateDiff=Date.Index - Date.Index.Lag1) %>%
	filter(!is.na(DateDiff)) %>% 
	mutate(LEQ1=case_when(DateDiff<=1~1, TRUE~0)) %>%
	mutate(LEQ4=case_when(DateDiff<=4~1, TRUE~0)) %>%
	mutate(GreaterThan12=case_when(DateDiff>12~1, TRUE~0)) %>%
	ungroup() %>% 
	summarise(N=n(), 
		NLEQ1=sum(LEQ1), PropLEQ1=NLEQ1/N, 
		NLEQ4=sum(LEQ4), PropLEQ4=NLEQ4/N, 
		NOver12=sum(GreaterThan12), PropOver12=NOver12/N)
print(interval_summary) 

