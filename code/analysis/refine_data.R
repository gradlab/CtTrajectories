ct_dat_refined <- with(as.list(c(global_pars, current_pars)),{
	ct_dat_clean %>% 
		printnrows("Rows in ct_dat_clean:") %>%
		printnpeople("People in ct_dat_clean:") %>%
		# ---------------------------------------------------------------------
		# Extract useful columns:
		select(Person.ID, Symptomatic, Date.Index, CT.Mean, Adjusted, Novel.Persistent.Infection) %>% 
		# ---------------------------------------------------------------------
		# Keep only people whose Ct reaches below 35:
		group_by(Person.ID) %>%
		mutate(Min.CT = min(CT.Mean, na.rm=TRUE)) %>% 
		filter(Min.CT <= 35) %>% 
		select(-Min.CT) %>%
		ungroup() %>%
		# ---------------------------------------------------------------------
		# Pull out people with persistent or recurrent infections if persistent_treatment is set to "exclude"
		filter(Novel.Persistent.Infection=="Novel") %>%
		# ---------------------------------------------------------------------
		# Filter out any individuals to exclude
		filter(!(Person.ID %in% current_pars[["toexclude"]])) %>%
		# ---------------------------------------------------------------------
		# Give people new ids: 
		clean_person_id()
		})