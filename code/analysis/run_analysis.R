# source('code/analysis/run_analysis.R')

# Import:
library(tidyverse) 
source('code/utilities/utils_analysis.R')
# In utils_private, define 'savdatdir' as the directory where you want to save the model output (it can get big) 
source('code/utilities/utils_private.R')

# Read data: 
source("code/data_parsing/parse_Ct_data.R")

# Generate essential values: 
source("code/analysis/set_global_pars.R")

# Generate list of parameter sets for sensitivity: 
source('code/analysis/make_masterpars.R')
final_fitlist <- list()

# Generate and save the posteriors: 
for(parset in 1:length(masterpars)){ #1:length(masterpars)

	print(paste0("STARTING PARSET ",parset," OF ",length(masterpars)))

	# set parameters: 
	current_pars <- masterpars[[parset]]

	# Refine data: 
	source('code/analysis/refine_data.R')

	# Fit posteriors: 
	source('code/analysis/fit_posteriors_preamble.R') 
	source('code/analysis/fit_posteriors.R') 	
	source('code/analysis/make_figures.R')
	source('code/analysis/save_figures.R')
	source('code/analysis/save_figures_png.R')

	if(parset==1){
		write_csv(params_df %>% 
			mutate(geml=convert_Ct_logGEML((global_pars[["lod"]] - dp))) %>%
			mutate(gemlA_trans=convert_Ct_logGEML((global_pars[["lod"]] - dpmeanA_trans))) %>%
			mutate(gemlS_trans=convert_Ct_logGEML((global_pars[["lod"]] - dpmeanS_trans))) %>% 
			select(id_clean, iteration, tp, dp, wp, wr, symptomatic, dpmeanA_trans, dpmeanS_trans, wpmeanA_trans, wpmeanS_trans, wrmeanA_trans, wrmeanS_trans), 
			file="output/params_df_split.csv")
	}

	if(parset==2){ ########
		write_csv(params_df %>%
			mutate(geml=convert_Ct_logGEML((global_pars[["lod"]] - dp))) %>%
			mutate(gemlA_trans=convert_Ct_logGEML((global_pars[["lod"]] - dpmeanA_trans))) %>%
			mutate(gemlS_trans=convert_Ct_logGEML((global_pars[["lod"]] - dpmeanS_trans))), 
			file="output/params_df_combined.csv")
	}
		
	# Save the figures
	final_fitlist[[parset]] <- ct_fit
	}

}

# save(final_fitlist, file=paste0(savedatdir,"/revisions/final_fitlist.RData"))

source('code/analysis/report_results.R')

# launch_shinystan_nonblocking(final_fitlist[[1]])

source('code/data_parsing/viz_ct_data.R')



# For post-hoc analysis: ======================================================
# First run up to start of loop, then:
load(paste0(savedatdir,"/revisions/final_fitlist.RData"))
parset <- 1
current_pars <- masterpars[[parset]]
source('code/analysis/refine_data.R')
source('code/analysis/fit_posteriors_preamble.R') 
source('code/analysis/report_results.R')

