# source('code/analysis/run_analysis.R')

# Import:
library(tidyverse) 
source('code/utilities/utils_analysis.R')

# Generate list of parameter sets for sensitivity: 
source('code/analysis/make_masterpars.R')
final_fitlist <- list()

# Read data: 
source("code/data_parsing/parse_Ct_data.R")

# Generate essential values: 
source("code/analysis/set_global_pars.R")

# Generate and save the posteriors: 
for(parset in 3:4){ #length(masterpars)

	print(paste0("STARTING PARSET ",parset," OF ",length(masterpars)))

	# set parameters: 
	current_pars <- masterpars[[parset]]

	# Refine data: 
	source('code/analysis/refine_data.R')

	# Fit posteriors: 
	if(current_pars[["parametrization"]]=="waittime"){

		source('code/analysis/fit_posteriors_preamble.R') 
		source('code/analysis/fit_posteriors.R') 	
		source('code/analysis/make_figures.R')
		source('code/analysis/save_figures.R')
		source('code/analysis/save_figures_png.R')

	} else if(current_pars[["parametrization"]]=="slope"){

		source('code/analysis/fit_posteriors_preamble_slope.R') 
		source('code/analysis/fit_posteriors_slope.R') 	
		source('code/analysis/make_figures_slope.R')
		source('code/analysis/save_figures_slope.R')
		source('code/analysis/save_figures_png_slope.R')
	}	
		
	# # Save the figures

	final_fitlist[[parset]] <- ct_fit

}

# save(final_fitlist, file="output/final_fitlist.RData")

source('code/analysis/report_results.R')

# launch_shinystan_nonblocking(final_fitlist[[1]])