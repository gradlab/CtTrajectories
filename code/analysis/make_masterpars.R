library(tidyverse) 
# Define the set of sensitivities to run: 

masterpars <- list()
masterpars[[1]] <- c(symptom_treatment="split", parametrization="waittime")
masterpars[[2]] <- c(symptom_treatment="combined", parametrization="waittime")