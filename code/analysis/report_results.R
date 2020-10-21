# source('code/analysis/report_results.R')

# Report for symptomatic/asymptomatic split: 

conflevel <- 0.95

parset <- 1;
source('code/analysis/loadfit.R')

reportsummary(
	"Symptomatic peak Ct:",
	global_pars[["lod"]] - shared_params_df$dpmeanS,
	conflevel)

reportsummary(
	"Symptomatic proliferation duration:",
	shared_params_df$wpmeanS,
	conflevel)

reportsummary(
	"Symptomatic clearance duration:",
	shared_params_df$wrmeanS,
	conflevel)

reportsummary(
	"Symptomatic acute infection duration:",
	shared_params_df$wpmeanS+shared_params_df$wrmeanS,
	conflevel)

print("")

reportsummary(
	"Asymptomatic peak Ct:",
	global_pars[["lod"]] - shared_params_df$dpmeanA,
	conflevel)

reportsummary(
	"Asymptomatic proliferation duration:",
	shared_params_df$wpmeanA,
	conflevel)

reportsummary(
	"Asymptomatic clearance duration:",
	shared_params_df$wrmeanA,
	conflevel)

reportsummary(
	"Asymptomatic acute infection duration:",
	shared_params_df$wpmeanA+shared_params_df$wrmeanA,
	conflevel)

print("")

parset <- 2;
source('code/analysis/loadfit.R')

reportsummary(
	"Overall peak Ct:",
	global_pars[["lod"]] - shared_params_df$dpmean,
	conflevel)

reportsummary(
	"Overall proliferation duration:",
	shared_params_df$wpmean,
	conflevel)

reportsummary(
	"Overall clearance duration:",
	shared_params_df$wrmean,
	conflevel)

reportsummary(
	"Overall acute infection duration:",
	shared_params_df$wpmean+shared_params_df$wrmean,
	conflevel)