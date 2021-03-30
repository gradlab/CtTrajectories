# source('code/analysis/report_results.R')

# Report for symptomatic/asymptomatic split: 

conflevel <- 0.95

parset <- 1;
source('code/analysis/loadfit.R')

reportsummary(
	"Symptomatic peak Ct:",
	global_pars[["lod"]] - shared_params_df$dpmeanS_trans,
	conflevel)

reportsummary(
	"Symptomatic proliferation duration:",
	shared_params_df$wpmeanS_trans,
	conflevel)

reportsummary(
	"Symptomatic clearance duration:",
	shared_params_df$wrmeanS_trans,
	conflevel)

reportsummary(
	"Symptomatic acute infection duration:",
	shared_params_df$wpmeanS_trans+shared_params_df$wrmeanS_trans,
	conflevel)

print("")

reportsummary(
	"Asymptomatic peak Ct:",
	global_pars[["lod"]] - shared_params_df$dpmeanA_trans,
	conflevel)

reportsummary(
	"Asymptomatic proliferation duration:",
	shared_params_df$wpmeanA_trans,
	conflevel)

reportsummary(
	"Asymptomatic clearance duration:",
	shared_params_df$wrmeanA_trans,
	conflevel)

reportsummary(
	"Asymptomatic acute infection duration:",
	shared_params_df$wpmeanA_trans+shared_params_df$wrmeanA_trans,
	conflevel)

print("")

parset <- 2;
source('code/analysis/loadfit.R')

reportsummary(
	"Overall peak Ct:",
	global_pars[["lod"]] - shared_params_df$dpmean_trans,
	conflevel)

reportsummary(
	"Overall proliferation duration:",
	shared_params_df$wpmean_trans,
	conflevel)

reportsummary(
	"Overall clearance duration:",
	shared_params_df$wrmean_trans,
	conflevel)

reportsummary(
	"Overall acute infection duration:",
	shared_params_df$wpmean_trans+shared_params_df$wrmean_trans,
	conflevel)