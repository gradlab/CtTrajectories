# source('code/analysis/report_results.R')

# Report for symptomatic/asymptomatic split: 

conflevel <- 0.95

parset <- 1; # ----------------------------------------------------------------
source('code/analysis/loadfit.R')

reportsummary(
	"Symptomatic peak Ct:",
	global_pars[["lod"]] - shared_params_df$dpmeanS_trans,
	conflevel)

reportsummary(
	"Symptomatic peak GE/ml:",
	convert_Ct_logGEML(global_pars[["lod"]] - shared_params_df$dpmeanS_trans),
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

reportsummary(
	"Symptomatic proliferation rate:",
	(shared_params_df$dpmeanS_trans)/(shared_params_df$wpmeanS_trans),
	conflevel)

reportsummary(
	"Symptomatic proliferation rate (GE/ml):",
	(convert_Ct_logGEML((global_pars[["lod"]] - shared_params_df$dpmeanS_trans))-(convert_Ct_logGEML(global_pars[["lod"]])))/(shared_params_df$wpmeanS_trans),
	conflevel)

reportsummary(
	"Symptomatic clearance rate:",
	(shared_params_df$dpmeanS_trans)/(shared_params_df$wrmeanS_trans),
	conflevel)

reportsummary(
	"Symptomatic clearance rate (GE/ml):",
	(convert_Ct_logGEML((global_pars[["lod"]] - shared_params_df$dpmeanS_trans))-(convert_Ct_logGEML(global_pars[["lod"]])))/(shared_params_df$wrmeanS_trans),
	conflevel)

print("")

reportsummary(
	"Asymptomatic peak Ct:",
	global_pars[["lod"]] - shared_params_df$dpmeanA_trans,
	conflevel)

reportsummary(
	"Asymptomatic peak GE/ml:",
	convert_Ct_logGEML(global_pars[["lod"]] - shared_params_df$dpmeanA_trans),
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

reportsummary(
	"Asymptomatic proliferation rate):",
	(shared_params_df$dpmeanA_trans)/(shared_params_df$wpmeanA_trans),
	conflevel)

reportsummary(
	"Asymptomatic proliferation rate (ge/ml):",
	(convert_Ct_logGEML((global_pars[["lod"]] - shared_params_df$dpmeanA_trans))-(convert_Ct_logGEML(global_pars[["lod"]])))/(shared_params_df$wpmeanA_trans),
	conflevel)

reportsummary(
	"Asymptomatic clearance rate:",
	(shared_params_df$dpmeanA_trans)/(shared_params_df$wrmeanA_trans),
	conflevel)

reportsummary(
	"Asymptomatic clearance rate (ge/ml):",
	(convert_Ct_logGEML((global_pars[["lod"]] - shared_params_df$dpmeanA_trans))-(convert_Ct_logGEML(global_pars[["lod"]])))/(shared_params_df$wrmeanA_trans),
	conflevel)

print("")

parset <- 2; # ----------------------------------------------------------------
source('code/analysis/loadfit.R')

reportsummary(
	"Overall peak Ct:",
	global_pars[["lod"]] - shared_params_df$dpmean_trans,
	conflevel)

reportsummary(
	"Overall peak Ct (ge/ml):",
	convert_Ct_logGEML(global_pars[["lod"]] - shared_params_df$dpmean_trans),
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

reportsummary(
	"Overall proliferation rate:",
	(shared_params_df$dpmean_trans)/(shared_params_df$wpmean_trans),
	conflevel)

reportsummary(
	"Overall proliferation rate (ge/ml):",
	(convert_Ct_logGEML((global_pars[["lod"]] - shared_params_df$dpmean_trans))-(convert_Ct_logGEML(global_pars[["lod"]])))/(shared_params_df$wpmean_trans),
	conflevel)

reportsummary(
	"Overall clearance rate:",
	(shared_params_df$dpmean_trans)/(shared_params_df$wrmean_trans),
	conflevel)

reportsummary(
	"Overall clearance rate (ge/ml):",
	(convert_Ct_logGEML((global_pars[["lod"]] - shared_params_df$dpmean_trans))-(convert_Ct_logGEML(global_pars[["lod"]])))/(shared_params_df$wrmean_trans),
	conflevel)

print("")

parset <- 3; # ----------------------------------------------------------------
source('code/analysis/loadfit.R')

reportsummary(
	"Symptomatic peak Ct, no 3047:",
	global_pars[["lod"]] - shared_params_df$dpmeanS_trans,
	conflevel)

reportsummary(
	"Symptomatic peak Ct, no 3047 (ge/ml):",
	convert_Ct_logGEML(global_pars[["lod"]] - shared_params_df$dpmeanS_trans),
	conflevel)

reportsummary(
	"Symptomatic proliferation duration, no 3047:",
	shared_params_df$wpmeanS_trans,
	conflevel)

reportsummary(
	"Symptomatic clearance duration, no 3047:",
	shared_params_df$wrmeanS_trans,
	conflevel)

reportsummary(
	"Symptomatic acute infection duration, no 3047:",
	shared_params_df$wpmeanS_trans+shared_params_df$wrmeanS_trans,
	conflevel)

reportsummary(
	"Symptomatic proliferation rate, no 3047:",
	(shared_params_df$dpmeanS_trans)/(shared_params_df$wpmeanS_trans),
	conflevel)

reportsummary(
	"Symptomatic proliferation rate, no 3047 (ge/ml):",
	(convert_Ct_logGEML((global_pars[["lod"]] - shared_params_df$dpmeanS_trans))-(convert_Ct_logGEML(global_pars[["lod"]])))/(shared_params_df$wpmeanS_trans),
	conflevel)

reportsummary(
	"Symptomatic clearance rate, no 3047:",
	(shared_params_df$dpmeanS_trans)/(shared_params_df$wrmeanS_trans),
	conflevel)

reportsummary(
	"Symptomatic clearance rate, no 3047 (ge/ml):",
	(convert_Ct_logGEML((global_pars[["lod"]] - shared_params_df$dpmeanS_trans))-(convert_Ct_logGEML(global_pars[["lod"]])))/(shared_params_df$wrmeanS_trans),
	conflevel)

print("")

reportsummary(
	"Asymptomatic peak Ct, no 3047:",
	global_pars[["lod"]] - shared_params_df$dpmeanA_trans,
	conflevel)

reportsummary(
	"Asymptomatic peak Ct, no 3047 (ge/ml):",
	convert_Ct_logGEML(global_pars[["lod"]] - shared_params_df$dpmeanA_trans),
	conflevel)

reportsummary(
	"Asymptomatic proliferation duration, no 3047:",
	shared_params_df$wpmeanA_trans,
	conflevel)

reportsummary(
	"Asymptomatic clearance duration, no 3047:",
	shared_params_df$wrmeanA_trans,
	conflevel)

reportsummary(
	"Asymptomatic acute infection duration, no 3047:",
	shared_params_df$wpmeanA_trans+shared_params_df$wrmeanA_trans,
	conflevel)

reportsummary(
	"Asymptomatic proliferation rate, no 3047:",
	(shared_params_df$dpmeanA_trans)/(shared_params_df$wpmeanA_trans),
	conflevel)

reportsummary(
	"Asymptomatic proliferation rate, no 3047 (ge/ml):",
	(convert_Ct_logGEML((global_pars[["lod"]] - shared_params_df$dpmeanA_trans))-(convert_Ct_logGEML(global_pars[["lod"]])))/(shared_params_df$wpmeanA_trans),
	conflevel)

reportsummary(
	"Asymptomatic clearance rate, no 3047:",
	(shared_params_df$dpmeanA_trans)/(shared_params_df$wrmeanA_trans),
	conflevel)

reportsummary(
	"Asymptomatic clearance rate, no 3047 (ge/ml):",
	(convert_Ct_logGEML((global_pars[["lod"]] - shared_params_df$dpmeanA_trans))-(convert_Ct_logGEML(global_pars[["lod"]])))/(shared_params_df$wrmeanA_trans),
	conflevel)

print("")

parset <- 4; # ----------------------------------------------------------------
source('code/analysis/loadfit.R')

reportsummary(
	"Symptomatic peak Ct, low sensitivity:",
	global_pars[["lod"]] - shared_params_df$dpmeanS_trans,
	conflevel)

reportsummary(
	"Symptomatic peak Ct, low sensitivity (ge/ml):",
	convert_Ct_logGEML((global_pars[["lod"]] - shared_params_df$dpmeanS_trans)),
	conflevel)

reportsummary(
	"Symptomatic proliferation duration, low sensitivity:",
	shared_params_df$wpmeanS_trans,
	conflevel)

reportsummary(
	"Symptomatic clearance duration, low sensitivity:",
	shared_params_df$wrmeanS_trans,
	conflevel)

reportsummary(
	"Symptomatic acute infection duration, low sensitivity:",
	shared_params_df$wpmeanS_trans+shared_params_df$wrmeanS_trans,
	conflevel)

reportsummary(
	"Symptomatic proliferation rate, low sensitivity:",
	(shared_params_df$dpmeanS_trans)/(shared_params_df$wpmeanS_trans),
	conflevel)

reportsummary(
	"Symptomatic proliferation rate, low sensitivity (ge/ml):",
	(convert_Ct_logGEML((global_pars[["lod"]] - shared_params_df$dpmeanS_trans))-(convert_Ct_logGEML(global_pars[["lod"]])))/(shared_params_df$wpmeanS_trans),
	conflevel)

reportsummary(
	"Symptomatic clearance rate, low sensitivity:",
	(shared_params_df$dpmeanS_trans)/(shared_params_df$wrmeanS_trans),
	conflevel)

reportsummary(
	"Symptomatic clearance rate, low sensitivity (ge/ml):",
	(convert_Ct_logGEML((global_pars[["lod"]] - shared_params_df$dpmeanS_trans))-(convert_Ct_logGEML(global_pars[["lod"]])))/(shared_params_df$wrmeanS_trans),
	conflevel)

print("")

reportsummary(
	"Asymptomatic peak Ct, low sensitivity:",
	global_pars[["lod"]] - shared_params_df$dpmeanA_trans,
	conflevel)

reportsummary(
	"Asymptomatic peak Ct, low sensitivity (ge/ml):",
	convert_Ct_logGEML((global_pars[["lod"]] - shared_params_df$dpmeanA_trans)),
	conflevel)

reportsummary(
	"Asymptomatic proliferation duration, low sensitivity:",
	shared_params_df$wpmeanA_trans,
	conflevel)

reportsummary(
	"Asymptomatic clearance duration, low sensitivity:",
	shared_params_df$wrmeanA_trans,
	conflevel)

reportsummary(
	"Asymptomatic acute infection duration, low sensitivity:",
	shared_params_df$wpmeanA_trans+shared_params_df$wrmeanA_trans,
	conflevel)

reportsummary(
	"Asymptomatic proliferation rate, low sensitivity:",
	(shared_params_df$dpmeanA_trans)/(shared_params_df$wpmeanA_trans),
	conflevel)

reportsummary(
	"Asymptomatic proliferation rate, low sensitivity (ge/ml):",
	(convert_Ct_logGEML((global_pars[["lod"]] - shared_params_df$dpmeanA_trans))-(convert_Ct_logGEML(global_pars[["lod"]])))/(shared_params_df$wpmeanA_trans),
	conflevel)

reportsummary(
	"Asymptomatic clearance rate, low sensitivity:",
	(shared_params_df$dpmeanA_trans)/(shared_params_df$wrmeanA_trans),
	conflevel)

reportsummary(
	"Asymptomatic clearance rate, low sensitivity (ge/ml):",
	(convert_Ct_logGEML((global_pars[["lod"]] - shared_params_df$dpmeanA_trans))-(convert_Ct_logGEML(global_pars[["lod"]])))/(shared_params_df$wrmeanA_trans),
	conflevel)

print("")

parset <- 5; # ----------------------------------------------------------------
source('code/analysis/loadfit.R')

reportsummary(
	"Symptomatic peak Ct, unconstrained wp/wr:",
	global_pars[["lod"]] - shared_params_df$dpmeanS_trans,
	conflevel)

reportsummary(
	"Symptomatic peak Ct, unconstrained wp/wr (ge/ml):",
	convert_Ct_logGEML(global_pars[["lod"]] - shared_params_df$dpmeanS_trans),
	conflevel)

reportsummary(
	"Symptomatic proliferation duration, unconstrained wp/wr:",
	shared_params_df$wpmeanS_trans,
	conflevel)

reportsummary(
	"Symptomatic clearance duration, unconstrained wp/wr:",
	shared_params_df$wrmeanS_trans,
	conflevel)

reportsummary(
	"Symptomatic acute infection duration, unconstrained wp/wr:",
	shared_params_df$wpmeanS_trans+shared_params_df$wrmeanS_trans,
	conflevel)

reportsummary(
	"Symptomatic proliferation rate, unconstrained wp/wr:",
	(shared_params_df$dpmeanS_trans)/(shared_params_df$wpmeanS_trans),
	conflevel)

reportsummary(
	"Symptomatic proliferation rate, unconstrained wp/wr (ge/ml):",
	(convert_Ct_logGEML((global_pars[["lod"]] - shared_params_df$dpmeanS_trans))-(convert_Ct_logGEML(global_pars[["lod"]])))/(shared_params_df$wpmeanS_trans),
	conflevel)

reportsummary(
	"Symptomatic clearance rate, unconstrained wp/wr:",
	(shared_params_df$dpmeanS_trans)/(shared_params_df$wrmeanS_trans),
	conflevel)

reportsummary(
	"Symptomatic clearance rate, unconstrained wp/wr (ge/ml):",
	(convert_Ct_logGEML((global_pars[["lod"]] - shared_params_df$dpmeanS_trans))-(convert_Ct_logGEML(global_pars[["lod"]])))/(shared_params_df$wrmeanS_trans),
	conflevel)

print("")

reportsummary(
	"Asymptomatic peak Ct, unconstrained wp/wr:",
	global_pars[["lod"]] - shared_params_df$dpmeanA_trans,
	conflevel)

reportsummary(
	"Asymptomatic peak Ct, unconstrained wp/wr (ge/ml):",
	convert_Ct_logGEML(global_pars[["lod"]] - shared_params_df$dpmeanA_trans),
	conflevel)

reportsummary(
	"Asymptomatic proliferation duration, unconstrained wp/wr:",
	shared_params_df$wpmeanA_trans,
	conflevel)

reportsummary(
	"Asymptomatic clearance duration, unconstrained wp/wr:",
	shared_params_df$wrmeanA_trans,
	conflevel)

reportsummary(
	"Asymptomatic acute infection duration, unconstrained wp/wr:",
	shared_params_df$wpmeanA_trans+shared_params_df$wrmeanA_trans,
	conflevel)

reportsummary(
	"Asymptomatic proliferation rate, unconstrained wp/wr:",
	(shared_params_df$dpmeanA_trans)/(shared_params_df$wpmeanA_trans),
	conflevel)

reportsummary(
	"Asymptomatic proliferation rate, unconstrained wp/wr (ge/ml):",
	(convert_Ct_logGEML((global_pars[["lod"]] - shared_params_df$dpmeanA_trans))-(convert_Ct_logGEML(global_pars[["lod"]])))/(shared_params_df$wpmeanA_trans),
	conflevel)

reportsummary(
	"Asymptomatic clearance rate, unconstrained wp/wr:",
	(shared_params_df$dpmeanA_trans)/(shared_params_df$wrmeanA_trans),
	conflevel)

reportsummary(
	"Asymptomatic clearance rate, unconstrained wp/wr (ge/ml):",
	(convert_Ct_logGEML((global_pars[["lod"]] - shared_params_df$dpmeanA_trans))-(convert_Ct_logGEML(global_pars[["lod"]])))/(shared_params_df$wrmeanA_trans),
	conflevel)

print("")

parset <- 6; # ----------------------------------------------------------------
source('code/analysis/loadfit.R')

reportsummary(
	"Symptomatic peak Ct, low wp/wr prior:",
	global_pars[["lod"]] - shared_params_df$dpmeanS_trans,
	conflevel)

reportsummary(
	"Symptomatic peak Ct, low wp/wr prior (ge/ml):",
	convert_Ct_logGEML(global_pars[["lod"]] - shared_params_df$dpmeanS_trans),
	conflevel)

reportsummary(
	"Symptomatic proliferation duration, low wp/wr prior:",
	shared_params_df$wpmeanS_trans,
	conflevel)

reportsummary(
	"Symptomatic clearance duration, low wp/wr prior:",
	shared_params_df$wrmeanS_trans,
	conflevel)

reportsummary(
	"Symptomatic acute infection duration, low wp/wr prior:",
	shared_params_df$wpmeanS_trans+shared_params_df$wrmeanS_trans,
	conflevel)

reportsummary(
	"Symptomatic proliferation rate, low wp/wr prior:",
	(shared_params_df$dpmeanS_trans)/(shared_params_df$wpmeanS_trans),
	conflevel)

reportsummary(
	"Symptomatic proliferation rate, low wp/wr prior (ge/ml):",
	(convert_Ct_logGEML((global_pars[["lod"]] - shared_params_df$dpmeanS_trans))-(convert_Ct_logGEML(global_pars[["lod"]])))/(shared_params_df$wpmeanS_trans),
	conflevel)

reportsummary(
	"Symptomatic clearance rate, low wp/wr prior:",
	(shared_params_df$dpmeanS_trans)/(shared_params_df$wrmeanS_trans),
	conflevel)

reportsummary(
	"Symptomatic clearance rate, low wp/wr prior (ge/ml):",
	(convert_Ct_logGEML((global_pars[["lod"]] - shared_params_df$dpmeanS_trans))-(convert_Ct_logGEML(global_pars[["lod"]])))/(shared_params_df$wrmeanS_trans),
	conflevel)

print("")

reportsummary(
	"Asymptomatic peak Ct, low wp/wr prior:",
	global_pars[["lod"]] - shared_params_df$dpmeanA_trans,
	conflevel)

reportsummary(
	"Asymptomatic peak Ct, low wp/wr prior (ge/ml):",
	convert_Ct_logGEML(global_pars[["lod"]] - shared_params_df$dpmeanA_trans),
	conflevel)

reportsummary(
	"Asymptomatic proliferation duration, low wp/wr prior:",
	shared_params_df$wpmeanA_trans,
	conflevel)

reportsummary(
	"Asymptomatic clearance duration, low wp/wr prior:",
	shared_params_df$wrmeanA_trans,
	conflevel)

reportsummary(
	"Asymptomatic acute infection duration, low wp/wr prior:",
	shared_params_df$wpmeanA_trans+shared_params_df$wrmeanA_trans,
	conflevel)

reportsummary(
	"Asymptomatic proliferation rate, low wp/wr prior:",
	(shared_params_df$dpmeanA_trans)/(shared_params_df$wpmeanA_trans),
	conflevel)

reportsummary(
	"Asymptomatic proliferation rate, low wp/wr prior (ge/ml):",
	(convert_Ct_logGEML((global_pars[["lod"]] - shared_params_df$dpmeanA_trans))-(convert_Ct_logGEML(global_pars[["lod"]])))/(shared_params_df$wpmeanA_trans),
	conflevel)

reportsummary(
	"Asymptomatic clearance rate, low wp/wr prior:",
	(shared_params_df$dpmeanA_trans)/(shared_params_df$wrmeanA_trans),
	conflevel)

reportsummary(
	"Asymptomatic clearance rate, low wp/wr prior (ge/ml):",
	(convert_Ct_logGEML((global_pars[["lod"]] - shared_params_df$dpmeanA_trans))-(convert_Ct_logGEML(global_pars[["lod"]])))/(shared_params_df$wrmeanA_trans),
	conflevel)

print("")

parset <- 7; # ----------------------------------------------------------------
source('code/analysis/loadfit.R')

reportsummary(
	"Symptomatic peak Ct, high wp/wr prior:",
	global_pars[["lod"]] - shared_params_df$dpmeanS_trans,
	conflevel)

reportsummary(
	"Symptomatic peak Ct, high wp/wr prior (ge/ml):",
	convert_Ct_logGEML(global_pars[["lod"]] - shared_params_df$dpmeanS_trans),
	conflevel)

reportsummary(
	"Symptomatic proliferation duration, high wp/wr prior:",
	shared_params_df$wpmeanS_trans,
	conflevel)

reportsummary(
	"Symptomatic clearance duration, high wp/wr prior:",
	shared_params_df$wrmeanS_trans,
	conflevel)

reportsummary(
	"Symptomatic acute infection duration, high wp/wr prior:",
	shared_params_df$wpmeanS_trans+shared_params_df$wrmeanS_trans,
	conflevel)

reportsummary(
	"Symptomatic proliferation rate, high wp/wr prior:",
	(shared_params_df$dpmeanS_trans)/(shared_params_df$wpmeanS_trans),
	conflevel)

reportsummary(
	"Symptomatic proliferation rate, high wp/wr prior (ge/ml):",
	(convert_Ct_logGEML((global_pars[["lod"]] - shared_params_df$dpmeanS_trans))-(convert_Ct_logGEML(global_pars[["lod"]])))/(shared_params_df$wpmeanS_trans),
	conflevel)

reportsummary(
	"Symptomatic clearance rate, high wp/wr prior:",
	(shared_params_df$dpmeanS_trans)/(shared_params_df$wrmeanS_trans),
	conflevel)

reportsummary(
	"Symptomatic clearance rate, high wp/wr prior (ge/ml):",
	(convert_Ct_logGEML((global_pars[["lod"]] - shared_params_df$dpmeanS_trans))-(convert_Ct_logGEML(global_pars[["lod"]])))/(shared_params_df$wrmeanS_trans),
	conflevel)

print("")

reportsummary(
	"Asymptomatic peak Ct, high wp/wr prior:",
	global_pars[["lod"]] - shared_params_df$dpmeanA_trans,
	conflevel)

reportsummary(
	"Asymptomatic peak Ct, high wp/wr prior (ge/ml):",
	convert_Ct_logGEML(global_pars[["lod"]] - shared_params_df$dpmeanA_trans),
	conflevel)

reportsummary(
	"Asymptomatic proliferation duration, high wp/wr prior:",
	shared_params_df$wpmeanA_trans,
	conflevel)

reportsummary(
	"Asymptomatic clearance duration, high wp/wr prior:",
	shared_params_df$wrmeanA_trans,
	conflevel)

reportsummary(
	"Asymptomatic acute infection duration, high wp/wr prior:",
	shared_params_df$wpmeanA_trans+shared_params_df$wrmeanA_trans,
	conflevel)

reportsummary(
	"Asymptomatic proliferation rate, high wp/wr prior:",
	(shared_params_df$dpmeanA_trans)/(shared_params_df$wpmeanA_trans),
	conflevel)

reportsummary(
	"Asymptomatic proliferation rate, high wp/wr prior (ge/ml):",
	(convert_Ct_logGEML((global_pars[["lod"]] - shared_params_df$dpmeanA_trans))-(convert_Ct_logGEML(global_pars[["lod"]])))/(shared_params_df$wpmeanA_trans),
	conflevel)

reportsummary(
	"Asymptomatic clearance rate, high wp/wr prior:",
	(shared_params_df$dpmeanA_trans)/(shared_params_df$wrmeanA_trans),
	conflevel)

reportsummary(
	"Asymptomatic clearance rate, high wp/wr prior (ge/ml):",
	(convert_Ct_logGEML((global_pars[["lod"]] - shared_params_df$dpmeanA_trans))-(convert_Ct_logGEML(global_pars[["lod"]])))/(shared_params_df$wrmeanA_trans),
	conflevel)

print("")