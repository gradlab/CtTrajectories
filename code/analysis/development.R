# =============================================================================
# UNDER DEVELOPMENT: 
# =============================================================================

# Load fit and extract parameters: 
load("output/final_fitlist.RData")
source('code/analysis/fit_posteriors_preamble.R') 
ct_fit <- final_fitlist[[1]]
params <- rstan::extract(ct_fit)
indiv_params_df <- make_indiv_params_df(params, c("tp","dp","wp","wr"), n_indiv) %>% 
	rename(id_clean=id) %>% 
	left_join(id_map, by="id_clean") %>%
	left_join(symptom_map, by="id")
if(current_pars[["symptom_treatment"]]=="split"){
	shared_params_df <- make_shared_params_df(params, c("dpmeanA","wpmeanA","wrmeanA","dpmeanS","wpmeanS","wrmeanS","dpsd","wpsd","wrsd")) 
	} else {
	shared_params_df <- make_shared_params_df(params, c("dpmean","wpmean","wrmean","dpsd","wpsd","wrsd")) 
	}
params_df <- indiv_params_df %>% 
	left_join(shared_params_df, by="iteration") %>% 
	select(-iteration) 




# Procedure: 
# 	o Function should take "raw parameters" and "bc parameters" - things that don't need a box-cox transform and things that do 
# 	o Extract a vector of box-cox lambda for the bc parameters
# 	o Generate a new data frame with transformed variables
# 	o Calculate mean and covariance vector for the transformed data 
# 	o New function for generating a draw given mean, covariance, and an indicator of which parameters are box-coxified (with lambda for the ones that are and NA for the ones that aren't, maybe). 
#	o 


make_transformed_MVN <- function(data, rawcols, bccols){

	nraw <- length(rawcols)
	nbc <- length(bccols) 

	# Initialize vector for box-cox lambdas: 
	lambdavals <- rep(NA, nraw+nbc) 
	names(lambdavals) <- c(rawcols,bccols)

	data_subset <- select(data, all_of(c(rawcols, bccols))) 

	for(colind in 1:nbc){
		vec <- data[[bccols[colind]]]
		lambdavals[nraw+colind] <- boxcoxfit(vec)$lambda
		vec_bc <- bctransform(vec, lambdavals[nraw+colind])
		newcol <- tibble(vec_bc)
		names(newcol) <- paste0(bccols[colind],"_bc")
		data_subset <- bind_cols(data_subset, newcol)
		data_subset <- select(data_subset, -bccols[colind])
	}

	bc_meanvec <- colMeans(data_subset)
	bc_covmat <- cov(data_subset)

	return(list(meanvec=bc_meanvec, covmat=bc_covmat, lambdavals=lambdavals))

}
# bc_struct <- make_transformed_MVN(params_df, c("dp"), c("wp","wr"))


r_transformed_MVN <- function(n, bc_struct){

	varnames <- names(bc_struct$lambdavals) 
	draws <- rnorm(n*length(varnames)) 
	draws <- matrix(draws, ncol=length(varnames))

	meanmat <- t(as.matrix(reduce(rep(list(bc_struct$meanvec),n), bind_rows)))
	bc_L <- t(chol(bc_struct$covmat))
	
	draws <- data.frame(t(meanmat + bc_L %*% t(draws)))
	
	for(colind in 1:length(varnames)){
		lambdaval <- bc_struct$lambdavals[colind]
		if(!is.na(lambdaval)){
			draws[colind] <- invbctransform(draws[colind],lambdaval)
		}
	}

	names(draws) <- varnames
	draws <- as_tibble(draws)

	return(draws)

}
# draws <- r_transformed_MVN(1000, bc_struct) 


# draws <- r_transformed_MVN(1000, bc_struct)
# draws_filt <- filter(draws, !is.na(wp) & !is.na(wr))
# cov(draws_filt)
# cov(select(params_df, dp,wp,wr))










temp <- params_df
temp$wp_bc <- boxcoxify(temp, "wp")
temp$wr_bc <- boxcoxify(temp, "wr")
bc_meanvec <- temp %>% 
	select(dp, wp_bc, wr_bc) %>% 
	summarise(dp=mean(dp), wp_bc=mean(wp_bc), wr_bc=mean(wr_bc))
bc_covmat <- cov(select(temp, dp, wp_bc, wr_bc))
bc_L <- t(chol(bc_covmat))


temp <- boxcoxify(params_df, "wp")
	

params_df %>% 
	mutate(wp_bc = boxcoxify(., wp))


boxcoxfit(params_df$dp)


fig_dpapscatter <- params_df %>% 
	ggplot(aes(x=dp, y=ap)) + 
		geom_point(size=0.1, alpha=0.1) + 
		theme_minimal()  + 
		facet_wrap(~id)

params_df %>% 
	ggplot(aes(x=bccustom(wr,.168))) + 
		geom_histogram(bins=50) + 
		theme_minimal()


boxcoxfit(params_df$wr)$lambda




fig_dpwpscatter <- params_df %>% 
	ggplot(aes(x=dp, y=wp)) + 
		geom_point(size=0.1, alpha=0.1) + 
		theme_minimal()  + 
		facet_wrap(~id)


nperm <- 1000
trueval <- cor(params_df$dp, params_df$wp)
permvals <- unlist(pmap(list(
	rep(list(params_df$dp),nperm),
	rep(list(params_df$wp),nperm)),
	function(x,y) cor(sample(x), sample(y))
))
ggplot(data=tibble(x=permvals), aes(x=x)) + 
	geom_histogram() + 
	geom_vline(xintercept=trueval) + 
	theme_minimal()


# hist(log(params_df$wp), breaks=100)
# cov(select(params_df,dp,wp,wr))

cor(sample(params_df$wp), sample(params_df$dp))



