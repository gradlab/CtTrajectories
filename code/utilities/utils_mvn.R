library(tidyverse) 
library(geoR)

# Custom box-cox transform: 
bctransform <- function(y, lambda){
	if(lambda != 0){
		out <- (y^lambda - 1)/lambda
		} else {
		out <- log(y)
		}
	return(out)
}

# Custom inverse box-cox transform:
invbctransform <- function(y, lambda){
	if(lambda != 0){
		out <- (lambda*y + 1)^(1/lambda)
		} else {
		out <- exp(y)
		}
	return(out)
}

# Run optimal box-cox transform on a vector: 
boxcoxify <- function(data, col){
	vec <- data[[col]]
	lambda <- boxcoxfit(vec)$lambda
	vec_bc <- bctransform(vec, lambda)
	return(vec_bc)
}


# Generate a mean vector, a covariance matrix, and the associated box-cox lambda values for a set of data with some columns undergoing a box-cox transform. Define the untransformed columns in "rawcols" and the columns to be box-cox transformed in "bccols". Example: 
# bc_struct <- make_transformed_MVN(params_df, c("dp"), c("wp","wr"))
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


# Using the output from make_transformed_MVN, generate random draws from the target distribution. This takes as input a structure bc_struct with a mean vector, covariance structure, and box-cox lambda values. Example: 
# draws <- r_transformed_MVN(1000, bc_struct) 
r_transformed_MVN <- function(n, bc_struct){

	n_input <- n
	n <- 2*n # Draw twice as many as asked for to deal with complex issues

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
	draws <- na.omit(draws)
	draws <- draws[1:min(n_input, nrow(draws)),]
	if(nrow(draws) < n_input){
		print("Warning: Inverse Box-Cox transform led to many complex values. Returned tibble has fewer draws than requested.")
	}
	return(draws)

}


# // Some validation of the previous two functions: 

# draws <- r_transformed_MVN(10000, bc_struct)

# temp <- params_df[1:10000,] %>%
# 	select(dp,wp,wr) %>%
# 	bind_cols(draws %>% rename(dp_synth=dp, wp_synth=wp, wr_synth=wr))	

# temp %>% 
# 	select(wr, wr_synth) %>% 
# 	pivot_longer(everything()) %>% 
# 	ggplot(aes(x=value)) + 
# 		geom_histogram(aes(fill=name, y=..density..), position="identity", alpha=0.2) + 
# 		geom_density(aes(col=name), adjust=2)

# draws <- r_transformed_MVN(1000, bc_struct)
# cov(draws_filt)
# cov(select(params_df, dp,wp,wr))

# dev.new(); 
# params_df[1:10000,] %>% 
# 	ggplot(aes(x=dp, y=wr)) + 
# 		geom_point(size=0.1, alpha=0.2)

# dev.new();
# draws %>% 
# 	ggplot(aes(x=dp, y=wr)) + 
# 		geom_point(size=0.1, alpha=0.2)


