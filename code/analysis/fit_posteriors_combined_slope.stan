functions {
  // mufun returns dCt, the delta Ct below the LOD:
  real mufun(real t, real tp, real wp, real wr, real dp){
    // Viral load rises between onset and peak: 
    if(t>(tp-wp) && t<=tp)
      return(dp/wp*(t-(tp-wp)));
    // Viral load falls between peak and recovery: 
    else if(t>tp && t<=(tp+wr))
      return(dp - dp/wr*(t-tp));
    // Ct = lod after recovery: 
    else 
      return(0);
  }
}

data {
  int<lower=0> N;             // Number of concatenated data points
  int<lower=0> n_id;          // Number of individuals 
  real<lower=0> lod;          // Limit of detection
  int<lower=0> id[N];         // Vector marking which datum belongs to which id
  int<lower=0> symp[n_id];    // Vector marking symptomatic ids
  real t[N];                  // Vector marking the time for each data point 
  real<lower=0, upper=lod> y[N];  // Concatenated data vector 
  real<lower=0> tpsd;         // Prior sd for the onset time (days)
  real<lower=0> dpmean_prior; // Prior mean peak Ct (delta from lod)
  real<lower=0> dpsd_prior;   // Prior sd peak Ct (delta from lod)
  real<lower=0> wpmax;
  real<lower=0> wrmax;
  real<lower=0> apmean_prior; // Prior mean proliferation duration
  real<lower=0> apsd_prior;   // Prior sd proliferation duration
  real<upper=0> armean_prior; // Prior mean clearance duration 
  real<lower=0> arsd_prior;   // Prior sd clearance duration 
  real<lower=0> sigma_max;    // Max allowable value for observation noise
  real<lower=0> sigma_prior_scale;  // Prior observation noise Cauchy scale
  real<lower=0, upper=1> lambda; // Mixture probability (~1-sensitivity)
  real<lower=0> fpmean;        // False positive mean Ct
  real<lower=0> epsilon[N];    // Obs. sd adjustment from Yale/FL regression
}

transformed data {
  real<lower=0, upper=lod> ydrop[N];  // Concatenated deviation from LOD 

  real loglambda;
  real log1mlambda;

  real dpcauchypriorscale;
  // real wpcauchypriorscale;
  // real wrcauchypriorscale;

  // real sigma;

  for(i in 1:N){
    ydrop[i] = lod-y[i];
  }

  loglambda = log(lambda);
  log1mlambda = log1m(lambda);

  // Define cauchy prior scales so thatt 90% of the half-distribution mass lies below the max cutoff for that parameter. 
  dpcauchypriorscale = lod/tan(pi()*(0.95-0.5));
  // wpcauchypriorscale = wpmax/tan(pi()*(0.95-0.5));
  // wrcauchypriorscale = wrmax/tan(pi()*(0.95-0.5));

}

parameters {

  real<lower=0, upper=lod> dpmean;   // Poplation peak Ct drop mean
  real<lower=dpmean/wpmax> apmean;              // Population onset-to-peak slope
  real<upper=-dpmean/wrmax> armean;              // Population peak-to-recovery slope 

  real<lower=0> dpsd;          // Poplation peak Ct drop sd
  real<lower=0> apsd;          // Population onset-to-peak time sd
  real<lower=0> arsd;          // Population peak-to-recovery time sd

  real tp[n_id];                        // Peak time
  real<lower=0, upper=lod> dp[n_id];    // Peak Ct drop


  real<lower=0> aptilde[n_id];               // Proliferation slope
  real<upper=0> artilde[n_id];               // Clearance slope
  
  real<lower=0, upper=sigma_max> sigma;    // Process noise during infection
}

transformed parameters {

  real<lower=0> wpmean;
  real<lower=0> wrmean;

  real<lower=0> wp[n_id];  // Onset-to-peak time
  real<lower=0> wr[n_id];  // Peak-to-recovery time 

  real<lower=0> ap[n_id];  // Proliferation slope
  real<upper=0> ar[n_id];  // Clearance slope

  real process_sd[N];      // Process noise
  real mu[N];              // Mean Ct trajectory

  wpmean = dpmean/apmean;
  wrmean = -dpmean/armean;

  for(i in 1:n_id){
    ap[i] = aptilde[i] + dp[i]/wpmax;
    ar[i] = artilde[i] - dp[i]/wrmax;
  }
  // No Jacobian adjustment needed since wpmax and wrmax are constants; see https://mc-stan.org/docs/2_25/stan-users-guide/vectors-with-varying-bounds.html

  for(i in 1:n_id){
    wp[i] = dp[i]/ap[i];                // Proliferation duration
    wr[i] = -dp[i]/ar[i];               // Clearance duration
  }

  for(i in 1:N){
    mu[i]=mufun(t[i], tp[id[i]], wp[id[i]], wr[id[i]], dp[id[i]]);
    process_sd[i]=sqrt((sigma*sigma) + (epsilon[i]*epsilon[i]));
  };
}


model {

  // Hierarchical priors:
  dpmean ~ normal(dpmean_prior,dpsd_prior) T[0,lod];
  apmean ~ normal(apmean_prior, apsd_prior) T[dpmeanS/wpmax,];
  armean ~ normal(armean_prior, arsd_prior) T[,-dpmeanS/wrmax];

  dpsd ~ cauchy(0,dpcauchypriorscale) T[0,];
  apsd ~ cauchy(0,5) T[0,];
  arsd ~ cauchy(0,5) T[0,];  

  sigma ~ cauchy(0,sigma_prior_scale) T[0,10];

  // Individual parameter specifications:
  tp ~ normal(0,tpsd);
  for(i in 1:n_id){
      dp[i] ~ normal(dpmean, dpsd) T[0, lod];
      ap[i] ~ normal(apmean, apsd) T[dp[i]/wpmax,];
      ar[i] ~ normal(armean, arsd) T[,-dp[i]/wrmax];
  } 

  // Main model specification: 
  for(i in 1:N){

    target += log_sum_exp(
      log1mlambda + normal_lpdf(ydrop[i] | mu[i], process_sd[i]),
      loglambda + exponential_lpdf(ydrop[i] | 1/fpmean));

    if (ydrop[i] < 0 || ydrop[i] > lod)
      target += negative_infinity();
    else
      target += -log_diff_exp(
        log1mlambda + normal_lcdf(lod | mu[i], process_sd[i]),
        log1mlambda + normal_lcdf(0 | mu[i], process_sd[i]));
    // see https://mc-stan.org/docs/2_18/reference-manual/sampling-statements-section.html

    }

}


