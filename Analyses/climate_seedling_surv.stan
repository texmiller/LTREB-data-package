
data {
  // indices
  int<lower=0> nYear;                       // number of years
  int<lower=0> nPlot;                         // number of plots
  int<lower=0> nSpp;                        // number of host species
  int<lower=0> nEndo;                        // number of endophyte levels
  // vital rate data
  int<lower=0> N;                       // number of observations for surv model
  int<lower=0, upper=nYear> year_t[N];         // year of observation for surv model
  int<lower=0> plot[N];                   // plot of observation for surv model
  int<lower=0, upper=nSpp> spp[N];         // year of observation for surv model
  int<lower=0, upper=1> y[N];      // plant survival at time t+1 or flowering at time t
  // vector<lower=0>[N] logsize;             // plant size at time t for surv model (we are fitting this model for all seedlings size == 1)
  int<lower=0,upper=1> endo_01[N];            // plant endophyte status for surv model
  int<lower=0,upper=1> origin_01[N];          // plant origin status for surv model
  vector[N] spei;                           // 12 or 3 month spei for each climate year
  // vector[N] spei_nl;
}

parameters {
  // vr params
  vector[nSpp] beta0;                  // predictor parameters as grand means and spp rfx
  vector[nSpp] betaendo;                  // spp specific endophyt effect 

  // climate interactions
  matrix[nSpp,nEndo] betaspei_endo;              // endo by climate interaction
  // matrix[nSpp,nEndo] betaspei_nl_endo;              // endo by climate interaction non-linear

  real tau_year[nSpp,nEndo,nYear];      // random year effect, unique to species and endo
  
  vector[nSpp] sigma0;                 // year variance
  vector[nSpp] sigmaendo;              // endo effect on variance
  
  vector[nPlot] tau_plot;        // random plot effect
  real<lower=0> sigma_plot;          // plot variance effect
}

transformed parameters {
  real p[N];                           
  real sigma_year[nSpp,nEndo];
  
  // surv Linear Predictor
  for(n in 1:N){
    p[n] = beta0[spp[n]] + betaendo[spp[n]]*endo_01[n]
    + betaspei_endo[spp[n],endo_01[n]+1]*spei[n]
    // + betaspei_nl_endo[spp[n],endo_01[n]+1]*spei_nl[n]
    + tau_year[spp[n],(endo_01[n]+1),year_t[n]] + tau_plot[plot[n]];
  }
  
  // endo effect on variance
  for(s in 1:nSpp){
    for(d in 1:nEndo){
      sigma_year[s,d] = exp(sigma0[s] + sigmaendo[s]*(d-1));
    }
  }
}

model {
  // priors
  //this is plot variance
  tau_plot ~ normal(0,sigma_plot);
  sigma_plot ~ normal(0, 1);
  
  //fixed effect priors
  beta0 ~ normal(0,5);
  betaendo ~ normal(0,5);
  to_vector(betaspei_endo) ~ normal(0,5);
  // to_vector(betaspei_nl_endo) ~ normal(0,5);
  sigma0 ~ normal(0,1);
  sigmaendo ~ normal(0,1);
  
  //species endo year priors
    for(s in 1:nSpp){
          to_vector(tau_year[s,1,]) ~ normal(0,sigma_year[s,1]); // sample year effects for each species for each endo status
          to_vector(tau_year[s,2,]) ~ normal(0,sigma_year[s,2]);
    }
    
  y ~ bernoulli_logit(p);
}
