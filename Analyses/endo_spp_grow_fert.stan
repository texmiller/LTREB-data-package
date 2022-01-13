
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
    int<lower=0> y[N];                  // plant growth at time t+1 or fert at time t+1
    vector<lower=0>[N] logsize;             // plant size at time t for growth model or size t1 for fert model
    int<lower=0,upper=1> endo_01[N];            // plant endophyte status for surv model
    int<lower=0,upper=1> origin_01[N];          // plant origin status for surv model
}

parameters {
    // vr params
    vector[nSpp] beta0; //spp-specific      // predictor parameters as grand means and spp rfx
    vector[nSpp] betasize;                  //   spp specific size slope 
    vector[nSpp] betaendo;                  // spp specific endophyt effect 
    real betaorigin;               //  origin effect
    
    real tau_year[nSpp,nEndo,nYear];      // random year effect, unique to species and endo

    vector[nSpp] sigma0;                 // year variance
    vector[nSpp] sigmaendo;              // endo effect on variance
    
    vector[nPlot] tau_plot;        // random plot effect
    real<lower=0> sigma_plot;          // plot variance effect
    
    //neg bin overdispersion
    vector[nSpp] phi; 
}

transformed parameters {
    real lambda[N]; 
    real od[N];     // overdispersion parameter
    real sigma_year[nSpp,nEndo]; 

    // surv Linear Predictor
    for(n in 1:N){
    lambda[n] = beta0[spp[n]] + betasize[spp[n]]*logsize[n] + betaendo[spp[n]]*endo_01[n] +
    betaorigin*origin_01[n]
    + tau_year[spp[n],(endo_01[n]+1),year_t[n]] + tau_plot[plot[n]];
    
    od[n] = exp(phi[spp[n]]);
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
      sigma_plot ~ normal(0,.1);
    
 
      // species specific fixed effects

      beta0 ~ normal(0,5); 
      betasize ~ normal(0,5); 
      betaendo ~ normal(0,5); 
      betaorigin ~ normal(0,5); 
      sigma0 ~ normal(0,1); 
      sigmaendo ~ normal(0,1); 
      phi ~ normal(0,1);    
      
 //species endo year priors
    for(s in 1:nSpp){
          to_vector(tau_year[s,1,]) ~ normal(0,sigma_year[s,1]); // sample year effects for each species for each endo status
          to_vector(tau_year[s,2,]) ~ normal(0,sigma_year[s,2]);
    }

    for(n in 1:N){
      y[n] ~ neg_binomial_2_log(lambda[n],od[n]);
      target += -log1m(neg_binomial_2_log_lpmf(0 | lambda[n], od[n])); 
     }
}
    
