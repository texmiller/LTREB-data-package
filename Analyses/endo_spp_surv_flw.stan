
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
    vector<lower=0>[N] logsize_t;             // plant size at time t for surv model
    int<lower=0,upper=1> endo_01[N];            // plant endophyte status for surv model
    int<lower=0,upper=1> origin_01[N];          // plant origin status for surv model
}

parameters {
    // vr params
    vector[nSpp] beta0;                  // predictor parameters as grand means and spp rfx
    vector[nSpp] betasize;                  //   spp specific size slope 
    vector[nSpp] betaendo;                  // spp specific endophyt effect 
    vector[nSpp] betaorigin;               // spp specific origin effect
    
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
    p[n] = beta0[spp[n]] + betasize[spp[n]]*logsize_t[n] + betaendo[spp[n]]*endo_01[n] + betaorigin[spp[n]]*origin_01[n]
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
      beta0 ~ normal(0,10);
      betasize ~ normal(0,10);
      betaendo ~ normal(0,10);
      betaorigin ~ normal(0,10);      
      sigma0 ~ normal(0,1);
      sigmaendo ~ normal(0,1);
      
    //species endo year priors
          to_vector(tau_year[1,1,]) ~ normal(0,sigma_year[1,1]); // sample year effects
          to_vector(tau_year[2,1,]) ~ normal(0,sigma_year[2,1]); 
          to_vector(tau_year[3,1,]) ~ normal(0,sigma_year[3,1]); 
          to_vector(tau_year[4,1,]) ~ normal(0,sigma_year[4,1]); 
          to_vector(tau_year[5,1,]) ~ normal(0,sigma_year[5,1]); 
          to_vector(tau_year[6,1,]) ~ normal(0,sigma_year[6,1]); 
          to_vector(tau_year[7,1,]) ~ normal(0,sigma_year[7,1]); 
          
          to_vector(tau_year[1,2,]) ~ normal(0,sigma_year[1,2]); 
          to_vector(tau_year[2,2,]) ~ normal(0,sigma_year[2,2]); 
          to_vector(tau_year[3,2,]) ~ normal(0,sigma_year[3,2]); 
          to_vector(tau_year[4,2,]) ~ normal(0,sigma_year[4,2]); 
          to_vector(tau_year[5,2,]) ~ normal(0,sigma_year[5,2]); 
          to_vector(tau_year[6,2,]) ~ normal(0,sigma_year[6,2]); 
          to_vector(tau_year[7,2,]) ~ normal(0,sigma_year[7,2]); 
    y ~ bernoulli_logit(p);
}
    
