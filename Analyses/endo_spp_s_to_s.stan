
  data { 
    int<lower=0> N;                       // number of observations
    int<lower=0> nYear;                       // number of years
    int<lower=0, upper=nYear> year_t[N];         // year of observation
    int<lower=0> nEndo;                       // number of endo treatments
    int<lower=1, upper=nEndo> endo_index[N];          // index for endophyte effect
    int<lower=0> nPlot;                         // number of plots
    int<lower=0, upper=nPlot> plot[N];                   // plot of observation
    int<lower=0> nSpp;                 // number of species
    int<lower=0, upper=nSpp> spp[N];   // species of observation
    int<lower=0> tot_recruit_t1[N];      // total recruits into the plot (response)
    int<lower=0> tot_seed_t[N];             // total seeds in plot (predictor)
    int<lower=0,upper=1> endo_01[N];            // plant endophyte status (predictor)
    }
    
    parameters {
  // surv params
  vector[nSpp] beta0; 
  vector[nSpp] betaendo; 

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
    p[n] = beta0[spp[n]] + betaendo[spp[n]]*endo_01[n] +
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
      sigma_plot ~ normal(0,1);
    
 
      // species specific fixed effects

      beta0 ~ normal(0,5); 
      betaendo ~ normal(0,5); 
      sigma0 ~ normal(0,1); 
      sigmaendo ~ normal(0,1); 

 //species endo year priors
    for(s in 1:nSpp){
          to_vector(tau_year[s,1,]) ~ normal(0,sigma_year[s,1]); // sample year effects for each species for each endo status
          to_vector(tau_year[s,2,]) ~ normal(0,sigma_year[s,2]);
    }
          
 // Likelihood
    tot_recruit_t1 ~ binomial_logit(tot_seed_t, p);

}
    

