
data {
  // indices
  int<lower=0> nYear;                       // number of years
  int<lower=0> nPlot;                         // number of plots
  int<lower=0> nSpp;                        // number of host species
  int<lower=0> nEndo;                        // number of endophyte levels
  // spike data
  int<lower=0> N;                       // number of observations for spike model
  int<lower=0, upper=nYear> year_t[N];         // year of observation for spike model
  int<lower=0> plot[N];                   // plot of observation for spike model
  int<lower=0, upper=nSpp> spp[N];         // year of observation for spike model
  int<lower=0> y[N];      // plant spikelet time t
  vector<lower=0>[N] logsize_t;             // plant size at time t for spike model
  int<lower=0,upper=1> endo_01[N];            // plant endophyte status for spike model
  int<lower=0,upper=1> origin_01[N];          // plant origin status for spike model
}

parameters {
  // surv params
  vector[nSpp] beta0; 
  vector[nSpp] betasize;
  vector[nSpp] betaendo; 
  vector[nSpp] betaorigin;  
  
  vector[nYear] tau_year;      // random year effect, 
  real<lower=0> sigma_year;          // year variance effect
  
  vector[nPlot] tau_plot;        // random plot effect
  real<lower=0> sigma_plot;          // plot variance effect
  
}

transformed parameters {
  real lambda[N]; 
  
  // surv Linear Predictor
  for(n in 1:N){
    lambda[n] = beta0[spp[n]] + 
      betasize[spp[n]]*logsize_t[n] + 
      betaendo[spp[n]]*endo_01[n] +
      betaorigin[spp[n]]*origin_01[n] + 
      tau_year[year_t[n]] + tau_plot[plot[n]];
  }
  
}

model {
  // priors
  
  sigma_plot ~ normal(0, 0.1);
  tau_plot ~ normal(0,sigma_plot);
  
  sigma_year ~ normal(0, 0.1);
  for(i in 1:nYear){
    tau_year[i] ~ normal(0,sigma_year);
  }
  

    beta0 ~ normal(0,10); 
    betasize ~ normal(0,10); 
    betaendo ~ normal(0,10); 
    betaorigin ~ normal(0,10); 

  y ~ poisson_log(lambda);
}
