
    data { 
    int<lower=0> N;                       // number of observations of seed/spikelet
    
    int<lower=0> nSpecies;                    // number of Species
    int<lower=0> species[N];                     // Species
    
    //int<lower=1,upper=2> endo_index[N];          // Endophyte index
    int<lower=0,upper=1> endo_01[N];            // plant endophyte status
    //int<lower=0> nEndo;                             // number of endophyte statuses
    real<lower=0> seed[N];               // number of seeds per spikelet
    }
    
    parameters {
    real<lower=0> sigma0;               // overall variance
    real betaspp[nSpecies];             // seed mean intercept for each species
    real<lower=0> sigmaspp[nSpecies];  // seed variance for each species
    
    real betaendo[nSpecies];           // species specific endophyte effect on mean
    }
    
    transformed parameters{
    real mu_seed[N];
    for(n in 1:N){
    mu_seed[n] = betaspp[species[n]] + betaendo[species[n]]*endo_01[n];
    }
    }
    
    model {
    // Priors
    sigma0 ~ normal(0,1);
    to_vector(betaspp) ~ normal(0,sigmaspp);
    to_vector(sigmaspp) ~ normal(0,10);
    
    to_vector(betaendo) ~ normal(0,10);
    // Likelihood
      seed ~ normal(mu_seed,sigma0);
    }
    
