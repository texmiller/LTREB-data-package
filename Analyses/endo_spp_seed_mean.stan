
    data { 
    int<lower=0> nSeed;                       // number of observations of seed/spikelet
    
    int<lower=0> nSpecies;                    // number of Species
    int<lower=0> species[nSeed];                     // Species
    
    //int<lower=1,upper=2> endo_index[nSeed];          // Endophyte index
    int<lower=0,upper=1> endo_01[nSeed];            // plant endophyte status
    //int<lower=0> nEndo;                             // number of endophyte statuses

    real<lower=0> seed[nSeed];               // number of seeds per spikelet
    }
    
    parameters {
    real<lower=0> sigma0;               // overall variance
    real betaspp[nSpecies];             // seed mean intercept for each species
    real<lower=0> sigmaspp[nSpecies];  // seed variance for each species
    
    real betaendo[nSpecies];           // species specific endophyte effect on mean
    }
    
    transformed parameters{
    real mu_seed[nSeed];

    for(n in 1:nSeed){
    mu_seed[n] = betaspp[species[n]] + betaendo[species[n]]*endo_01[n];
    }
    }
    
    model {
    // Priors
    sigma0 ~ normal(0,1);

    to_vector(betaspp) ~ normal(0,sigmaspp);
    to_vector(sigmaspp) ~ normal(0,1);
    
    to_vector(betaendo) ~ normal(0,10);
    // Likelihood
      seed ~ normal(mu_seed,sigma0);
    }
    
