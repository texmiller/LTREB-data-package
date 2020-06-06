
    data { 
    int<lower=0> nSeed;                       // number of observations of seed/spikelet
    
    int<lower=0> nSpecies;                    // number of Species
    int<lower=0> species[nSeed];                     // Species
    
    int<lower=1,upper=2> endo_index[nSeed];          // Endophyte index
    int<lower=0,upper=1> endo_01[nSeed];            // plant endophyte status
    int<lower=0> nEndo;                             // number of endophyte statuses

    real<lower=0> seed[nSeed];               // number of seeds per spikelet
    }
    
    parameters {
    real b_seed;            // mean intercept
    real<lower=0> s_seed;   // 

    vector[nSpecies] b_spp[nEndo];
    real<lower=0> s_spp[nEndo]; // endo-species intercept
    real b_endo;
    }
    
    transformed parameters{
    real mu_seed[nSeed];
    for(n in 1:nSeed){
    mu_seed[n] = b_seed + b_endo*endo_01[n] + b_spp[endo_index[n],species[n]];
    }
    }
    
    model {
    // Priors
    b_seed ~ normal(0,10);
    b_endo ~ normal(0,10);
    to_vector(b_spp[1]) ~ normal(0,s_spp[1]);
    to_vector(b_spp[2]) ~ normal(0,s_spp[2]);
    // Likelihood
      seed ~ normal(mu_seed,s_seed);
    }
    
