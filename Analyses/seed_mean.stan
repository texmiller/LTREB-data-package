
data { 
int<lower=0> N;                       // number of observations of seed/spikelet

int<lower=0> nSpecies;                    // number of Species
int<lower=0> species[N];                     // Species

//int<lower=1,upper=2> endo_index[N];          // Endophyte index
int<lower=0,upper=1> endo_01[N];            // plant endophyte status
//int<lower=0> nEndo;                             // number of endophyte statuses
vector<lower=0>[N] seed;               // number of seeds per spikelet
}

parameters {
vector[nSpecies] beta0;             // seed mean intercept for each species

vector[nSpecies] betaendo;           // species specific endophyte effect on mean
real<lower=0> sigma0; 
}

transformed parameters{
real mu_seed[N];

for(n in 1:N){
// mu_seed[n] = beta0 + betaspp[species[n]] + betaendo[species[n]]*endo_01[n];
mu_seed[n] = beta0[species[n]] + betaendo[species[n]]*endo_01[n];
}

}

model {
// Priors
// beta0 ~ normal(0,1);
beta0 ~ normal(0,5);
betaendo ~ normal(0,5);
sigma0 ~ normal(0,1);
// Likelihood

  seed ~ normal(mu_seed,sigma0);

}

