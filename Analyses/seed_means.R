## Title: Grass endophyte population model with a bayesian framework
## Purpose: Creates seed production kernel written in STAN, 
## and does visualisation of posterior predictive checks
## Authors: Joshua and Tom
#############################################################

library(tidyverse)
library(rstan)
library(StanHeaders)
library(shinystan)
library(bayesplot)
library(devtools)

invlogit<-function(x){exp(x)/(1+exp(x))}
logit = function(x) { log(x/(1-x)) }


#############################################################################################
####### Data manipulation to prepare data as lists for Stan models------------------
#############################################################################################

# seed data lists are generated in the endodemog_data_processing.R file, 
source("Analyses/endodemog_data_processing.R")



##############################################################################
####### Preparing datalists for Seed Means Kernel ------------------------------
##############################################################################

LTREB_data_for_seedmeans <- LTREB_repro1 %>% 
  mutate(seed_per_spike = seed/spikelets) %>% 
  mutate(SEED_PER_SPIKE= case_when(species != "AGPE" ~ seed_per_spike,
                                   species == "AGPE" & tillerid_fixed == "multitillermean" ~ seed, # AGPE has some of its seeds data already recorded as seed/spike
                                   species == "AGPE" & tillerid_fixed != "multitillermean" ~ seed_per_spike)) %>% 
  mutate(species_index = as.integer(recode(species, !!!species_factor_key))) %>% 
  mutate(endo_index = as.integer(as.factor(endo_01+1)))  %>% 
  filter(!is.na(SEED_PER_SPIKE))

dim(LTREB_data_for_seedmeans)
# View(LTREB_data_for_seedmeans)

# Creating data list to be passed to the model
seed_means_data_list <- list(seed = LTREB_data_for_seedmeans$SEED_PER_SPIKE,
                             endo_01 = LTREB_data_for_seedmeans$endo_01,
                             endo_index = LTREB_data_for_seedmeans$endo_index,
                             year = LTREB_data_for_seedmeans$year,
                             plot = LTREB_data_for_seedmeans$plot_fixed,
                             species = LTREB_data_for_seedmeans$species_index,
                             N = length(na.omit(LTREB_data_for_seedmeans$SEED_PER_SPIKE)),
                             K = 2L,
                             nYear = length(unique(LTREB_data_for_seedmeans$year)),
                             nPlot = length(unique(LTREB_data_for_seedmeans$plot_fixed)),
                             nSpecies = length(unique(LTREB_data_for_seedmeans$species_index)),
                             nEndo =   length(unique(LTREB_data_for_seedmeans$endo_01)))
str(seed_means_data_list)



#########################################################################################################
# Stan model for mean of seed production per spikelet ------------------------------
#########################################################################################################
## run this code to optimize computer system settings for MCMC
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## MCMC settings
ni <-1000
nb <- 500
nc <- 3

# Stan model -------------
## here is the Stan model ##

sink("endo_spp_seed_mean.stan")
cat("
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
    to_vector(sigmaspp) ~ normal(0,1);
    
    to_vector(betaendo) ~ normal(0,10);
    // Likelihood
      seed ~ normal(mu_seed,sigma0);
    }
    ", fill = T)
sink()

stanmodel <- stanc("endo_spp_seed_mean.stan")

## Run the model by calling stan()
## Save the outputs as rds files

sm_seed_mean <- stan(file = "endo_spp_seed_mean.stan", data = seed_means_data_list,
                     iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(sm_seed_mean, file = "seed_means.rds")




