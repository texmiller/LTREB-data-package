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
# within the section titled "Preparing datalists for Seed Means Kernel"
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
                             nSeed = length(na.omit(LTREB_data_for_seedmeans$SEED_PER_SPIKE)),
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

sink("endodemog_seed_mean.stan")
cat("
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
    ", fill = T)
sink()

stanmodel <- stanc("endodemog_seed_mean.stan")

## Run the model by calling stan()
## Save the outputs as rds files

sm_seed_mean <- stan(file = "endodemog_seed_mean.stan", data = seed_means_data_list,
                     iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(sm_seed_mean, file = "endodemog_seed_mean_all.rds")








smFESU <- stan(file = "endodemog_seed_mean.stan", data = FESU_seed_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smFESU, file = "endodemog_seed_mean_FESU.rds")

smLOAR <- stan(file = "endodemog_seed_mean.stan", data = LOAR_seed_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smLOAR, file = "endodemog_seed_mean_LOAR.rds")

smPOAL <- stan(file = "endodemog_seed_mean.stan", data = POAL_seed_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smPOAL, file = "endodemog_seed_mean_POAL.rds")

smPOSY <- stan(file = "endodemog_seed_mean.stan", data = POSY_seed_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smPOSY, file = "endodemog_seed_mean_POSY.rds")

# AGPE had data recorded as seed/spikelet already, but we are using this model to calculate seed/spikelet on average, so this is using the seed/spikelet info from AGPE
smAGPE<- stan(file = "endodemog_seed_mean.stan", data = AGPE_seed_data_list,
              iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smAGPE, file = "endodemog_seed_mean_AGPE.rds")


# ELRI and ELVI had data collected in a slightly different way. They recorded seeds/inflorescence, so this is our calculation for the mean seeds/infl
smELRI <- stan(file = "endodemog_seed_mean.stan", data = ELRI_seed_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smELRI, file = "endodemog_seed_mean_ELRI.rds")

smELVI <- stan(file = "endodemog_seed_mean.stan", data = ELVI_seed_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smELVI, file = "endodemog_seed_mean_ELVI.rds")

