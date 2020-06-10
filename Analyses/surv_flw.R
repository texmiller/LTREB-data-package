## Title: Grass endophyte population model with a bayesian framework
## Purpose: Creates Survival, and, Flowering kernels written in STAN with mixed effects, 
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

#  data are prepared in the endodemog_data_processing.R file, 
source("Analyses/endodemog_data_processing.R")


#############################################################################################
####### Preparing data lists for vital rate kernels ------------------
#############################################################################################

## Clean up the main data frame for NA's, other small data entry errors
LTREB_data_forsurv <- LTREB_full %>% 
  filter(!is.na(surv_t1)) %>% 
  filter(!is.na(logsize_t)) %>% 
  filter(!is.na(endo_01)) # There are a few LOAR that don't have a plot level endo assigned
dim(LTREB_data_forsurv)

LTREB_data_forflw <- LTREB_full %>% 
  filter(!is.na(FLW_STAT_T)) %>% 
  filter(!is.na(logsize_t)) %>% 
  filter(!is.na(endo_01))
dim(LTREB_data_forflw)


# Create data lists to be used for the Stan model

surv_data_list <- list(y = LTREB_data_forsurv$surv_t1,
                            logsize_t = LTREB_data_forsurv$logsize_t,
                            origin_01 = LTREB_data_forsurv$origin_01,
                            endo_01 = as.integer(LTREB_data_forsurv$endo_01),
                            endo_index = as.integer(LTREB_data_forsurv$endo_index),
                            spp = as.integer(LTREB_data_forsurv$species_index),
                            year_t = as.integer(LTREB_data_forsurv$year_t_index),
                            plot = as.integer(LTREB_data_forsurv$plot_index),
                            N = nrow(LTREB_data_forsurv),
                            nSpp = length(unique(LTREB_data_forsurv$species_index)),
                            nYear = max(unique(LTREB_data_forsurv$year_t_index)),
                            nPlot = length(unique(LTREB_data_forsurv$plot_index)),
                            nEndo =   length(unique(LTREB_data_forsurv$endo_01)))
str(surv_data_list)

flw_data_list <- list(y = LTREB_data_forflw$FLW_STAT_T,
                           logsize_t = LTREB_data_forflw$logsize_t,
                           origin_01 = LTREB_data_forflw$origin_01,
                           endo_01 = as.integer(LTREB_data_forflw$endo_01),
                           endo_index = as.integer(LTREB_data_forflw$endo_index),
                           spp = as.integer(LTREB_data_forflw$species_index),
                           year_t = as.integer(LTREB_data_forflw$year_t_index),
                           plot = as.integer(LTREB_data_forflw$plot_index),
                           N = nrow(LTREB_data_forflw),
                           nSpp = length(unique(LTREB_data_forflw$species_index)),
                           nYear = max(unique(LTREB_data_forflw$year_t_index)),
                           nPlot = length(unique(LTREB_data_forflw$plot_index)),
                           nEndo =   length(unique(LTREB_data_forflw$endo_01)))
str(flw_data_list)




#########################################################################################################
# Stan model for Survival and Flowering ------------------------------
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
sink("Analyses/endo_spp_surv_flw.stan")
cat("
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
      sigma_plot ~ inv_gamma(0.001, 0.001);
      to_vector(tau_plot) ~ normal(0,sigma_plot);
      
    //fixed effect priors
      to_vector(beta0) ~ normal(0,100);
      to_vector(betasize) ~ normal(0,100);
      to_vector(betaendo) ~ normal(0,100);
      to_vector(betaorigin) ~ normal(0,100);      
      to_vector(sigma0) ~ normal(0,100);
      to_vector(sigmaendo) ~ normal(0,100);
      
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
    ", fill = T)
sink()

stanmodel <- stanc("Analyses/endo_spp_surv_flw.stan")

## Run the model by calling stan()
## Save the outputs as rds files

sm_surv <- stan(file = "Analyses/endo_spp_surv_flw.stan", data = surv_data_list,
                     iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(sm_surv, file = "endo_spp_surv.rds")

sm_flw <- stan(file = "endo_spp_surv_flw.stan", data = flw_data_list,
                iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(sm_flw, file = "endo_spp_flw.rds")


