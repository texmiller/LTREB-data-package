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
library(moments)
library(gridExtra)

invlogit<-function(x){exp(x)/(1+exp(x))}
logit = function(x) { log(x/(1-x)) }
Lkurtosis=function(x) log(kurtosis(x)); 


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
seed_mean_data_list <- list(seed = LTREB_data_for_seedmeans$SEED_PER_SPIKE,
                             endo_01 = LTREB_data_for_seedmeans$endo_01,
                             endo_index = LTREB_data_for_seedmeans$endo_index,
                             year = LTREB_data_for_seedmeans$year,
                             plot = LTREB_data_for_seedmeans$plot_fixed,
                             species = LTREB_data_for_seedmeans$species_index,
                             N = length(na.omit(LTREB_data_for_seedmeans$SEED_PER_SPIKE)),
                             nYear = length(unique(LTREB_data_for_seedmeans$year)),
                             nPlot = length(unique(LTREB_data_for_seedmeans$plot_fixed)),
                             nSpecies = length(unique(LTREB_data_for_seedmeans$species_index)),
                             nEndo =   length(unique(LTREB_data_for_seedmeans$endo_01)))
str(seed_mean_data_list)



#########################################################################################################
# Stan model for mean of seed production per spikelet ------------------------------
#########################################################################################################
## run this code to optimize computer system settings for MCMC
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## MCMC settings
mcmc_pars <- list(
  warmup = 5000, 
  iter = 10000, 
  thin = 3, 
  chains = 1
)
# Stan model -------------

## Run the model by calling stan()
## Save the outputs as rds files

sm_seed_mean <- stan(file = "Analyses/endo_spp_seed_mean.stan", data = seed_mean_data_list,
                     iter = mcmc_pars$iter,
                     warmup = mcmc_pars$warmup,
                     chains = mcmc_pars$chains, 
                     thin = mcmc_pars$thin)

sm_seed_mean_wo_beta0 <- stan(file = "Analyses/endo_spp_seed_mean.stan", data = seed_mean_data_list,
                     iter = mcmc_pars$iter,
                     warmup = mcmc_pars$warmup,
                     chains = mcmc_pars$chains, 
                     thin = mcmc_pars$thin)

sm_seed_mean_centered <- stan(file = "Analyses/endo_spp_seed_mean_centered.stan", data = seed_mean_data_list,
                              iter = mcmc_pars$iter,
                              warmup = mcmc_pars$warmup,
                              chains = mcmc_pars$chains, 
                              thin = mcmc_pars$thin)

sm_seed_mean_centered_w0_beta0 <- stan(file = "Analyses/endo_spp_seed_mean_centered.stan", data = seed_mean_data_list,
                              iter = mcmc_pars$iter,
                              warmup = mcmc_pars$warmup,
                              chains = mcmc_pars$chains, 
                              thin = mcmc_pars$thin)
saveRDS(sm_seed_mean, file = "~/Dropbox/EndodemogData/Model_Runs/seed_means.rds")


# Model test of non-centered Jun 18,

saveRDS(sm_seed_mean, file = "~/Dropbox/EndodemogData/Model_Runs/seed_mean.rds")
saveRDS(sm_seed_mean_wo_beta0, file = "~/Dropbox/EndodemogData/Model_Runs/seed_mean_noncentered_wo_beta0.rds")
saveRDS(sm_seed_mean_centered, file = "~/Dropbox/EndodemogData/Model_Runs/seed_mean_centered.rds")
saveRDS(sm_seed_mean_centered_w0_beta0, file = "~/Dropbox/EndodemogData/Model_Runs/seed_mean_centered_wo_beta0.rds")




#########################################################################################################
# Model Diagnostics ------------------------------
#########################################################################################################

seedmean_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/seed_means.rds") #this was our original shorter run that does not included beta0 and doesn't have spp as random effects
seedmean_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/seed_mean_noncentered.rds")
seedmean_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/seed_mean_noncentered_wo_beta0.rds")
seedmean_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/seed_mean_centered.rds")
seedmean_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/seed_mean_centered_wo_beta0.rds")

predSeedmean <- rstan::extract(seedmean_fit, pars = c("mu_seed"))$mu_seed
sdSeedmean <- rstan::extract(seedmean_fit, pars = c("sigma0"))$sigma0
n_post_draws <- 100
post_draws <- sample.int(dim(predSeedmean)[1], n_post_draws)
y_seedmean_sim <- matrix(NA,n_post_draws,length(seed_mean_data_list$seed))
for(i in 1:n_post_draws){
  y_seedmean_sim[i,] <- rnorm(n=length(seed_mean_data_list$seed), mean = invlogit(predSeedmean[post_draws[i],]), sd = sdSeedmean[post_draws[i]])
}
ppc_dens_overlay(seed_mean_data_list$seed, y_seedmean_sim)
 
mean_sm_plot <-   ppc_stat(seed_mean_data_list$seed, y_seedmean_sim, stat = "mean")
sd_sm_plot <- ppc_stat(seed_mean_data_list$seed, y_seedmean_sim, stat = "sd")
skew_sm_plot <- ppc_stat(seed_mean_data_list$seed, y_seedmean_sim, stat = "skewness")
kurt_sm_plot <- ppc_stat(seed_mean_data_list$seed, y_seedmean_sim, stat = "Lkurtosis")
grid.arrange(mean_sm_plot,sd_sm_plot,skew_sm_plot,kurt_sm_plot)



