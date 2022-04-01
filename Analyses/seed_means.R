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
library(patchwork)

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
  filter(!is.na(SEED_PER_SPIKE), SEED_PER_SPIKE >0)

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
  warmup =2500, 
  iter = 5000, 
  thin = 1, 
  chains = 3
)
# Stan model -------------

## Run the model by calling stan()
## Save the outputs as rds files

# Stan model to fit a mean for each species. The model fits each with a shared sd. Because the data is wonky for some species (ELVI and ELRI have all 1s because they are measured as seed/infl and spikelets/infl are the same, the model has a hard time sampling)
sm_seed_mean <- stan(file = "Analyses/seed_mean.stan", data = seed_mean_data_list,
                     iter = mcmc_pars$iter,
                     warmup = mcmc_pars$warmup,
                     chains = mcmc_pars$chains, 
                     thin = mcmc_pars$thin)
saveRDS(sm_seed_mean, file = "~/Dropbox/EndodemogData/Model_Runs/endo_spp_seed_mean.rds")


#########################################################################################################
# Model Diagnostics ------------------------------
#########################################################################################################
print(sm_seed_mean)
traceplot(sm_seed_mean)


seedmean_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_seed_mean.rds") #this was our original shorter run that does not included beta0 and doesn't have spp as random effects

predSeedmean <- rstan::extract(seedmean_fit, pars = c("mu_seed"))$mu_seed
sdSeedmean <- rstan::extract(seedmean_fit, pars = c("sigma0"))$sigma0
n_post_draws <- 500
post_draws <- sample.int(dim(predSeedmean)[1], n_post_draws)
y_seedmean_sim <- matrix(NA,n_post_draws,length(seed_mean_data_list$seed))
for(i in 1:n_post_draws){
  y_seedmean_sim[i,] <- rnorm(n=length(seed_mean_data_list$seed), mean = predSeedmean[post_draws[i],], sd = sdSeedmean[post_draws[i]])
}
saveRDS(y_seedmean_sim, file = "yrep_seedmeanmodel.rds")
y_seedmean_sim <- readRDS(file = "yrep_seedmeanmodel.rds")

ppc_dens_overlay(seed_mean_data_list$seed, y_seedmean_sim)

seedmean_densplot <-ppc_dens_overlay(seed_mean_data_list$seed, y_seedmean_sim)+ theme_classic() + labs(title = "Seed Means", x = "Seeds per Infl.", y = "Density") 
ggsave(seedmean_densplot, filename = "seedmean_densplot.png", width = 4, height = 4)

mean_sm_plot <-   ppc_stat(seed_mean_data_list$seed, y_seedmean_sim, stat = "mean")
sd_sm_plot <- ppc_stat(seed_mean_data_list$seed, y_seedmean_sim, stat = "sd")
skew_sm_plot <- ppc_stat(seed_mean_data_list$seed, y_seedmean_sim, stat = "skewness")
kurt_sm_plot <- ppc_stat(seed_mean_data_list$seed, y_seedmean_sim, stat = "Lkurtosis")
seedmean_moments <- mean_sm_plot + sd_sm_plot + skew_sm_plot + kurt_sm_plot
seedmean_moments
ggsave(seedmean_moments, filename = "seedmean_moments.png", width = 4, height = 4)



