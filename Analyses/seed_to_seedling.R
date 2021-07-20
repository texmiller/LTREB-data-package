## Title: Grass endophyte population model with a bayesian framework
## Purpose: Creates seed to seedling recruitment kernel written in STAN, 
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

source("Analyses/endodemog_data_processing.R")

# Loading in our vital rate model samples to estimate seed production
sm_seed_means <- read_rds("~/Dropbox/EndodemogData/Model_Runs/seed_mean.rds")
sm_spikelet <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_spike_year_plot.rds")


# I'm going to use LTREB_full to generate estimates of seed production and calculate our numbers of recruits
seed_mean_pars <- rstan::extract(sm_seed_means, pars = c("beta0", "betaendo", "sigma0"))
spikelet_pars <- rstan::extract(sm_spikelet, pars = c("beta0", "betasize", "betaendo", "betaorigin", 
                                                      "tau_plot", "sigma_plot", 
                                                      "tau_year", "sigma_year"))

LTREB_rec <- LTREB_full
seedmean_prob <- c(NA)
spikeperinf_prob <- c(NA)

for(i in 1:nrow(LTREB_rec)){
      seedmean_prob[i] <- sample(seed_mean_pars$beta0[,LTREB_rec$species_index[i]], size = 1) +
                            sample(seed_mean_pars$betaendo[,LTREB_rec$species_index[i]], size = 1)*LTREB_rec$endo_01[i]
      
      spikeperinf_prob[i] <- exp(sample(spikelet_pars$beta0[,LTREB_rec$species_index[i]], size = 1) +
                               sample(spikelet_pars$betasize[,LTREB_rec$species_index[i]], size = 1)*LTREB_rec$logsize_t[i] +
                                 sample(spikelet_pars$betaendo[,LTREB_rec$species_index[i]], size = 1)*LTREB_rec$endo_01[i] +
                                   sample(spikelet_pars$betaorigin[,LTREB_rec$species_index[i]], size = 1)*LTREB_rec$origin_01[i] +
                                     sample(spikelet_pars$betasize[,LTREB_rec$species_index[i]], size = 1)*LTREB_rec$logsize_t[i] +
                                       sample(spikelet_pars$tau_plot[,LTREB_rec$plot_index[i]], size = 1) +
                                         sample(spikelet_pars$tau_year[,LTREB_rec$year_t_index[i]], size = 1))
        
}
LTREB_rec$spikeperinf_pred <- rpois(n = nrow(LTREB_rec), lambda = spikeperinf_prob)
LTREB_rec$seedmean_prob <- seedmean_prob

LTREB_s_to_s_data <- LTREB_rec %>% 
  mutate(seed_est = as.integer(FLW_STAT_T*FLW_COUNT_T*spikeperinf_pred*seedmean_prob)) %>% 
  group_by(species, species_index, plot_index, year_t, year_t_index, year_t1, year_t1_index, endo_01, endo_index) %>% 
  summarize(tot_seed_t = as.integer(round(sum(seed_est, na.rm = TRUE))),
            tot_recruit_t1 = length(origin_01 == "1" & year_t == birth),
            samplesize = n()) %>% 
  filter(tot_seed_t>=tot_recruit_t1) #There are 373 rows where there are recruits but no seeds from the previous year, so we filter these out for now. Think about a seed bank
dim(LTREB_s_to_s_data)

# Create data lists to be used for the Stan model
s_to_s_data_list <- list(tot_recruit_t1 = LTREB_s_to_s_data$tot_recruit_t1,
                              tot_seed_t = LTREB_s_to_s_data$tot_seed_t,
                              endo_01 = as.integer(LTREB_s_to_s_data$endo_01),
                              endo_index = as.integer(LTREB_s_to_s_data$endo_index),
                              year_t = as.integer(LTREB_s_to_s_data$year_t_index),
                              plot = as.integer(LTREB_s_to_s_data$plot_index),
                              spp = LTREB_s_to_s_data$species_index,
                              N = nrow(LTREB_s_to_s_data),
                              nYear = as.integer(max(unique(LTREB_s_to_s_data$year_t_index))),
                              nPlot = max(unique(LTREB_s_to_s_data$plot_index)),
                              nSpp = length(unique(LTREB_s_to_s_data$species_index)),
                              nEndo = length(unique(LTREB_s_to_s_data$endo_01)))
str(s_to_s_data_list)

#########################################################################################################
# Stan model for seed to seedling recruitment rate ------------------------------
#########################################################################################################
## run this code to optimize computer system settings for MCMC
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## MCMC settings
mcmc_pars <- list(
  warmup = 5000, 
  iter = 10000, 
  thin = 1, 
  chains = 3
)

# Stan model -------------

## Run the model by calling stan()
## Save the outputs as rds files
sm_s_to_s <- stanc(file = "Analyses/endo_spp_s_to_s.stan")
sm_s_to_s <- stan(file = "Analyses/endo_spp_s_to_s.stan", data = s_to_s_data_list,
                     iter = mcmc_pars$iter,
                     warmup = mcmc_pars$warmup,
                     chains = mcmc_pars$chains, 
                     thin = mcmc_pars$thin,
                     control = list(max_treedepth = 15))
saveRDS(sm_s_to_s, file = "~/Dropbox/EndodemogData/Model_Runs/endo_spp_s_to_s.rds")




#########################################################################################################
# Model Diagnostics ------------------------------
#########################################################################################################

s_to_s_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_s_to_s.rds") 
s_to_s_par <- rstan::extract(s_to_s_fit, pars = c("p"))$p
predRecruit<- s_to_s_par

n_post_draws <- 100
post_draws <- sample.int(dim(predRecruit)[1], n_post_draws)
y_recruit_sim <- matrix(NA,n_post_draws,length(s_to_s_data_list$tot_recruit_t1))

for(i in 1:n_post_draws){
  y_recruit_sim[i,] <- rbinom(n = length(s_to_s_data_list$tot_recruit_t1), size = s_to_s_data_list$tot_seed_t, prob = invlogit(predRecruit[post_draws[i],]))
}
ppc_dens_overlay(s_to_s_data_list$tot_recruit_t1, y_recruit_sim)
ppc_dens_overlay(s_to_s_data_list$tot_recruit_t1, y_recruit_sim) +xlim(0,100) # doesn't seem to fit super great for some of the lower values
# I think I'm simulating the yrep correctly but I might need to check that
# overall mean and sd look okay.

mean_sm_plot <-   ppc_stat(s_to_s_data_list$tot_recruit_t1, y_recruit_sim, stat = "mean")
sd_sm_plot <- ppc_stat(s_to_s_data_list$tot_recruit_t1, y_recruit_sim, stat = "sd")
skew_sm_plot <- ppc_stat(s_to_s_data_list$tot_recruit_t1, y_recruit_sim, stat = "skewness")
kurt_sm_plot <- ppc_stat(s_to_s_data_list$tot_recruit_t1, y_recruit_sim, stat = "Lkurtosis")
grid.arrange(mean_sm_plot,sd_sm_plot,skew_sm_plot,kurt_sm_plot)
