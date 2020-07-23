## Title: Grass endophyte population model with a bayesian framework
## Purpose: Builds matrix population model from vital rate estimates
## and runs stochastic LTRE
## Authors: Joshua and Tom
#############################################################

library(tidyverse)
library(scales)
library(bayesplot)
library(popbio)
library(countreg)
library(actuar)
library(rstan)

quote_bare <- function( ... ){
  substitute( alist(...) ) %>% 
    eval( ) %>% 
    sapply( deparse )
}

#############################################################################################
####### Read in Data and creating size bins------------------
#############################################################################################

source("Analyses/endodemog_data_processing.R")

max_size <- LTREB_full %>% 
  dplyr::select(species,size_t) %>% 
  filter(!is.na(size_t)) %>% 
  group_by(species) %>% 
  summarise(max_size = quantile(size_t,probs=0.975))


#############################################################################################
####### Read in matrix population functions ------------------
#############################################################################################

source("Analyses/MPM_functions.R")


#############################################################################################
####### Read in Stan vital rate model outputs ------------------
#############################################################################################

surv_fit_seedling <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_seedling_surv.rds")
surv_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_surv_woseedling.rds")
grow_fit_seedling <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_seedling_grow.rds")
grow_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_grow_PIG.rds")
flw_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_flw.rds")
fert_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_fert_noplot.rds")
spike_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_spike_year_plot.rds")
seedmean_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/seed_mean.rds")
stos_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_s_to_s.rds") 

surv_par <- rstan::extract(surv_fit, pars =quote_bare(beta0,betasize,betaendo,betaorigin,
                                                                   tau_year, tau_plot))
surv_sdlg_par <- rstan::extract(surv_fit_seedling, pars =quote_bare(beta0,betaendo,
                                                           tau_year, tau_plot))
grow_par <- rstan::extract(grow_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin,
                                                                    tau_year, tau_plot,
                                                       sigma))
grow_sdlg_par <- rstan::extract(grow_fit_seedling, pars = quote_bare(beta0,betaendo,
                                                       tau_year, tau_plot,
                                                       phi))
flow_par <- rstan::extract(flw_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin,
                                                       tau_year, tau_plot))
fert_par <- rstan::extract(fert_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin,
                                                      tau_year)) # plot effects didn't converge
spike_par <- rstan::extract(spike_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin,
                                                                         tau_year, tau_plot))
seed_par <- rstan::extract(seedmean_fit, pars = quote_bare(beta0,betaendo)) #no plot or year effect
recruit_par <- rstan::extract(stos_fit, pars = quote_bare(beta0,betaendo,
                                                          tau_year, tau_plot))

#############################################################################################
####### Run the MPM ------------------
#############################################################################################

# make the list of parameters and calculate mean lambdas
n_draws <- 100
post_draws <- sample.int(1300,size=n_draws) # this is smaller because of the flowering model that was run for a shorter number of iterations
lambda_mean <- array(dim = c(8,2,n_draws))

for(i in 1:length(post_draws)){
  for(e in 1:2){
    for(s in 1:7){
      lambda_mean[s,e,i] <- lambda2(bigmatrix(make_params(species=s,
                                                         endo_mean=(e-1),
                                                         endo_var=(e-1),
                                                         draw=post_draws[i],
                                                         max_size=max_size,
                                                         rfx=F,
                                                         surv_par=surv_par,
                                                         surv_sdlg_par = surv_sdlg_par,
                                                         grow_par=grow_par,
                                                         grow_sdlg_par = grow_sdlg_par,
                                                         flow_par=flow_par,
                                                         fert_par=fert_par,
                                                         spike_par=spike_par,
                                                         seed_par=seed_par,
                                                         recruit_par=recruit_par))$MPMmat)
    }
    lambda_mean[8,e,i] <- mean(lambda_mean[1:7,e,i])
  }
}

# Need to double check PIG shape parameter (sigma) which should be positive. If shape is negative, it gives NA's. Right now, we adjusted this from grow_mean*params$sigma to be exp(grow_mean)*params$sigma
# Double check this with Tom's POAR code
test_params <- make_params(species=3,
            endo_mean=(2-1),
            endo_var=(2-1),
            draw=post_draws[61],
            max_size=max_size,
            rfx=F,
            surv_par=surv_par,
            surv_sdlg_par = surv_sdlg_par,
            grow_par=grow_par,
            grow_sdlg_par = grow_sdlg_par,
            flow_par=flow_par,
            fert_par=fert_par,
            spike_par=spike_par,
            seed_par=seed_par,
            recruit_par=recruit_par)
sx(1, test_params)
sx_sdlg(test_params)
gxy(1,2, test_params)
gxy_sdlg(1,2, test_params)
pxy(1,2, test_params)

fx(1,test_params)
bigmatrix(test_params)
lambda(bigmatrix(test_params)$MPMmat)

# Mean endophyte difference and quantiles
lambda_mean_diff <- matrix(NA,8,7)
for(s in 1:8){
  lambda_mean_diff[s,1] = mean(lambda_mean[s,2,] - lambda_mean[s,1,])
  lambda_mean_diff[s,2:7] = quantile(lambda_mean[s,2,] - lambda_mean[s,1,],probs=c(0.05,0.125,0.25,0.75,0.875,0.95))
}

## now do variance in lambda 

##- need to figure out why this is giving NA's
lambda_hold <- array(dim = c(11,7,2,n_draws))
lambda_var <- array(dim = c(8,2,n_draws))
for(i in 1:length(post_draws)){
  for(e in 1:2){
    for(s in 1:7){
      for(y in 1:11){
        
        lambda_hold[y,s,e,i] <- lambda2(bigmatrix(make_params(species=s,
                                                             endo_mean=(e-1),
                                                             endo_var=(e-1),
                                                             draw=post_draws[i],
                                                             max_size=max_size,
                                                             rfx=T,
                                                             year=y,
                                                             surv_par=surv_par,
                                                             surv_sdlg_par = surv_sdlg_par,
                                                             grow_par=grow_par,
                                                             grow_sdlg_par = grow_sdlg_par,
                                                             flow_par=flow_par,
                                                             fert_par=fert_par,
                                                             spike_par=spike_par,
                                                             seed_par=seed_par,
                                                             recruit_par=recruit_par))$MPMmat)
      }
      lambda_var[s,e,i] <- sd(lambda_hold[,s,e,i])
    }
    lambda_var[8,e,i] <- mean(lambda_var[1:7,e,i])
  }
}

lambda_var_diff <- matrix(NA,8,7)
for(s in 1:8){
  lambda_var_diff[s,1] = mean(lambda_var[s,2,] - lambda_var[s,1,])
  lambda_var_diff[s,2:7] = quantile(lambda_var[s,2,] - lambda_var[s,1,],probs=c(0.05,0.125,0.25,0.75,0.875,0.95))
}



################################################################

mean_surv_coef <- lapply(rstan::extract(surv_fit, pars =quote_bare(beta0,betasize,betaendo,betaorigin,
                                                                    tau_year, tau_plot))
                         , colMeans)

mean_surv_seedling_coef <- lapply(rstan::extract(surv_fit_seedling, pars = quote_bare(beta0,betaendo,
                                                                             tau_year, tau_plot))
                                  , colMeans)

mean_grow_coef <- lapply(rstan::extract(grow_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin,
                                                                    tau_year, tau_plot))
                         , colMeans)

mean_grow_seedling_coef <- lapply(rstan::extract(grow_fit_seedling, pars = quote_bare(beta0,betaendo,
                                                                             tau_year, tau_plot))
                                  , colMeans)


mean_flw_coef <- lapply(rstan::extract(flw_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin,
                                                                    tau_year, tau_plot))
                         , colMeans)

mean_fert_coef <- lapply(rstan::extract(fert_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin,
                                                                    tau_year)) # No plot effect
                         , colMeans)

mean_spikelet_coef <- lapply(rstan::extract(spike_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin,
                                                                        tau_year, tau_plot))
                             , colMeans)

mean_s_to_s_coef <- lapply(rstan::extract(stos_fit, pars = quote_bare(beta0,betaendo,
                                                                      tau_year, tau_plot))
                           , colMeans)

mean_seedmean_coef <- lapply(rstan::extract(seedmean_fit, pars = quote_bare(beta0,betaendo)) #no plot or year effect
                              , colMeans)





