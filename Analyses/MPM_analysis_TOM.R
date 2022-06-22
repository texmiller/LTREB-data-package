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
library(patchwork) # for putting plots together

# library(gridExtra)
# library(grid)
# library(cowplot) # for pulling legend from ggplots
# library(cubelyr) # for working between lists of matrixes and dataframes


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
  dplyr::select(species,species_index, size_t) %>% 
  filter(!is.na(size_t)) %>% 
  group_by(species, species_index) %>% 
  summarise(actual_max_size = max(size_t),
            max_size = quantile(size_t,probs=0.975),
            max_size_99 = quantile(size_t,probs=0.99)) # The mean and sd effects plots look basically identical with either max size

# ggplot(data = LTREB_full)+
#   geom_histogram(aes(size_t)) +
#   geom_vline(data = max_size, aes(xintercept = max_size))+
#   geom_vline(data = max_size, aes(xintercept = max_size_99), col = "red")+
#   facet_wrap(~species, scales = "free")

#############################################################################################
####### Read in matrix population functions ------------------
#############################################################################################

source("Analyses/MPM_functions.R")


#############################################################################################
####### Read in Stan vital rate model outputs ------------------
#############################################################################################
tompath <- "C:/Users/tm9/Dropbox/EndodemogData/"
joshpath <- NULL
path<-tompath

surv_fit_seedling <- read_rds(paste0(path,"/Model_Runs/endo_seedling_surv.rds"))
surv_fit <- read_rds(paste0(path,"/Model_Runs/endo_spp_surv_woseedling.rds"))
grow_fit_seedling <- read_rds(paste0(path,"/Model_Runs/endo_seedling_grow_PIG_10000iterations.rds"))
grow_fit <- read_rds(paste0(path,"/Model_Runs/endo_spp_grow_PIG.rds"))
flw_fit <- read_rds(paste0(path,"/Model_Runs/endo_spp_flw.rds"))
fert_fit <- read_rds(paste0(path,"/Model_Runs/endo_spp_fert_PIG.rds"))
spike_fit <- read_rds(paste0(path,"/Model_Runs/endo_spp_spike_year_plot_nb.rds"))
seedmean_fit <- read_rds(paste0(path,"/Model_Runs/seed_mean.rds"))
stos_fit <- read_rds(paste0(path,"/Model_Runs/endo_spp_s_to_s.rds")) 

# surv_fit_seedling <- readRDS(url("https://www.dropbox.com/s/vf1mju5u4c4fs3t/endo_seedling_surv.rds?dl=1"))
# surv_fit <- readRDS(url("https://www.dropbox.com/s/00bor35inv5dypd/endo_spp_surv_woseedling.rds?dl=1"))
# grow_fit_seedling <- readRDS(url("https://www.dropbox.com/s/m0mw5z29slpm4p7/endo_seedling_grow_PIG_10000iterations.rds?dl=1"))
# grow_fit <- readRDS(url("https://www.dropbox.com/s/0ze8aooi9axj3oq/endo_spp_grow_PIG.rds?dl=1"))
# flw_fit <- readRDS(url("https://www.dropbox.com/s/ej65pn5k0km0z9c/endo_spp_flw.rds?dl=1"))
# fert_fit <- readRDS(url("https://www.dropbox.com/s/pk4x1j97kazu6pb/endo_spp_fert_pig.rds?dl=1"))
# spike_fit <- readRDS(url("https://www.dropbox.com/s/pjgui0n9tng6427/endo_spp_spike_year_plot_nb.rds?dl=1"))
# seedmean_fit <- readRDS(url("https://www.dropbox.com/s/3ma5yc8iusu8bh0/endo_spp_seed_mean.rds?dl=1"))
# stos_fit <- readRDS(url("https://www.dropbox.com/s/nf50hd76iw3hucw/endo_spp_s_to_s.rds?dl=1"))

surv_par <- rstan::extract(surv_fit, pars =quote_bare(beta0,betasize,betaendo,betaorigin,sigma_year,
                                                                   tau_year, tau_plot))
surv_sdlg_par <- rstan::extract(surv_fit_seedling, pars =quote_bare(beta0,betaendo,sigma_year,
                                                           tau_year, tau_plot))
grow_par <- rstan::extract(grow_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin,sigma_year,
                                                                    tau_year, tau_plot,
                                                       sigma))
grow_sdlg_par <- rstan::extract(grow_fit_seedling, pars = quote_bare(beta0,betaendo,sigma_year,
                                                       tau_year, tau_plot,
                                                       sigma))
flow_par <- rstan::extract(flw_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin,sigma_year,
                                                       tau_year, tau_plot))
fert_par <- rstan::extract(fert_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin,sigma_year,
                                                      tau_year, tau_plot))
spike_par <- rstan::extract(spike_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin,sigma_year,
                                                                         tau_year, tau_plot,
                                                         phi))
seed_par <- rstan::extract(seedmean_fit, pars = quote_bare(beta0,betaendo)) #no plot or year effect
recruit_par <- rstan::extract(stos_fit, pars = quote_bare(beta0,betaendo,sigma_year,tau_year, tau_plot))
## There are two possible ways of running the analysis:
## 1. Sample observation years (10-12 of these depending how you count)
## 2. Sample many hypothetical years from posterior of year effects

## SET-UP
# make the list of parameters and calculate mean lambdas
n_draws <- 100 # the means are the same whether we do 500 or 1000 draws
post_draws <- sample.int(7500,size=n_draws) # The models except for seedling growth have 7500 iterations. That one has more (15000 iterations) to help it converge.
n_spp <- length(unique(LTREB_full$species))
n_endo <- 2
# observation years
years_obs <- unique(LTREB_full$year_t1_index)
lambda_year_obs <- array(dim = c(length(years_obs),(n_spp+1),n_endo,n_draws))
# sampled years
n_years_samp <- 50
lambda_year_samp <- array(dim = c(n_years_samp,(n_spp+1),n_endo,n_draws))
# store mean and variance outputs for observed and sampled years
lambda_mean_obs <- lambda_mean_samp <- array(dim = c((n_spp+1),n_endo,n_draws))
lambda_sd_obs <- lambda_sd_samp <- array(dim = c((n_spp+1),n_endo,n_draws))
lambda_cv_obs <- lambda_cv_samp <- array(dim = c((n_spp+1),n_endo,n_draws))

for(i in 1:n_draws){
  for(e in 1:n_endo){
    for(s in 1:n_spp){
      for(y in years_obs){
        # 1. Sample observation years, calculate mean and variance
        lambda_year_obs[y,s,e,i] <- lambda(bigmatrix(make_params(species=s,
                                                             endo_mean=(e-1),
                                                             endo_var=(e-1),
                                                             original = 1, # should be =1 to represent recruit
                                                             draw=post_draws[i],
                                                             max_size=max_size,
                                                             rfx=T,
                                                             year=y, ## this pulls the t to t1 growth-surv transition and t1 reproduction -- see make_params()
                                                             surv_par=surv_par,
                                                             surv_sdlg_par = surv_sdlg_par,
                                                             grow_par=grow_par,
                                                             grow_sdlg_par = grow_sdlg_par,
                                                             flow_par=flow_par,
                                                             fert_par=fert_par,
                                                             spike_par=spike_par,
                                                             seed_par=seed_par,
                                                             recruit_par=recruit_par),
                                                 extension = 100)$MPMmat) # the extension parameter is used to fit the growth kernel to sizes larger than max size without losing probability density
      }
      lambda_mean_obs[s,e,i] <- mean(lambda_year_obs[,s,e,i])
      lambda_sd_obs[s,e,i] <- sd(lambda_year_obs[,s,e,i])
      lambda_cv_obs[s,e,i] <- lambda_mean_obs[s,e,i]/lambda_cv_obs[s,e,i]
      
      # 2.Sample many hypothetical years from fitted year variances  
      for(y in 1:n_years_samp){
        lambda_year_samp[y,s,e,i] <- lambda(bigmatrix(make_params(species=s,
                                                                  endo_mean=(e-1),
                                                                  endo_var=(e-1),
                                                                  original = 1, # should be =1 to represent recruit
                                                                  draw=post_draws[i],
                                                                  max_size=max_size,
                                                                  rfx=T,
                                                                  samp=T,
                                                                  year=y, ## this pulls the t to t1 growth-surv transition and t1 reproduction -- see make_params()
                                                                  surv_par=surv_par,
                                                                  surv_sdlg_par = surv_sdlg_par,
                                                                  grow_par=grow_par,
                                                                  grow_sdlg_par = grow_sdlg_par,
                                                                  flow_par=flow_par,
                                                                  fert_par=fert_par,
                                                                  spike_par=spike_par,
                                                                  seed_par=seed_par,
                                                                  recruit_par=recruit_par),
                                                      extension = 100)$MPMmat)
      }
      lambda_mean_samp[s,e,i] <- mean(lambda_year_samp[,s,e,i])
      lambda_sd_samp[s,e,i] <- sd(lambda_year_samp[,s,e,i])
      lambda_cv_samp[s,e,i] <- lambda_mean_samp[s,e,i]/lambda_sd_samp[s,e,i]
    }
    
    lambda_mean_obs[8,e,i] <- mean(lambda_year_obs[1:7,e,i])
    lambda_sd_obs[8,e,i] <- sd(lambda_year_obs[1:7,e,i])
    lambda_cv_obs[8,e,i] <- lambda_mean_obs[8,e,i]/lambda_cv_obs[8,e,i]
    
    lambda_mean_samp[8,e,i] <- mean(lambda_year_samp[1:7,e,i])
    lambda_sd_samp[8,e,i] <- sd(lambda_year_samp[1:7,e,i])
    lambda_cv_samp[8,e,i] <- lambda_mean_samp[8,e,i]/lambda_sd_samp[8,e,i]
  }
}


## there is something wrong with the way years are being indexed but I cannot deal with it right now
