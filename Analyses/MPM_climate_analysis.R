## Title: Grass endophyte population model with a bayesian framework
## Purpose: Builds matrix population model from vital rate estimates with SPEI effect
## Authors: Joshua and Tom
#############################################################


library(tidyverse)
library(scales)
library(bayesplot)
library(popbio)
# library(countreg)
library(actuar)
library(rstan)
library(patchwork) # for putting plots together


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

spei_range_df <- LTREB_full %>% 
  dplyr::select(species,species_index, spei12) %>% 
  filter(!is.na(spei12)) %>% 
  group_by(species, species_index) %>% 
  summarise(min_spei = min(spei12),
            max_spei = max(spei12))
  
#############################################################################################
####### Read in matrix population functions ------------------
#############################################################################################

source("Analyses/MPM_functions.R")



#############################################################################################
####### Read in Stan vital rate model outputs ------------------
#############################################################################################
tompath <- "C:/Users/tm9/Dropbox/EndodemogData/"
joshpath <- "~/Dropbox/EndodemogData/"
path<-joshpath

spei_surv_fit_seedling <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_seedling_surv_withoutquadraticterm.rds")
spei_surv_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_spp_surv_woseedling_linear.rds")
spei_grow_fit_seedling <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_seedling_grow_linear_10000iterations.rds")
spei_grow_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_grow_PIG_linear.rds")
spei_flw_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_spp_flw.rds")
spei_fert_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_fert_PIG_linear.rds")
spei_spike_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_spp_spike_linear.rds")
seedmean_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_seed_mean.rds") # doesn't include spei predictor
spei_stos_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_spp_s_to_s.rds") 



# Pulling out the actual parameters
#survival
spei_surv_par <- rstan::extract(spei_surv_fit, pars =quote_bare(beta0,betasize,betaendo,betaorigin,
                                                                betaspei_endo,
                                                                tau_year, tau_plot, sigma_year))
# seedling survival
spei_surv_sdlg_par <- rstan::extract(spei_surv_fit_seedling, pars =quote_bare(beta0,betaendo, betaspei_endo, tau_year, tau_plot, sigma_year))
# growth
spei_grow_par <- rstan::extract(spei_grow_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin, betaspei_endo, tau_year, tau_plot, sigma, sigma_year))
# seedling growth
spei_grow_sdlg_par <- rstan::extract(spei_grow_fit_seedling, pars = quote_bare(beta0,betaendo, betaspei_endo, tau_year, tau_plot, sigma, sigma_year))
#flowering
spei_flow_par <- rstan::extract(spei_flw_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin, betaspei_endo, tau_year, tau_plot, sigma_year))
# number of flower tillers
spei_fert_par <- rstan::extract(spei_fert_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin, betaspei_endo, tau_year, tau_plot, sigma_year))
# spikelets/infl
spei_spike_par <- rstan::extract(spei_spike_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin, betaspei_endo, tau_year, tau_plot, sigma_year))
seed_par <- rstan::extract(seedmean_fit, pars = quote_bare(beta0,betaendo)) #no plot or year effect
# seedling recruitment
spei_recruit_par <- rstan::extract(spei_stos_fit, pars = quote_bare(beta0,betaendo, betaspei_endo, tau_year, tau_plot, sigma_year))


# endophyte effects on lambda mean across spei values ---------------------------


## SET-UP
# make the list of parameters and calculate mean lambdas
n_draws <- 10 # the means are the same whether we do 500 or 1000 draws
post_draws <- sample.int(7500,size=n_draws) # The models except for seedling growth have 7500 iterations. That one has more (15000 iterations) to help it converge.
n_spp <- length(unique(LTREB_full$species))
n_endo <- 2
# stepsize of spei values
spei_steps <- 10
spei_range <- array(dim = c(n_spp,spei_steps))
for(s in 1:n_spp){
spei_range[s,] <- seq(from = spei_range_df$min_spei[s], to = spei_range_df$max_spei[s], length.out = 10)
}
lambda_spei <- array(dim = c(spei_steps,(n_spp+1),n_endo,n_draws))


for(i in 1:n_draws){
  for(e in 1:n_endo){
    for(s in 1:n_spp){
      for(c in 1:spei_steps){
        # 1. Sample observation years, calculate mean and variance
        lambda_spei[c,s,e,i] <- lambda(bigmatrix(make_params(species=s,
                                                                 endo_mean=(e-1),
                                                                 original = 1, # should be =1 to represent recruit
                                                                 draw=post_draws[i],
                                                                 max_size=max_size,
                                                                 rfx=F,
                                                                 spei = T,
                                                                 surv_par=spei_surv_par,
                                                                 surv_sdlg_par = spei_surv_sdlg_par,
                                                                 grow_par=spei_grow_par,
                                                                 grow_sdlg_par = spei_grow_sdlg_par,
                                                                 flow_par=spei_flow_par,
                                                                 fert_par=spei_fert_par,
                                                                 spike_par=spei_spike_par,
                                                                 seed_par=seed_par,
                                                                 recruit_par=spei_recruit_par),
                                                     spei = spei_range[s,c],
                                                     extension = 100)$MPMmat) # the extension parameter is used to fit the growth kernel to sizes larger than max size without losing probability density
      } # endo of spei loop
    } # endo of species loop
    # calculating overall species mean
    lambda_spei[,8,e,i] <- mean(lambda_spei[,1:7,e,i])
  }# endo of endo loop
} # end of iteration loop

dimnames(lambda_spei) <- list(spei = paste0("spei",1:spei_steps),species = paste0("s",1:nspp+1), Endo = paste0("e",1:2), Iteration= paste0("i",1:n_draws))
lambda_spei_cube <- cubelyr::as.tbl_cube(lambda_spei)

lambda_spei_df <- as_tibble(lambda_spei_cube) %>% 
  pivot_wider(names_from = Endo, values_from = lambda_spei) %>% 
  mutate(lambda_diff = e2-e1) %>% 
  mutate(Species = case_when(species == "s1" ~ "Agrostis perennans",
                             species == "s2" ~ "Elymus villosus",
                             species == "s3" ~ "Elymus virginicus",
                             species == "s4" ~ "Festuca subverticillata",
                             species == "s5" ~ "Lolium arundinaceum",
                             species == "s6" ~ "Poa alsodes",
                             species == "s7" ~ "Poa sylvestris",
                             species == "s8" ~ "Species Mean")) 


