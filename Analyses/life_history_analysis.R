## Title: Grass endophyte life history metrics
## Purpose: calculates life histoory metrics from matrix population model
## and compares them to endophyte effects on mean and variance
## Authors: Joshua and Tom
#############################################################

library(tidyverse)
library(cubelyr)
library(scales)
library(popbio) # for MPM
library(Rage) # for gen time, longevity etc.
# library(countreg)
library(actuar)
library(rstan)
library(brms)

library(ape) # for phylogenetic contrasts

library(patchwork)
library(bbmle) # for AIC

quote_bare <- function( ... ){
  substitute( alist(...) ) %>% 
    eval( ) %>% 
    sapply( deparse )
}

#############################################################################################
####### Read in Data and creating size bins------------------
#############################################################################################

# source("Analyses/endodemog_data_processing.R")

tompath <- "C:/Users/tm9/Dropbox/EndodemogData/"
joshpath <- "~/Dropbox/EndodemogData/"
path<-joshpath

LTREB_full <- read_csv(paste0(path,"Fulldataplusmetadata/LTREB_full.csv"))

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

surv_par <- rstan::extract(surv_fit, pars =quote_bare(beta0,betasize,betaendo,betaorigin,
                                                      tau_year, tau_plot))
surv_sdlg_par <- rstan::extract(surv_fit_seedling, pars =quote_bare(beta0,betaendo,
                                                                    tau_year, tau_plot))
grow_par <- rstan::extract(grow_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin,
                                                       tau_year, tau_plot,
                                                       sigma))
grow_sdlg_par <- rstan::extract(grow_fit_seedling, pars = quote_bare(beta0,betaendo,
                                                                     tau_year, tau_plot,
                                                                     sigma))
flow_par <- rstan::extract(flw_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin,
                                                      tau_year, tau_plot))
fert_par <- rstan::extract(fert_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin,
                                                       tau_year, tau_plot))
spike_par <- rstan::extract(spike_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin,
                                                         tau_year, tau_plot,
                                                         phi))
seed_par <- rstan::extract(seedmean_fit, pars = quote_bare(beta0,betaendo)) #no plot or year effect
recruit_par <- rstan::extract(stos_fit, pars = quote_bare(beta0,betaendo,
                                                          tau_year, tau_plot))
# dim(surv_par$tau_year)
# plot(surv_par$tau_year[,1,1,], grow_par$tau_year[,1,1,], col = levels(as_factor(grow_par$tau_year[1,1,1,])))
#############################################################################################
####### Run the MPM ------------------
#############################################################################################

# make the list of parameters and calculate mean lambdas
n_draws <- 500 # the means are the same whether we do 500 or 1000 draws
post_draws <- sample.int(7500,size=n_draws) # The models except for seedling growth have 7500 iterations. That one has more (15000 iterations) to help it converge.

gen_time <- array(dim = c(7,2,n_draws))
longev <- array(dim = c(7,2,n_draws))
mean_life_expect <- array(dim = c(7,2,n_draws))
var_life_expect <- array(dim = c(7,2,n_draws))
R0 <- array(dim = c(7,2,n_draws))
repro_age <- array(dim = c(7,2,n_draws))
for(i in 1:length(post_draws)){
  for(e in 1:2){
    for(s in 1:7){
      gen_time[s,e,i] <- gen_time(matU = bigmatrix(make_params(species=s,
                                                               endo_mean=(e-1),
                                                               endo_var=(e-1),
                                                               original = 0, # should be =1 to represent recruit
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
                                                               recruit_par=recruit_par), 
                                                   extension = 100)$Tmat,
                                  matR = bigmatrix(make_params(species=s,
                                                               endo_mean=(e-1),
                                                               endo_var=(e-1),
                                                               original = 0, # should be =1 to represent recruit
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
                                                               recruit_par=recruit_par), 
                                                   extension = 100)$Fmat)
      longev[s,e,i] <- longevity(matU = bigmatrix(make_params(species=s,
                                                                        endo_mean=(e-1),
                                                                        endo_var=(e-1),
                                                                        original = 0, # should be =1 to represent recruit
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
                                                                        recruit_par=recruit_par), 
                                                            extension = 100)$Tmat, start = 1, lx_crit = 0.05)
      mean_life_expect[s,e,i] <- life_expect_mean(matU = bigmatrix(make_params(species=s,
                                                              endo_mean=(e-1),
                                                              endo_var=(e-1),
                                                              original = 0, # should be =1 to represent recruit
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
                                                              recruit_par=recruit_par), 
                                                  extension = 100)$Tmat, start = 1)
      var_life_expect[s,e,i] <- life_expect_var(matU = bigmatrix(make_params(species=s,
                                                                               endo_mean=(e-1),
                                                                               endo_var=(e-1),
                                                                               original = 0, # should be =1 to represent recruit
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
                                                                               recruit_par=recruit_par), 
                                                                   extension = 100)$Tmat, start = 1)
      R0[s,e,i] <- net_repro_rate(matU = bigmatrix(make_params(species=s,
                                                                         endo_mean=(e-1),
                                                                         endo_var=(e-1),
                                                                         original = 0, # should be =1 to represent recruit
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
                                                                         recruit_par=recruit_par), 
                                                             extension = 100)$Tmat,
                                            matR = bigmatrix(make_params(species=s,
                                                                         endo_mean=(e-1),
                                                                         endo_var=(e-1),
                                                                         original = 0, # should be =1 to represent recruit
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
                                                                         recruit_par=recruit_par), 
                                                             extension = 100)$Fmat)
      repro_age[s,e,i] <- mature_age(matU = bigmatrix(make_params(species=s,
                                                               endo_mean=(e-1),
                                                               endo_var=(e-1),
                                                               original = 0, # should be =1 to represent recruit
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
                                                               recruit_par=recruit_par), 
                                                   extension = 100)$Tmat,
                                  matR = bigmatrix(make_params(species=s,
                                                               endo_mean=(e-1),
                                                               endo_var=(e-1),
                                                               original = 0, # should be =1 to represent recruit
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
                                                               recruit_par=recruit_par), 
                                                   extension = 100)$Fmat)
    }
  }
}

#saving the generation time calculation and calculating mean for each species
saveRDS(gen_time, file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/gen_time.rds")
saveRDS(longev, file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/longev.rds")
saveRDS(mean_life_expect, file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/mean_life_expect.rds")
saveRDS(var_life_expect, file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/var_life_expect.rds")
saveRDS(R0, file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/R0.rds")
saveRDS(repro_age, file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/repro_age.rds")

gen_time <- read_rds(file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/gen_time.rds")
longev <- read_rds(file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/longev.rds")
mean_life_expect <- read_rds(file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/mean_life_expect.rds")
var_life_expect <- read_rds(file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/var_life_expect.rds")
R0 <- read_rds(file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/R0.rds")
# repro_age <- read_rds(repro_age, file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/repro_age.rds")

gen_time_summary <- longev_summary <- mean_life_expect_summary <- var_life_expect_summary <-  R0_summary <- repro_age_summary <- matrix(NA,7,4)
for(s in 1:7){
  gen_time_summary[s,1] <- mean(gen_time[s,1,])
  gen_time_summary[s,2] <- mean(gen_time[s,2,])
  gen_time_summary[s,3] <- mean(gen_time[s,,])
  gen_time_summary[s,4] <- sd(gen_time[s,1,])
  
  longev_summary[s,1] <- mean(longev[s,1,])
  longev_summary[s,2] <- mean(longev[s,2,])
  longev_summary[s,3] <- mean(longev[s,,])
  longev_summary[s,4] <- sd(longev[s,1,])
  
  
  mean_life_expect_summary[s,1] <- mean(mean_life_expect[s,1,])
  mean_life_expect_summary[s,2] <- mean(mean_life_expect[s,2,])
  mean_life_expect_summary[s,3] <- mean(mean_life_expect[s,,])
  mean_life_expect_summary[s,4] <- sd(mean_life_expect[s,1,])
  
  
  var_life_expect_summary[s,1] <- mean(var_life_expect[s,1,])
  var_life_expect_summary[s,2] <- mean(var_life_expect[s,2,])
  var_life_expect_summary[s,3] <- mean(var_life_expect[s,,])
  
  R0_summary[s,1] <- mean(R0[s,1,])
  R0_summary[s,2] <- mean(R0[s,2,])
  R0_summary[s,3] <- mean(R0[s,,])
  R0_summary[s,4] <- sd(R0[s,1,])
  
  
  # repro_age_summary[s,1] <- mean(repro_age[s,1,])
  # repro_age_summary[s,2] <- mean(repro_age[s,2,])
  # repro_age_summary[s,3] <- mean(repro_age[s,,])
}


# reading in lambda_mean and lambda_var with 500 post draws from dropbox, derived from MPM_analysis script
lambda_mean <- read_rds(file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/lambda_mean.rds")

lambda_hold <- read_rds(file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/lambda_hold.rds")
lambda_var <- read_rds(file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/lambda_var.rds")

# Mean endophyte difference and quantiles
lambda_means <- matrix(NA,8,2)
lambda_mean_diff <- matrix(NA,8,8)
for(s in 1:8){
  lambda_means[s,1] <- mean(lambda_mean[s,1,])
  lambda_means[s,2] <- mean(lambda_mean[s,2,])
  lambda_mean_diff[s,1] = mean(lambda_mean[s,2,] - lambda_mean[s,1,])
  lambda_mean_diff[s,2:7] = quantile(lambda_mean[s,2,] - lambda_mean[s,1,],probs=c(0.05,0.125,0.25,0.75,0.875,0.95))
  lambda_mean_diff[s,8] = sd(lambda_mean[s,2,] - lambda_mean[s,1,])
  
}

# Calculating endophyte effect on sd and variance
lambda_sds <- matrix(NA,8,2)
lambda_vars <- matrix(NA,8,2)

lambda_cv <- (lambda_var^2)/(lambda_mean) 
lambda_cvs <- matrix(NA,8,2)

lambda_sd_diff <- matrix(NA,8,7)
lambda_var_diff <- matrix(NA,8,7)
lambda_cv_diff <-  matrix(NA,8,8)
for(s in 1:8){
  lambda_sds[s,1] <- mean(lambda_var[s,1,])
  lambda_sds[s,2] <- mean(lambda_var[s,2,])
  
  lambda_vars[s,1] <- mean(lambda_var[s,1,])^2
  lambda_vars[s,2] <- mean(lambda_var[s,2,])^2
  # 
  lambda_cvs[s,1] <- mean(lambda_cv[s,1,])
  lambda_cvs[s,2] <- mean(lambda_cv[s,2,])
  
  lambda_sd_diff[s,1] = mean(lambda_var[s,2,]) - mean(lambda_var[s,1,])
  lambda_sd_diff[s,2:7] = quantile(lambda_var[s,2,] - lambda_var[s,1,],probs=c(0.05,0.125,0.25,0.75,0.875,0.95))
  
  lambda_var_diff[s,1] = mean(lambda_var[s,2,]^2 - lambda_var[s,1,]^2)
  lambda_var_diff[s,2:7] = quantile(lambda_var[s,2,]^2 - lambda_var[s,1,]^2,probs=c(0.05,0.125,0.25,0.75,0.875,0.95))
  
  lambda_cv_diff[s,1] = mean(lambda_cv[s,2,] - lambda_cv[s,1,])
  lambda_cv_diff[s,2:7] = quantile(lambda_cv[s,2,] - lambda_cv[s,1,],probs=c(0.05,0.125,0.25,0.75,0.875,0.95))
  lambda_cv_diff[s,8] = sd(lambda_cv[s,2,] - lambda_cv[s,1,])
  
}
#seed length measurements taken from Rudgers et al. 2009. from Flora of North America
seed_size <- c(1.75,7.25,8,3.75,7,3.45,2.6)
#imperfect transmission measurements from LTREB proposal, filling in LOAR with 100 for now, and POSY rate seems really low though?
imperfect_trans <- c(69.8,100,100,42.7,100,99.9,16.6)
# imperfect_trans <- c(69.8,100,100,42.7,100,99.9,99)

# count of branches in phylogeny from J of Sytematics Evolution, Volume: 53, Issue: 2, Pages: 117-137, First published: 02 March 2015, DOI: (10.1111/jse.12150) 
# relatedness = c(7,5,5,10,7,7,10)
subtribes <- c("Agrostinae", "Tritaceae","Tritaceae", "Loliinae", "Loliinae", "Poinae", "Poinae")


empirical_ages <- LTREB_full %>% 
  dplyr::select(species,species_index, pos, id) %>% 
  # filter(!is.na(size_t)) %>% 
  group_by(species, species_index, id) %>% 
  summarise(age = n()) %>% 
  group_by(species, species_index) %>% 
  summarize(observed_max_age = max(age),
            max_age_97.5 = quantile(age,probs=0.975),
            max_age_99 = quantile(age,probs=0.99))

traits_df <- empirical_ages
traits_df$gen_time <- gen_time_summary[,1]
traits_df$sd_gen_time <- gen_time_summary[,4]
traits_df$longev <- longev_summary[,1]
traits_df$sd_longev <- longev_summary[,4]

traits_df$mean_life_expect <- mean_life_expect_summary[,1]
traits_df$var_life_expect <- var_life_expect_summary[,1]
traits_df$R0 <- R0_summary[,1]
# traits_df$repro_age <- repro_age_summary[,3]
traits_df$seed_size <- seed_size
traits_df$imperfect_trans <- imperfect_trans

traits_df$sd_effect <- lambda_sd_diff[1:7,1]
traits_df$cv_effect <- lambda_cv_diff[1:7,1]
traits_df$sdof_cv_effect <- lambda_cv_diff[1:7,8]
traits_df$mean_effect <- lambda_mean_diff[1:7,1]
traits_df$sdof_mean_effect <- lambda_mean_diff[1:7,8]

traits_df_long <- traits_df %>% 
  pivot_longer(cols = c("observed_max_age", "max_age_97.5", "max_age_99", "gen_time", 
                        "longev", "mean_life_expect", "var_life_expect", 
                        "R0", "seed_size", "imperfect_trans")) %>% 
  mutate(name_label = case_when(name == "observed_max_age" ~ "Obs. Max Age", name == "max_age_99" ~ "99th Percentile Max Age", name == "R0" ~ "R0",
                                name == "longev" ~ "Longevity", name == "mean_life_expect" ~ "Mean Life Expectancy",
                                name == "gen_time" ~ "Generation Time", name == "seed_size" ~ "Seed Length (mm)", name == "imperfect_trans" ~ "Imperfect Transmission Rate")) 


ggplot(data = traits_df_long)+
  geom_point(aes(y = sd_effect, x = value))+
  geom_smooth(aes(y = sd_effect, x = value), method = "glm")+
  facet_wrap(~name, scales = "free_x")

ggplot(data = filter(traits_df_long, species != "LOAR"))+
  geom_point(aes(y = sd_effect, x = value))+
  geom_smooth(aes(y = sd_effect, x = value), method = "glm")+
  facet_wrap(~name, scales = "free_x")

ggplot(data = traits_df_long)+
  geom_point(aes(y = cv_effect, x = value))+
  geom_smooth(aes(y = cv_effect, x = value), method = "glm")+
  facet_wrap(~name, scales = "free")

ggplot(data = traits_df_long)+
  geom_point(aes(y = mean_effect, x = value))+
  geom_smooth(aes(y = mean_effect, x = value), method = "glm")+
  facet_wrap(~name, scales = "free")

ggplot(data = traits_df)+
  geom_point(aes(x = mean_life_expect, y = var_life_expect))


############################################################################################
#####  Fitting regressions of the traits on the mean and variance effects        ###########
############################################################################################

# model selection
lh_models <- list()
# lh_models[[1]] <- lm(sd_effect ~ observed_max_age, data = traits_df)
# lh_models[[2]] <- lm(sd_effect ~ max_age_97.5, data = traits_df)
# lh_models[[3]] <- lm(sd_effect ~ max_age_99, data = traits_df)
# lh_models[[4]] <- lm(sd_effect ~ gen_time, data = traits_df)
# lh_models[[5]] <- lm(sd_effect ~ longev, data = traits_df)
# lh_models[[6]] <- lm(sd_effect ~ mean_life_expect, data = traits_df)
# lh_models[[7]] <- lm(sd_effect ~ var_life_expect, data = traits_df)
# lh_models[[8]] <- lm(sd_effect ~ R0, data = traits_df)
# lh_models[[9]] <- lm(sd_effect ~ seed_size, data = traits_df)
# lh_models[[10]] <- lm(sd_effect ~ imperfect_trans, data = traits_df)
#
# lh_models[[11]] <- lm(mean_effect ~ observed_max_age, data = traits_df)
# lh_models[[12]] <- lm(mean_effect ~ max_age_97.5, data = traits_df)
# lh_models[[13]] <- lm(mean_effect ~ max_age_99, data = traits_df)
# lh_models[[14]] <- lm(mean_effect ~ gen_time, data = traits_df)
# lh_models[[15]] <- lm(mean_effect ~ longev, data = traits_df)
# lh_models[[16]] <- lm(mean_effect ~ mean_life_expect, data = traits_df)
# lh_models[[17]] <- lm(mean_effect ~ var_life_expect, data = traits_df)
# lh_models[[18]] <- lm(mean_effect ~ R0, data = traits_df)
# lh_models[[19]] <- lm(mean_effect ~ seed_size, data = traits_df)
# lh_models[[20]] <- lm(mean_effect ~ imperfect_trans, data = traits_df)

summary(lh_models[[1]] <- lm(cv_effect ~ observed_max_age, data = traits_df))
summary(lh_models[[2]] <- lm(cv_effect ~ max_age_99, data = traits_df))
summary(lh_models[[3]] <- lm(cv_effect ~ gen_time, data = traits_df))
summary(lh_models[[4]] <- lm(cv_effect ~ longev, data = traits_df))
summary(lh_models[[5]] <- lm(cv_effect ~ mean_life_expect, data = traits_df))
summary(lh_models[[6]] <- lm(cv_effect ~ R0, data = traits_df))
summary(lh_models[[7]] <- lm(cv_effect ~ seed_size, data = traits_df))
summary(lh_models[[8]] <- lm(cv_effect ~ imperfect_trans, data = traits_df))
# 
# AICtab(lh_models)
# summary(lh_models[[1]]) #*
# summary(lh_models[[2]])#*
# summary(lh_models[[3]])
# summary(lh_models[[4]])#.06
# summary(lh_models[[5]])#.09
# summary(lh_models[[6]])#*
# summary(lh_models[[7]])#.05
# summary(lh_models[[8]])
# 
# # making data ranges for predictions
# newdata <- data.frame(observed_max_age = seq(from = min(traits_df$observed_max_age), to = max(traits_df$observed_max_age), length.out = 10),
#                       max_age_99 = seq(from = min(traits_df$max_age_99), to = max(traits_df$max_age_99), length.out = 10),
#                       gen_time = seq(from = min(traits_df$gen_time), to = max(traits_df$gen_time), length.out = 10),
#                       longev = seq(from = min(traits_df$longev), to = max(traits_df$longev), length.out = 10),
#                       mean_life_expect = seq(from = min(traits_df$mean_life_expect), to = max(traits_df$mean_life_expect), length.out = 10),
#                       R0 = seq(from = min(traits_df$R0), to = max(traits_df$R0), length.out = 10),
#                       seed_size = seq(from = min(traits_df$seed_size), to = max(traits_df$seed_size), length.out = 10),
#                       imperfect_trans = seq(from = min(traits_df$imperfect_trans), to = max(traits_df$imperfect_trans), length.out = 10)) 
# newdata_fit <- newdata %>%
#   rename_with(.fn = str_c, pattern = ".x") %>% 
#   mutate(row_id = row_number()) %>% 
#   mutate(observed_max_age.fit = predict(lh_models[[1]], newdata, type = "response"),
#          max_age_99.fit = predict(lh_models[[2]], newdata = newdata,  type = "response"),
#          gen_time.fit = predict(lh_models[[3]], newdata = newdata,  type = "response"),
#          longev.fit = predict(lh_models[[4]], newdata = newdata,   type = "response"),
#          mean_life_expect.fit = predict(lh_models[[5]], newdata = newdata,  type = "response"),
#          R0.fit =predict(lh_models[[6]], newdata = newdata,  type = "response"),
#          seed_size.fit =predict(lh_models[[7]], newdata = newdata,  type = "response"),
#          imperfect_trans.fit =predict(lh_models[[8]], newdata = newdata,  type = "response")) %>%
#   mutate(observed_max_age.lwr = predict(lh_models[[1]], newdata = newdata,  interval = "confidence", type = "response")[,2],
#          max_age_99.lwr = predict(lh_models[[2]], newdata = newdata, interval = "confidence", type = "response")[,2],
#          gen_time.lwr = predict(lh_models[[3]], newdata = newdata, interval = "confidence", type = "response")[,2],
#          longev.lwr = predict(lh_models[[4]], newdata = newdata,  interval = "confidence", type = "response")[,2],
#          mean_life_expect.lwr = predict(lh_models[[5]], newdata = newdata, interval = "confidence", type = "response")[,2],
#          R0.lwr =predict(lh_models[[6]], newdata = newdata, interval = "confidence", type = "response")[,2],
#          seed_size.lwr =predict(lh_models[[7]], newdata = newdata, interval = "confidence", type = "response")[,2],
#          imperfect_trans.lwr =predict(lh_models[[8]], newdata = newdata, interval = "confidence", type = "response")[,2]) %>%
#   mutate(observed_max_age.upr = predict(lh_models[[1]], newdata = newdata,  interval = "confidence", type = "response")[,3],
#          max_age_99.upr = predict(lh_models[[2]], newdata = newdata, interval = "confidence", type = "response")[,3],
#          gen_time.upr = predict(lh_models[[3]], newdata = newdata, interval = "confidence", type = "response")[,3],
#          longev.upr = predict(lh_models[[4]], newdata = newdata,  interval = "confidence", type = "response")[,3],
#          mean_life_expect.upr = predict(lh_models[[5]], newdata = newdata, interval = "confidence", type = "response")[,3],
#          R0.upr =predict(lh_models[[6]], newdata = newdata, interval = "confidence", type = "response")[,3],
#          seed_size.upr =predict(lh_models[[7]], newdata = newdata, interval = "confidence", type = "response")[,3],
#          imperfect_trans.upr =predict(lh_models[[8]], newdata = newdata, interval = "confidence", type = "response")[,3]) %>%
#   pivot_longer(cols = -row_id, names_to = c("name", "interval"), names_sep = "\\.") %>% 
#   pivot_wider(id_cols = c(row_id,name), names_from = interval, values_from = value) %>% 
#   mutate(significance = case_when(name == "observed_max_age" | name == "max_age_99" | name == "R0" ~ "<.05",
#                                   name == "longev" | name == "mean_life_expect" | name == "seed_size" ~ "<.1",
#                                   TRUE ~ ">.1")) %>% 
#   mutate(name_label = case_when(name == "observed_max_age" ~ "Obs. Max Age", name == "max_age_99" ~ "99th Percentile Max Age", name == "R0" ~ "R0",
#          name == "longev" ~ "Longevity", name == "mean_life_expect" ~ "Mean Life Expectancy",
#          name == "gen_time" ~ "Generation Time", name == "seed_size" ~ "Seed Length (mm)", name == "imperfect_trans" ~ "Imperfect Transmission Rate")) 
# 
# lh_facetplot <- ggplot()+
#   geom_ribbon(data = newdata_fit, aes(ymin = lwr, ymax = upr, x = x, alpha = significance))+
#   geom_line(data = newdata_fit, aes(y = fit, x = x, alpha = significance))+
#   geom_point(data = filter(traits_df_long, name !="var_life_expect" & name != "repro_age" & name != "max_age_97.5"), aes(y = cv_effect, x = value))+
#   facet_wrap(~factor(name_label, levels = c("Obs. Max Age","99th Percentile Max Age","R0","Longevity", "Mean Life Expectancy","Generation Time","Seed Length (mm)","Imperfect Transmission Rate")),nrow = 2, scales = "free", strip.position = "bottom", labeller = labeller(name_label = label_wrap_gen(10)))+
#   scale_alpha_manual(values = c(.7,.3,.1))+
#   theme_classic()+
#   theme(strip.background = element_blank(),
#         strip.placement = "outside")+
#   labs(x = "", y = expression(paste("Effect on CV", (lambda))), alpha = "p-value")
# lh_facetplot
# 
# lh_legend <- cowplot::get_legend(lh_facetplot) 
# ggsave(lh_legend, filename = "lh_legend.png")
# #making the facets into separate plots to annotate
# 
# oba_plot <- lh_plot <- ggplot()+
#   geom_ribbon(data = filter(newdata_fit, name == "observed_max_age"), aes(ymin = lwr, ymax = upr, x = x), alpha = .7)+
#   geom_line(data = filter(newdata_fit, name == "observed_max_age"), aes(y = fit, x = x))+
#   geom_point(data = filter(traits_df_long, name == "observed_max_age"), aes(y = cv_effect, x = value))+
#   geom_text(aes(x = Inf, y = Inf), label = "**",  hjust   = 2, vjust   = 1, size = 3)+
#   theme_classic()+
#   theme(strip.background = element_blank(),
#         strip.placement = "outside",
#         axis.title.x = element_text(size=10))+
#   labs(x = "Obs. Max Age", y = expression(paste("Effect on CV", (lambda))), alpha = "p-value")
# # oba_plot
# 
# ma99_plot <- lh_plot <- ggplot()+
#   geom_ribbon(data = filter(newdata_fit, name == "max_age_99"), aes(ymin = lwr, ymax = upr, x = x), alpha = .7)+
#   geom_line(data = filter(newdata_fit, name == "max_age_99"), aes(y = fit, x = x))+
#   geom_point(data = filter(traits_df_long, name == "max_age_99"), aes(y = cv_effect, x = value))+
#   geom_text(aes(x = Inf, y = Inf), label = "**",  hjust   = 2, vjust   = 1, size = 3)+
#   theme_classic()+
#   theme(strip.background = element_blank(),
#         strip.placement = "outside",
#         axis.title.x = element_text(size=10))+
#   labs(x = "99th Percentile Max Age", y = "", alpha = "p-value")
# # ma99_plot
# 
# R0_plot <- lh_plot <- ggplot()+
#   geom_ribbon(data = filter(newdata_fit, name == "R0"), aes(ymin = lwr, ymax = upr, x = x), alpha = .7)+
#   geom_line(data = filter(newdata_fit, name == "R0"), aes(y = fit, x = x))+
#   geom_point(data = filter(traits_df_long, name == "R0"), aes(y = cv_effect, x = value))+
#   geom_text(aes(x = Inf, y = Inf), label = "**",  hjust   = 2, vjust   = 1, size = 3)+
#   theme_classic()+
#   theme(strip.background = element_blank(),
#         strip.placement = "outside",
#         axis.title.x = element_text(size=10))+
#   labs(x = "R0", y = "", alpha = "p-value")
# # R0_plot
# 
# longev_plot <- lh_plot <- ggplot()+
#   geom_ribbon(data = filter(newdata_fit, name == "longev"), aes(ymin = lwr, ymax = upr, x = x), alpha = .3)+
#   geom_line(data = filter(newdata_fit, name == "longev"), aes(y = fit, x = x), alpha = .4)+
#   geom_point(data = filter(traits_df_long, name == "longev"), aes(y = cv_effect, x = value))+
#   geom_text(aes(x = Inf, y = Inf), label = "*",  hjust   = 2, vjust   = 1, size = 3)+
#   theme_classic()+
#   theme(strip.background = element_blank(),
#         strip.placement = "outside",
#         axis.title.x = element_text(size=10))+
#   labs(x = "Longevity", y = "", alpha = "p-value")
# # longev_plot
# 
# mle_plot <- lh_plot <- ggplot()+
#   geom_ribbon(data = filter(newdata_fit, name == "mean_life_expect"), aes(ymin = lwr, ymax = upr, x = x), alpha = .3)+
#   geom_line(data = filter(newdata_fit, name == "mean_life_expect"), aes(y = fit, x = x), alpha = .4)+
#   geom_point(data = filter(traits_df_long, name == "mean_life_expect"), aes(y = cv_effect, x = value))+
#   geom_text(aes(x = Inf, y = Inf), label = "*",  hjust   = 2, vjust   = 1, size = 3)+
#   theme_classic()+
#   theme(strip.background = element_blank(),
#         strip.placement = "outside",
#         axis.title.x = element_text(size=10))+
#   labs(x = "Mean Life Expectancy", y = expression(paste("Effect on CV", (lambda))), alpha = "p-value")
# # mle_plot
# 
# gt_plot <- lh_plot <- ggplot()+
#   geom_ribbon(data = filter(newdata_fit, name == "gen_time"), aes(ymin = lwr, ymax = upr, x = x), alpha = .1)+
#   geom_line(data = filter(newdata_fit, name == "gen_time"), aes(y = fit, x = x), alpha = .2)+
#   geom_point(data = filter(traits_df_long, name == "gen_time"), aes(y = cv_effect, x = value))+
#   geom_text(aes(x = Inf, y = Inf), label = "",  hjust   = 2, vjust   = 1, size = 3)+
#   scale_alpha_manual(values = c(.1))+
#   theme_classic()+
#   theme(strip.background = element_blank(),
#         strip.placement = "outside",
#         axis.title.x = element_text(size=10))+
#   labs(x = "Generation Time", y = "", alpha = "p-value")
# # gt_plot
# 
# ss_plot <- lh_plot <- ggplot()+
#   geom_ribbon(data = filter(newdata_fit, name == "seed_size"), aes(ymin = lwr, ymax = upr, x = x), alpha = .1)+
#   geom_line(data = filter(newdata_fit, name == "seed_size"), aes(y = fit, x = x),alpha =.2)+
#   geom_point(data = filter(traits_df_long, name == "seed_size"), aes(y = cv_effect, x = value))+
#   geom_text(aes(x = Inf, y = Inf), label = "",  hjust   = 2, vjust   = 1, size = 3)+
#   theme_classic()+
#   theme(strip.background = element_blank(),
#         strip.placement = "outside",
#         axis.title.x = element_text(size=10))+
#   labs(x = "Seed Length (mm)", y = "", alpha = "p-value")
# # ss_plot
# 
# it_plot <- lh_plot <- ggplot()+
#   geom_ribbon(data = filter(newdata_fit, name == "imperfect_trans"), aes(ymin = lwr, ymax = upr, x = x), alpha = .1)+
#   geom_line(data = filter(newdata_fit, name == "imperfect_trans"), aes(y = fit, x = x), alpha = .2)+
#   geom_point(data = filter(traits_df_long, name == "imperfect_trans"), aes(y = cv_effect, x = value))+
#   geom_text(aes(x = Inf, y = Inf), label = "",  hjust   = 2, vjust   = 1, size = 3)+
#   theme_classic()+
#   theme(strip.background = element_blank(),
#         strip.placement = "outside",
#         axis.title.x = element_text(size=10))+
#   labs(x = "Imperfect Transmission Rate", y = "", alpha = "p-value")
# # it_plot
# 
# lh_plot <- oba_plot + ma99_plot + R0_plot + longev_plot + mle_plot + gt_plot + ss_plot + it_plot + 
#   plot_layout(nrow = 2) + plot_annotation(tag_levels = "A") 
# ggsave(lh_plot, filename = "lh_plot.png", width = 10, height = 5)
# 
# 

# brms version
## run this code to optimize computer system settings for MCMC
rstan_options(auto_write = FALSE)
options(mc.cores = parallel::detectCores())
set.seed(123)


lh_models <- list()
lh_models[[1]] <- brm(cv_effect ~ observed_max_age, data = traits_df,
                      family = "gaussian",
                      prior = c(set_prior("normal(0,5)", class = "b")))
lh_models[[2]] <- brm(cv_effect ~ max_age_99, data = traits_df,
                      family = "gaussian",
                      prior = c(set_prior("normal(0,5)", class = "b")))
lh_models[[3]] <- brm(cv_effect ~ R0, data = traits_df,
                      family = "gaussian",
                      prior = c(set_prior("normal(0,5)", class = "b")))
lh_models[[4]] <- brm(cv_effect ~ longev, data = traits_df,
                      family = "gaussian",
                      prior = c(set_prior("normal(0,5)", class = "b")))
lh_models[[5]] <- brm(cv_effect ~ mean_life_expect, data = traits_df,
                      family = "gaussian",
                      prior = c(set_prior("normal(0,5)", class = "b")))
lh_models[[6]] <- brm(cv_effect ~ gen_time, data = traits_df,
                      family = "gaussian",
                      prior = c(set_prior("normal(0,5)", class = "b")))
lh_models[[7]] <- brm(cv_effect ~ seed_size, data = traits_df,
                      family = "gaussian",
                      prior = c(set_prior("normal(0,5)", class = "b")))
lh_models[[8]] <- brm(cv_effect ~ imperfect_trans, data = traits_df,
                      family = "gaussian",
                      prior = c(set_prior("normal(0,5)", class = "b")))



newdata_fit <- newdata %>%
  rename_with(.fn = str_c, pattern = ".x")  %>% 
  mutate(row_id = row_number()) 
  mutate(observed_max_age.fit = fitted(lh_models[[1]], newdata)[,1],
         max_age_99.fit = fitted(lh_models[[2]], newdata = newdata)[,1],
         R0.fit = fitted(lh_models[[3]], newdata = newdata)[,1],
         longev.fit = fitted(lh_models[[4]], newdata = newdata)[,1],
         mean_life_expect.fit = fitted(lh_models[[5]], newdata = newdata)[,1],
         gen_time.fit =fitted(lh_models[[6]], newdata = newdata)[,1],
         seed_size.fit =fitted(lh_models[[7]], newdata = newdata)[,1],
         imperfect_trans.fit =fitted(lh_models[[8]], newdata = newdata)[,1]) %>% 
  mutate(observed_max_age.lwr = fitted(lh_models[[1]], newdata = newdata, probs = c(0.05, 0.95))[,3],
         max_age_99.lwr = fitted(lh_models[[2]], newdata = newdata, probs = c(0.05, 0.95))[,3],
         R0.lwr = fitted(lh_models[[3]], newdata = newdata, probs = c(0.05, 0.95))[,3],
         longev.lwr = fitted(lh_models[[4]], newdata = newdata, probs = c(0.05, 0.95))[,3],
         mean_life_expect.lwr = fitted(lh_models[[5]], newdata = newdata, probs = c(0.05, 0.95))[,3],
         gen_time.lwr =fitted(lh_models[[6]], newdata = newdata, probs = c(0.05, 0.95))[,3],
         seed_size.lwr =fitted(lh_models[[7]], newdata = newdata, probs = c(0.05, 0.95))[,3],
         imperfect_trans.lwr =fitted(lh_models[[8]], newdata = newdata, probs = c(0.05, 0.95))[,3]) %>%
  mutate(observed_max_age.upr = fitted(lh_models[[1]], newdata = newdata, probs = c(0.05, 0.95))[,4],
         max_age_99.upr = fitted(lh_models[[2]], newdata = newdata, probs = c(0.05, 0.95))[,4],
         R0.upr = fitted(lh_models[[3]], newdata = newdata, probs = c(0.05, 0.95))[,4],
         longev.upr = fitted(lh_models[[4]], newdata = newdata, probs = c(0.05, 0.95))[,4],
         mean_life_expect.upr = fitted(lh_models[[5]], newdata = newdata, probs = c(0.05, 0.95))[,4],
         gen_time.upr =fitted(lh_models[[6]], newdata = newdata, probs = c(0.05, 0.95))[,4],
         seed_size.upr =fitted(lh_models[[7]], newdata = newdata, probs = c(0.05, 0.95))[,4],
         imperfect_trans.upr =fitted(lh_models[[8]], newdata = newdata, probs = c(0.05, 0.95))[,4]) %>%
  pivot_longer(cols = -row_id, names_to = c("name", "interval"), names_sep = "\\.") %>% 
  pivot_wider(id_cols = c(row_id,name), names_from = interval, values_from = value) %>% 
  mutate(name_label = case_when(name == "observed_max_age" ~ "Obs. Max Age", name == "max_age_99" ~ "99th Percentile Max Age", name == "R0" ~ "R0",
                                name == "longev" ~ "Longevity", name == "mean_life_expect" ~ "Mean Life Expectancy",
                                name == "gen_time" ~ "Generation Time", name == "seed_size" ~ "Seed Length (mm)", name == "imperfect_trans" ~ "Imperfect Transmission Rate")) 


oma_plot <- lh_plot <- ggplot()+
  geom_ribbon(data = filter(newdata_fit, name == "observed_max_age"), aes(ymin = lwr, ymax = upr, x = x), alpha = .4)+
  geom_line(data = filter(newdata_fit, name == "observed_max_age"), aes(y = fit, x = x))+
  geom_point(data = filter(traits_df_long, name == "observed_max_age"), aes(y = cv_effect, x = value))+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.x = element_text(size=10))+
  labs(x = "Obs. Max Age", y = expression(paste("Effect on CV", (lambda))), alpha = "p-value")
# oma_plot


ma99_plot <- lh_plot <- ggplot()+
  geom_ribbon(data = filter(newdata_fit, name == "max_age_99"), aes(ymin = lwr, ymax = upr, x = x), alpha = .4)+
  geom_line(data = filter(newdata_fit, name == "max_age_99"), aes(y = fit, x = x))+
  geom_point(data = filter(traits_df_long, name == "max_age_99"), aes(y = cv_effect, x = value))+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.x = element_text(size=10))+
  labs(x = "99th Percentile Max Age", y = "", alpha = "p-value")
# ma99_plot

R0_plot <- lh_plot <- ggplot()+
  geom_ribbon(data = filter(newdata_fit, name == "R0"), aes(ymin = lwr, ymax = upr, x = x), alpha = .4)+
  geom_line(data = filter(newdata_fit, name == "R0"), aes(y = fit, x = x))+
  geom_point(data = filter(traits_df_long, name == "R0"), aes(y = cv_effect, x = value))+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.x = element_text(size=10))+
  labs(x = "R0", y = "", alpha = "p-value")
# R0_plot

longev_plot <- lh_plot <- ggplot()+
  geom_ribbon(data = filter(newdata_fit, name == "longev"), aes(ymin = lwr, ymax = upr, x = x), alpha = .4)+
  geom_line(data = filter(newdata_fit, name == "longev"), aes(y = fit, x = x))+
  geom_point(data = filter(traits_df_long, name == "longev"), aes(y = cv_effect, x = value))+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.x = element_text(size=10))+
  labs(x = "Longevity", y = "", alpha = "p-value")
# longev_plot

mle_plot <- lh_plot <- ggplot()+
  geom_ribbon(data = filter(newdata_fit, name == "mean_life_expect"), aes(ymin = lwr, ymax = upr, x = x), alpha = .4)+
  geom_line(data = filter(newdata_fit, name == "mean_life_expect"), aes(y = fit, x = x))+
  geom_point(data = filter(traits_df_long, name == "mean_life_expect"), aes(y = cv_effect, x = value))+
  scale_x_continuous(breaks = c(3,6,9,12))+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.x = element_text(size=10))+
  labs(x = "Mean Life Expectancy", y = expression(paste("Effect on CV", (lambda))), alpha = "p-value")
# mle_plot

gt_plot <- lh_plot <- ggplot()+
  geom_ribbon(data = filter(newdata_fit, name == "gen_time"), aes(ymin = lwr, ymax = upr, x = x), alpha = .4)+
  geom_line(data = filter(newdata_fit, name == "gen_time"), aes(y = fit, x = x))+
  geom_point(data = filter(traits_df_long, name == "gen_time"), aes(y = cv_effect, x = value))+
  scale_alpha_manual(values = c(.1))+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.x = element_text(size=10))+
  labs(x = "Generation Time", y = "", alpha = "p-value")
# gt_plot

ss_plot <- lh_plot <- ggplot()+
  geom_ribbon(data = filter(newdata_fit, name == "seed_size"), aes(ymin = lwr, ymax = upr, x = x), alpha = .4)+
  geom_line(data = filter(newdata_fit, name == "seed_size"), aes(y = fit, x = x))+
  geom_point(data = filter(traits_df_long, name == "seed_size"), aes(y = cv_effect, x = value))+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.x = element_text(size=10))+
  labs(x = "Seed Length (mm)", y = "", alpha = "p-value")
# ss_plot

it_plot <- lh_plot <- ggplot()+
  geom_ribbon(data = filter(newdata_fit, name == "imperfect_trans"), aes(ymin = lwr, ymax = upr, x = x), alpha = .4)+
  geom_line(data = filter(newdata_fit, name == "imperfect_trans"), aes(y = fit, x = x))+
  geom_point(data = filter(traits_df_long, name == "imperfect_trans"), aes(y = cv_effect, x = value))+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.x = element_text(size=10))+
  labs(x = "Imperfect Transmission Rate", y = "", alpha = "p-value")
# it_plot

lh_plot_brms <- oma_plot + ma99_plot + R0_plot + longev_plot + mle_plot + gt_plot + ss_plot + it_plot + 
  plot_layout(nrow = 2) + plot_annotation(tag_levels = "A") 
ggsave(lh_plot_brms, filename = "lh_plot_brms.png", width = 10.5, height = 5)


#############################################################################################
####### Adjusting the variables for phylogenetic comparison ------------------
#############################################################################################
# reading in the phylo tree of vascular plants (Zanne et al. 2014: )
vasc.plant <- ape::read.tree("/Users/joshuacfowler/Downloads/Vascular_Plants_rooted.dated.tre")

# reading in phylo tree of epichloe (Leuchtmann et al. 2014; https://doi.org/10.3852/13-251)
epichloe <- ape::read.nexus("/Users/joshuacfowler/Documents/R_projects/Grass-Endophyte-Stochastic-Demography/Analyses/T68066.nex")


# adding the names of the tip labels to the traits df. In some cases have to take the label for closely related species in same genera
print(vasc.plant$tip.label)

traits_df$plant_label <- c("Agrostis_hyemalis", "Elymus_hystrix", "Elymus_virginicus", "Festuca_subverticillata", "Festuca_arundinacea", "Poa_alsodes", "Poa_sylvestris")

pruned.plant.tree <- drop.tip(vasc.plant, setdiff(vasc.plant$tip.label, traits_df$plant_label))
plot(pruned.plant.tree)
pruned.plant.tree$edge.length
# Names are from tree, which include the hosts sampled from Leuchtmann et al. 2014. This list takes the names for our host species, or those which are closest to our host based on literature.
#Need to find FESU completely. Need to double check that POAL and ELVI are close enough.
# Now taking just the closest endophyte species from the epichloe tree. This is a bit more complicated in reality, but basing this off of XXX sources for XXX species
print(epichloe$tip.label)
traits_df$epichloe_label <- c("Epichloe_amarillans_ATCC_200744_ex._Agrostis_hyemalis",
                        "Epichloe_elymi_ATCC_201555_ex._Elymus_villosus", 
                        "Epichloe_amarillans_E1087_ex._Elymus_virginicus", 
                        "Epichloe_sp._FalTG1_e507_tubBbo_ex._Festuca_altissima",
                        "Epichloe_coenophiala_ATCC_90664_tubBty_ex._Schedonorus_arundinaceus_6x",
                        "Epichloe_sp._PauTG1_e55_tubBty_ex._Poa_autumnalis", 
                        "Epichloe_typhina_subsp._poae_e1097_ex._Poa_sylvestris") 

pruned.epichloe.tree<-drop.tip(epichloe, setdiff(epichloe$tip.label, traits_df$epichloe_label))
plot(pruned.epichloe.tree)

pruned.epichloe.tree

# plot(epichloe)

# Calculating phylogenetically independent contrasts for each life history trait and the effect on CV

plant_pic_values <- as_tibble(apply(column_to_rownames(traits_df, var = "plant_label")[3:16], MARGIN = 2, FUN = pic, phy = pruned.plant.tree))
epichloe_pic_values <- as_tibble(apply(column_to_rownames(traits_df, var = "epichloe_label")[3:15], MARGIN = 2, FUN = pic, phy = pruned.epichloe.tree))

#############################################################################################
####### Fitting regressions for the phylogenetically corrected traits and cv effects ------------------
#############################################################################################

summary(lm(cv_effect ~ observed_max_age - 1, data = plant_pic_values))
summary(lm(cv_effect ~ max_age_99 - 1, data = plant_pic_values))
summary(lm(cv_effect ~ R0 - 1, data = plant_pic_values))
summary(lm(cv_effect ~ longev - 1, data = plant_pic_values))
summary(lm(cv_effect ~ mean_life_expect - 1, data = plant_pic_values))
summary(lm(cv_effect ~ gen_time - 1, data = plant_pic_values))
summary(lm(cv_effect ~ seed_size - 1, data = plant_pic_values))
summary(lm(cv_effect ~ imperfect_trans - 1, data = plant_pic_values))

summary(lm(cv_effect ~ observed_max_age - 1, data = epichloe_pic_values))
summary(lm(cv_effect ~ max_age_99 - 1, data = epichloe_pic_values))
summary(lm(cv_effect ~ R0 - 1, data = epichloe_pic_values))
summary(lm(cv_effect ~ longev - 1, data = epichloe_pic_values))
summary(lm(cv_effect ~ mean_life_expect - 1, data = epichloe_pic_values))
summary(lm(cv_effect ~ gen_time - 1, data = epichloe_pic_values))
summary(lm(cv_effect ~ seed_size - 1, data = epichloe_pic_values))
summary(lm(cv_effect ~ imperfect_trans - 1, data = epichloe_pic_values))

##### Fitting regressions corrected for host plant phylogeny ####
plant_cv_models <- list()
plant_cv_models[[1]] <- brm(cv_effect ~ observed_max_age -1, data = plant_pic_values,
                         family = "gaussian",
                         prior = c(set_prior("normal(0,5)", class = "b")))
plant_cv_models[[2]] <- brm(cv_effect ~ max_age_99 -1, data = plant_pic_values,
                         family = "gaussian",
                         prior = c(set_prior("normal(0,5)", class = "b")))
plant_cv_models[[3]] <- brm(cv_effect ~ pic_R0 - 1, data = plant_pic_values,
                         family = "gaussian",
                         prior = c(set_prior("normal(0,5)", class = "b")))
plant_cv_models[[4]] <- brm(pic_cv_effect ~ pic_longev - 1, data = plant_pic_values,
                         family = "gaussian",
                         prior = c(set_prior("normal(0,5)", class = "b")))
plant_cv_models[[5]] <- brm(pic_cv_effect ~ pic_mle - 1, data = plant_pic_values,
                         family = "gaussian",
                         prior = c(set_prior("normal(0,5)", class = "b")))
plant_cv_models[[6]] <- brm(pic_cv_effect ~ pic_gen_time - 1, data = plant_pic_values,
                         family = "gaussian",
                         prior = c(set_prior("normal(0,5)", class = "b")))
plant_cv_models[[7]] <- brm(pic_cv_effect ~ pic_seed_size - 1, data = plant_pic_values,
                         family = "gaussian",
                         prior = c(set_prior("normal(0,5)", class = "b")))
# the imperf. transmission model has a few divergent transitions...
plant_cv_models[[8]] <- brm(pic_cv_effect ~ pic_imperfect_trans - 1, data = plant_pic_values,
                         family = "gaussian",
                         prior = c(set_prior("normal(0,1)", class = "b")))

# making data ranges for predictions
newdata_plant <- data.frame(pic_oma = seq(from = min(plant_pic_values$pic_oma), to = max(plant_pic_values$pic_oma), length.out = 100),
                      pic_ma99 = seq(from = min(plant_pic_values$pic_ma99), to = max(plant_pic_values$pic_ma99), length.out = 100),
                      pic_gen_time = seq(from = min(plant_pic_values$pic_gen_time), to = max(plant_pic_values$pic_gen_time), length.out = 100),
                      pic_longev = seq(from = min(plant_pic_values$pic_longev), to = max(plant_pic_values$pic_longev), length.out = 100),
                      pic_mle = seq(from = min(plant_pic_values$pic_mle), to = max(plant_pic_values$pic_mle), length.out = 100),
                      pic_R0 = seq(from = min(plant_pic_values$pic_R0), to = max(plant_pic_values$pic_R0), length.out = 100),
                      pic_seed_size = seq(from = min(plant_pic_values$pic_seed_size), to = max(plant_pic_values$pic_seed_size), length.out = 100),
                      pic_imperfect_trans = seq(from = min(plant_pic_values$pic_imperfect_trans), to = max(plant_pic_values$pic_imperfect_trans), length.out = 100))
newdata_plant_fit <- newdata_plant %>%
  rename_with(.fn = str_c, pattern = ".x")  %>% 
  mutate(row_id = row_number()) %>% 
mutate(pic_oma.fit = fitted(plant_cv_models[[1]], newdata = newdata_plant)[,1],
       pic_ma99.fit = fitted(plant_cv_models[[2]], newdata = newdata_plant)[,1],
       pic_R0.fit = fitted(plant_cv_models[[3]], newdata = newdata_plant)[,1],
       pic_longev.fit = fitted(plant_cv_models[[4]], newdata = newdata_plant)[,1],
       pic_mle.fit = fitted(plant_cv_models[[5]], newdata = newdata_plant)[,1],
       pic_gen_time.fit =fitted(plant_cv_models[[6]], newdata = newdata_plant)[,1],
       pic_seed_size.fit =fitted(plant_cv_models[[7]], newdata = newdata_plant)[,1],
       pic_imperfect_trans.fit =fitted(plant_cv_models[[8]], newdata = newdata_plant)[,1]) %>% 
  mutate(pic_oma.lwr = fitted(plant_cv_models[[1]], newdata = newdata_plant, probs = c(0.05, 0.95))[,3],
         pic_ma99.lwr = fitted(plant_cv_models[[2]], newdata = newdata_plant, probs = c(0.05, 0.95))[,3],
         pic_R0.lwr = fitted(plant_cv_models[[3]], newdata = newdata_plant, probs = c(0.05, 0.95))[,3],
         pic_longev.lwr = fitted(plant_cv_models[[4]], newdata = newdata_plant, probs = c(0.05, 0.95))[,3],
         pic_mle.lwr = fitted(plant_cv_models[[5]], newdata = newdata_plant, probs = c(0.05, 0.95))[,3],
         pic_gen_time.lwr =fitted(plant_cv_models[[6]], newdata = newdata_plant, probs = c(0.05, 0.95))[,3],
         pic_seed_size.lwr =fitted(plant_cv_models[[7]], newdata = newdata_plant, probs = c(0.05, 0.95))[,3],
         pic_imperfect_trans.lwr =fitted(plant_cv_models[[8]], newdata = newdata_plant, probs = c(0.05, 0.95))[,3]) %>%
  mutate(pic_oma.upr = fitted(plant_cv_models[[1]], newdata = newdata_plant, probs = c(0.05, 0.95))[,4],
         pic_ma99.upr = fitted(plant_cv_models[[2]], newdata = newdata_plant, probs = c(0.05, 0.95))[,4],
         pic_R0.upr = fitted(plant_cv_models[[3]], newdata = newdata_plant, probs = c(0.05, 0.95))[,4],
         pic_longev.upr = fitted(plant_cv_models[[4]], newdata = newdata_plant, probs = c(0.05, 0.95))[,4],
         pic_mle.upr = fitted(plant_cv_models[[5]], newdata = newdata_plant, probs = c(0.05, 0.95))[,4],
         pic_gen_time.upr =fitted(plant_cv_models[[6]], newdata = newdata_plant, probs = c(0.05, 0.95))[,4],
         pic_seed_size.upr =fitted(plant_cv_models[[7]], newdata = newdata_plant, probs = c(0.05, 0.95))[,4],
         pic_imperfect_trans.upr =fitted(plant_cv_models[[8]], newdata = newdata_plant, probs = c(0.05, 0.95))[,4]) %>% 
  pivot_longer(cols = -row_id, names_to = c("name", "interval"), names_sep = "\\.") %>% 
  pivot_wider(id_cols = c(row_id,name), names_from = interval, values_from = value)  


# plotting each relationship

oma_plant_plot <- ggplot()+
  geom_ribbon(data = filter(newdata_plant_fit, name == "pic_oma"), aes(ymin = lwr, ymax = upr, x = x), alpha = .4)+
  geom_line(data = filter(newdata_plant_fit, name == "pic_oma"), aes(y = fit, x = x))+
  geom_point(data = plant_pic_values, aes(y = pic_cv_effect, x = pic_oma))+
  theme_classic()+
  theme(axis.title.x = element_text(size=10))+
  labs(x = "PIC of Obs. Max Age", y = expression(paste("PIC of Effect on CV", (lambda))))
oma_plant_plot


ma99_plant_plot <- ggplot()+
  geom_ribbon(data = filter(newdata_plant_fit, name == "pic_ma99"), aes(ymin = lwr, ymax = upr, x = x), alpha = .4)+
  geom_line(data = filter(newdata_plant_fit, name == "pic_ma99"), aes(y = fit, x = x))+
  geom_point(data = plant_pic_values, aes(y = pic_cv_effect, x = pic_ma99))+
  theme_classic()+
  theme(axis.title.x = element_text(size=10))+
  labs(x = "PIC of 99th Percentile Max Age", y = "")
ma99_plant_plot

R0_plant_plot <- ggplot()+
  geom_ribbon(data = filter(newdata_plant_fit, name == "pic_R0"), aes(ymin = lwr, ymax = upr, x = x), alpha = .4)+
  geom_line(data = filter(newdata_plant_fit, name == "pic_R0"), aes(y = fit, x = x))+
  geom_point(data = plant_pic_values, aes(y = pic_cv_effect, x = pic_R0))+
  theme_classic()+
  theme(axis.title.x = element_text(size=10))+
  labs(x = "PIC of R0", y = "")
R0_plant_plot

longev_plant_plot <- ggplot()+
  geom_ribbon(data = filter(newdata_plant_fit, name == "pic_longev"), aes(ymin = lwr, ymax = upr, x = x), alpha = .4)+
  geom_line(data = filter(newdata_plant_fit, name == "pic_longev"), aes(y = fit, x = x))+
  geom_point(data = plant_pic_values, aes(y = pic_cv_effect, x = pic_longev))+
  theme_classic()+
  theme(axis.title.x = element_text(size=10))+
  labs(x = "PIC of Longevity", y = "")
longev_plant_plot

mle_plant_plot <-  ggplot()+
  geom_ribbon(data = filter(newdata_plant_fit, name == "pic_mle"), aes(ymin = lwr, ymax = upr, x = x), alpha = .4)+
  geom_line(data = filter(newdata_plant_fit, name == "pic_mle"), aes(y = fit, x = x))+
  geom_point(data = plant_pic_values, aes(y = pic_cv_effect, x = pic_mle))+
  theme_classic()+
  theme(axis.title.x = element_text(size=10))+
  labs(x = "PIC of Mean Life Expectancy", y = expression(paste("PIC of Effect on CV", (lambda))))
mle_plant_plot

gt_plant_plot <-  ggplot()+
  geom_ribbon(data = filter(newdata_plant_fit, name == "pic_gen_time"), aes(ymin = lwr, ymax = upr, x = x), alpha = .4)+
  geom_line(data = filter(newdata_plant_fit, name == "pic_gen_time"), aes(y = fit, x = x))+
  geom_point(data = plant_pic_values, aes(y = pic_cv_effect, x = pic_gen_time))+
  theme_classic()+
  theme(axis.title.x = element_text(size=10))+
  labs(x = "PIC of Generation Time", y = "")
gt_plant_plot

ss_plant_plot <- ggplot()+
  geom_ribbon(data = filter(newdata_plant_fit, name == "pic_seed_size"), aes(ymin = lwr, ymax = upr, x = x), alpha = .4)+
  geom_line(data = filter(newdata_plant_fit, name == "pic_seed_size"), aes(y = fit, x = x))+
  geom_point(data = plant_pic_values, aes(y = pic_cv_effect, x = pic_seed_size))+
  theme_classic()+
  theme(axis.title.x = element_text(size=10))+
  labs(x = "PIC of Seed Length (mm)", y = "")
ss_plant_plot

it_plant_plot <- ggplot()+
  geom_ribbon(data = filter(newdata_plant_fit, name == "pic_imperfect_trans"), aes(ymin = lwr, ymax = upr, x = x), alpha = .4)+
  geom_line(data = filter(newdata_plant_fit, name == "pic_imperfect_trans"), aes(y = fit, x = x))+
  geom_point(data = plant_pic_values, aes(y = pic_cv_effect, x = pic_imperfect_trans))+
  theme_classic()+
  theme(axis.title.x = element_text(size=10))+
  labs(x = "PIC of Imperfect Transmission Rate", y = "")
it_plant_plot

lh_plant_plot <- oma_plant_plot + ma99_plant_plot + R0_plant_plot + longev_plant_plot + mle_plant_plot + gt_plant_plot + ss_plant_plot + it_plant_plot + 
  plot_layout(nrow = 2) + plot_annotation(tag_levels = "A")
lh_plant_plot
ggsave(lh_plant_plot, filename = "lh_plant_plot.png", width = 10.5, height = 5)


###### Fitting regressions corrected for symbiont phylogeny #####
epichloe_cv_models <- list()
epichloe_cv_models[[1]] <- brm(pic_cv_effect ~ pic_oma -1, data = epichloe_pic_values,
                         family = "gaussian",
                         prior = c(set_prior("normal(0,5)", class = "b")))
epichloe_cv_models[[2]] <- brm(pic_cv_effect ~ pic_ma99 -1, data = epichloe_pic_values,
                         family = "gaussian",
                         prior = c(set_prior("normal(0,5)", class = "b")))
epichloe_cv_models[[3]] <- brm(pic_cv_effect ~ pic_R0 - 1, data = epichloe_pic_values,
                         family = "gaussian",
                         prior = c(set_prior("normal(0,5)", class = "b")))
# model 4 has a few div transitions
epichloe_cv_models[[4]] <- brm(pic_cv_effect ~ pic_longev - 1, data = epichloe_pic_values,
                         family = "gaussian",
                         prior = c(set_prior("normal(0,5)", class = "b")))
epichloe_cv_models[[5]] <- brm(pic_cv_effect ~ pic_mle - 1, data = epichloe_pic_values,
                         family = "gaussian",
                         prior = c(set_prior("normal(0,5)", class = "b")))
epichloe_cv_models[[6]] <- brm(pic_cv_effect ~ pic_gen_time - 1, data = epichloe_pic_values,
                         family = "gaussian",
                         prior = c(set_prior("normal(0,5)", class = "b")))
epichloe_cv_models[[7]] <- brm(pic_cv_effect ~ pic_seed_size - 1, data = epichloe_pic_values,
                         family = "gaussian",
                         prior = c(set_prior("normal(0,5)", class = "b")))
epichloe_cv_models[[8]] <- brm(pic_cv_effect ~ pic_imperfect_trans - 1, data = epichloe_pic_values,
                         family = "gaussian",
                         prior = c(set_prior("normal(0,5)", class = "b")))

# making data ranges for predictions
newdata_epichloe <- data.frame(pic_oma = seq(from = min(epichloe_pic_values$pic_oma), to = max(epichloe_pic_values$pic_oma), length.out = 100),
                               pic_ma99 = seq(from = min(epichloe_pic_values$pic_ma99), to = max(epichloe_pic_values$pic_ma99), length.out = 100),
                               pic_gen_time = seq(from = min(epichloe_pic_values$pic_gen_time), to = max(epichloe_pic_values$pic_gen_time), length.out = 100),
                               pic_longev = seq(from = min(epichloe_pic_values$pic_longev), to = max(epichloe_pic_values$pic_longev), length.out = 100),
                               pic_mle = seq(from = min(epichloe_pic_values$pic_mle), to = max(epichloe_pic_values$pic_mle), length.out = 100),
                               pic_R0 = seq(from = min(epichloe_pic_values$pic_R0), to = max(epichloe_pic_values$pic_R0), length.out = 100),
                               pic_seed_size = seq(from = min(epichloe_pic_values$pic_seed_size), to = max(epichloe_pic_values$pic_seed_size), length.out = 100),
                               pic_imperfect_trans = seq(from = min(epichloe_pic_values$pic_imperfect_trans), to = max(epichloe_pic_values$pic_imperfect_trans), length.out = 100))
newdata_epichloe_fit <- newdata_epichloe %>%
  rename_with(.fn = str_c, pattern = ".x")  %>% 
  mutate(row_id = row_number()) %>% 
  mutate(pic_oma.fit = fitted(epichloe_cv_models[[1]], newdata = newdata_epichloe)[,1],
         pic_ma99.fit = fitted(epichloe_cv_models[[2]], newdata = newdata_epichloe)[,1],
         pic_R0.fit = fitted(epichloe_cv_models[[3]], newdata = newdata_epichloe)[,1],
         pic_longev.fit = fitted(epichloe_cv_models[[4]], newdata = newdata_epichloe)[,1],
         pic_mle.fit = fitted(epichloe_cv_models[[5]], newdata = newdata_epichloe)[,1],
         pic_gen_time.fit =fitted(epichloe_cv_models[[6]], newdata = newdata_epichloe)[,1],
         pic_seed_size.fit =fitted(epichloe_cv_models[[7]], newdata = newdata_epichloe)[,1],
         pic_imperfect_trans.fit =fitted(epichloe_cv_models[[8]], newdata = newdata_epichloe)[,1]) %>% 
  mutate(pic_oma.lwr = fitted(epichloe_cv_models[[1]], newdata = newdata_epichloe, probs = c(0.05, 0.95))[,3],
         pic_ma99.lwr = fitted(epichloe_cv_models[[2]], newdata = newdata_epichloe, probs = c(0.05, 0.95))[,3],
         pic_R0.lwr = fitted(epichloe_cv_models[[3]], newdata = newdata_epichloe, probs = c(0.05, 0.95))[,3],
         pic_longev.lwr = fitted(epichloe_cv_models[[4]], newdata = newdata_epichloe, probs = c(0.05, 0.95))[,3],
         pic_mle.lwr = fitted(epichloe_cv_models[[5]], newdata = newdata_epichloe, probs = c(0.05, 0.95))[,3],
         pic_gen_time.lwr =fitted(epichloe_cv_models[[6]], newdata = newdata_epichloe, probs = c(0.05, 0.95))[,3],
         pic_seed_size.lwr =fitted(epichloe_cv_models[[7]], newdata = newdata_epichloe, probs = c(0.05, 0.95))[,3],
         pic_imperfect_trans.lwr =fitted(epichloe_cv_models[[8]], newdata = newdata_epichloe, probs = c(0.05, 0.95))[,3]) %>%
  mutate(pic_oma.upr = fitted(epichloe_cv_models[[1]], newdata = newdata_epichloe, probs = c(0.05, 0.95))[,4],
         pic_ma99.upr = fitted(epichloe_cv_models[[2]], newdata = newdata_epichloe, probs = c(0.05, 0.95))[,4],
         pic_R0.upr = fitted(epichloe_cv_models[[3]], newdata = newdata_epichloe, probs = c(0.05, 0.95))[,4],
         pic_longev.upr = fitted(epichloe_cv_models[[4]], newdata = newdata_epichloe, probs = c(0.05, 0.95))[,4],
         pic_mle.upr = fitted(epichloe_cv_models[[5]], newdata = newdata_epichloe, probs = c(0.05, 0.95))[,4],
         pic_gen_time.upr =fitted(epichloe_cv_models[[6]], newdata = newdata_epichloe, probs = c(0.05, 0.95))[,4],
         pic_seed_size.upr =fitted(epichloe_cv_models[[7]], newdata = newdata_epichloe, probs = c(0.05, 0.95))[,4],
         pic_imperfect_trans.upr =fitted(epichloe_cv_models[[8]], newdata = newdata_epichloe, probs = c(0.05, 0.95))[,4]) %>%
  pivot_longer(cols = -row_id, names_to = c("name", "interval"), names_sep = "\\.") %>% 
  pivot_wider(id_cols = c(row_id,name), names_from = interval, values_from = value) 

# plotting each relationship
oma_epichloe_plot <- ggplot()+
  geom_ribbon(data = filter(newdata_epichloe_fit, name == "pic_oma"), aes(ymin = lwr, ymax = upr, x = x), alpha = .4)+
  geom_line(data = filter(newdata_epichloe_fit, name == "pic_oma"), aes(y = fit, x = x))+
  geom_point(data = epichloe_pic_values, aes(y = pic_cv_effect, x = pic_oma))+
  theme_classic()+
  theme(axis.title.x = element_text(size=10))+
  labs(x = "PIC of Obs. Max Age", y = expression(paste("PIC of Effect on CV", (lambda))))
oma_epichloe_plot


ma99_epichloe_plot <- ggplot()+
  geom_ribbon(data = filter(newdata_epichloe_fit, name == "pic_ma99"), aes(ymin = lwr, ymax = upr, x = x), alpha = .4)+
  geom_line(data = filter(newdata_epichloe_fit, name == "pic_ma99"), aes(y = fit, x = x))+
  geom_point(data = epichloe_pic_values, aes(y = pic_cv_effect, x = pic_ma99))+
  theme_classic()+
  theme(axis.title.x = element_text(size=10))+
  labs(x = "PIC of 99th Percentile Max Age", y = "")
ma99_epichloe_plot

R0_epichloe_plot <- ggplot()+
  geom_ribbon(data = filter(newdata_epichloe_fit, name == "pic_R0"), aes(ymin = lwr, ymax = upr, x = x), alpha = .4)+
  geom_line(data = filter(newdata_epichloe_fit, name == "pic_R0"), aes(y = fit, x = x))+
  geom_point(data = epichloe_pic_values, aes(y = pic_cv_effect, x = pic_R0))+
  theme_classic()+
  theme(axis.title.x = element_text(size=10))+
  labs(x = "PIC of R0", y = "",)
R0_epichloe_plot

longev_epichloe_plot <- ggplot()+
  geom_ribbon(data = filter(newdata_epichloe_fit, name == "pic_longev"), aes(ymin = lwr, ymax = upr, x = x), alpha = .4)+
  geom_line(data = filter(newdata_epichloe_fit, name == "pic_longev"), aes(y = fit, x = x))+
  geom_point(data = epichloe_pic_values, aes(y = pic_cv_effect, x = pic_longev))+
  theme_classic()+
  theme(axis.title.x = element_text(size=10))+
  labs(x = "PIC of Longevity", y = "")
longev_epichloe_plot

mle_epichloe_plot <-  ggplot()+
  geom_ribbon(data = filter(newdata_epichloe_fit, name == "pic_mle"), aes(ymin = lwr, ymax = upr, x = x), alpha = .4)+
  geom_line(data = filter(newdata_epichloe_fit, name == "pic_mle"), aes(y = fit, x = x))+
  geom_point(data = epichloe_pic_values, aes(y = pic_cv_effect, x = pic_mle))+
  theme_classic()+
  theme(axis.title.x = element_text(size=10))+
  labs(x = "PIC of Mean Life Expectancy", y = expression(paste("PIC of Effect on CV", (lambda))))
mle_epichloe_plot

gt_epichloe_plot <-  ggplot()+
  geom_ribbon(data = filter(newdata_epichloe_fit, name == "pic_gen_time"), aes(ymin = lwr, ymax = upr, x = x), alpha = .4)+
  geom_line(data = filter(newdata_epichloe_fit, name == "pic_gen_time"), aes(y = fit, x = x))+
  geom_point(data = epichloe_pic_values, aes(y = pic_cv_effect, x = pic_gen_time))+
  theme_classic()+
  theme(axis.title.x = element_text(size=10))+
  labs(x = "PIC of Generation Time", y = "")
gt_epichloe_plot

ss_epichloe_plot <- ggplot()+
  geom_ribbon(data = filter(newdata_epichloe_fit, name == "pic_seed_size"), aes(ymin = lwr, ymax = upr, x = x), alpha = .4)+
  geom_line(data = filter(newdata_epichloe_fit, name == "pic_seed_size"), aes(y = fit, x = x))+
  geom_point(data = epichloe_pic_values, aes(y = pic_cv_effect, x = pic_seed_size))+
  theme_classic()+
  theme(axis.title.x = element_text(size=10))+
  labs(x = "PIC of Seed Length (mm)", y = "")
ss_epichloe_plot

it_epichloe_plot <- ggplot()+
  geom_ribbon(data = filter(newdata_epichloe_fit, name == "pic_imperfect_trans"), aes(ymin = lwr, ymax = upr, x = x), alpha = .4)+
  geom_line(data = filter(newdata_epichloe_fit, name == "pic_imperfect_trans"), aes(y = fit, x = x))+
  geom_point(data = epichloe_pic_values, aes(y = pic_cv_effect, x = pic_imperfect_trans))+
  theme_classic()+
  theme(axis.title.x = element_text(size=10))+
  labs(x = "PIC of Imperfect Transmission Rate", y = "")
it_epichloe_plot

lh_epichloe_plot <- oma_epichloe_plot + ma99_epichloe_plot + R0_epichloe_plot + longev_epichloe_plot + mle_epichloe_plot + gt_epichloe_plot + ss_epichloe_plot + it_epichloe_plot + 
  plot_layout(nrow = 2) + plot_annotation(tag_levels = "A")
lh_epichloe_plot
ggsave(lh_epichloe_plot, filename = "lh_epichloe_plot_notscaled.png", width = 10.5, height = 5)




cor.test(plant_civ_values$cv_effect, plant_civ_values$observed_max_age)
cor.test(plant_civ_values$cv_effect, plant_civ_values$max_age_99)
cor.test(plant_civ_values$cv_effect, plant_civ_values$imperfect_trans)


pic.cv <- pic(traits_df$cv_effect, pruned.plant.tree, scaled = TRUE)
pic.oma <- pic(traits_df$observed_max_age, pruned.plant.tree, scaled = FALSE)
plot( pic.oma, pic.cv)
plot(traits_df$observed_max_age, traits_df$cv_effect)
cor.test(pic.oma, pic.cv)
pic_df <- data.frame(pic.cv, pic.oma)

lh_models[[2]] <- brm(pic.cv ~ pic.oma, data = pic_df,
                      family = "gaussian",
                      prior = c(set_prior("normal(0,5)", class = "b")))
pic_model <- lm(pic_cv ~ pic.oma - 1, data = pic_df)


# creating a co-variance matrix of the plant phylogeny
A <- ape::vcv.phylo(pruned.plant.tree)

# trying out a phylogenetic mixed-effects model
## MCMC settings
mcmc_pars <- list(
  iter = 10000, 
  warmup = 7500, 
  thin = 1, 
  chains = 3
)
           
           
           
meanof_cv <- mean(traits_df$cv_effect)
sdof_cv <- sd(traits_df$cv_effect)
plant_models <- list()
plant_models[[1]] <- brm(cv_effect|mi(sdof_cv_effect) ~ observed_max_age + (1|gr(plant_label, cov = A)),
                         data = traits_df,
                         family = gaussian(),
                         data2 = list(A = A),
                         prior = c(
                           prior(normal(0, .1), "b"),
                           prior(normal(-0.0337147, .5), "Intercept"),
                           prior(normal(.0,.1), class = "sd", lb = 0),
                           prior(normal(.04417872,.01), class = "sigma", lb = 0)),
                         control = list(adapt_delta = 0.999,
                                        max_treedepth = 20),
                         iter = mcmc_pars$iter,
                         warmup = mcmc_pars$warmup,
                         save_pars = save_pars(latent = TRUE))
print(plant_models[[1]])
plot(conditional_effects(plant_models[[1]], method = "fitted"))
posterior_pm <- as.array(plant_models[[1]])
posterior_df <- as_tibble(cubelyr::as.tbl_cube(posterior_pm))
ggplot(data = filter(posterior_df, variable == "b_observed_max_age"))+
  geom_histogram(aes(x = posterior_pm), bins = 100)+
  facet_wrap(~variable, scales = "free_x")



plant_models[[2]] <- brm(cv_effect|mi(sdof_cv_effect) ~ max_age_99 + (1|gr(plant_label, cov = A)),
                         data = traits_df,
                         family = gaussian(),
                         data2 = list(A = A),
                         prior = c(
                           prior(normal(0, .1), "b"),
                           prior(normal(-0.0337147, .5), "Intercept"),
                           prior(normal(0,.1), class = "sd", lb = 0),
                           prior(normal(.04417872,.01), class = "sigma", lb = 0)),
                         control = list(adapt_delta = 0.999,
                                        max_treedepth = 20),
                         iter = mcmc_pars$iter,
                         warmup = mcmc_pars$warmup,
                         save_pars = save_pars(latent = TRUE))
plant_models[[3]] <- brm(cv_effect|mi(sdof_cv_effect) ~ R0 + (1|gr(plant_label, cov = A)),
                         data = traits_df,
                         family = gaussian(),
                         data2 = list(A = A),
                         prior = c(
                           prior(normal(0, .1), "b"),
                           prior(normal(-0.0337147, .5), "Intercept"),
                           prior(normal(0,.1), class = "sd", lb = 0),
                           prior(normal(.04417872,.01), class = "sigma", lb = 0)),
                         control = list(adapt_delta = 0.999,
                                        max_treedepth = 20),
                         iter = mcmc_pars$iter,
                         warmup = mcmc_pars$warmup,
                         save_pars = save_pars(latent = TRUE))
plant_models[[4]] <- brm(cv_effect|mi(sdof_cv_effect) ~ longev + (1|gr(plant_label, cov = A)),
                         data = traits_df,
                         family = gaussian(),
                         data2 = list(A = A),
                         prior = c(
                           prior(normal(0, .1), "b"),
                           prior(normal(-0.0337147, .5), "Intercept"),
                           prior(normal(0,.1), class = "sd", lb = 0),
                           prior(normal(.04417872,.01), class = "sigma", lb = 0)),
                         control = list(adapt_delta = 0.999,
                                        max_treedepth = 20),
                         iter = mcmc_pars$iter,
                         warmup = mcmc_pars$warmup,
                         save_pars = save_pars(latent = TRUE))
plant_models[[5]] <- brm(cv_effect|mi(sdof_cv_effect) ~ mean_life_expect + (1|gr(plant_label, cov = A)),
                         data = traits_df,
                         family = gaussian(),
                         data2 = list(A = A),
                         prior = c(
                           prior(normal(0, .1), "b"),
                           prior(normal(-0.0337147, .5), "Intercept"),
                           prior(normal(0,.1), class = "sd", lb = 0),
                           prior(normal(.04417872,.01), class = "sigma", lb = 0)),
                         control = list(adapt_delta = 0.999,
                                        max_treedepth = 20),
                         iter = mcmc_pars$iter,
                         warmup = mcmc_pars$warmup,
                         save_pars = save_pars(latent = TRUE))
plant_models[[6]] <- brm(cv_effect|mi(sdof_cv_effect) ~ gen_time + (1|gr(plant_label, cov = A)),
                       data = traits_df,
                       family = gaussian(),
                       data2 = list(A = A),
                       prior = c(
                         prior(normal(0, .05), "b"),
                         prior(normal(-0.0337147, .5), "Intercept"),
                         prior(normal(0,.1), class = "sd", lb = 0),
                         prior(normal(.04417872,.01), class = "sigma", lb = 0)),
                       control = list(adapt_delta = 0.999,
                                      max_treedepth = 20),
                       iter = mcmc_pars$iter,
                       warmup = mcmc_pars$warmup,
                       save_pars = save_pars(latent = TRUE))
plant_models[[7]] <- brm(cv_effect|mi(sdof_cv_effect) ~ seed_size + (1|gr(plant_label, cov = A)),
                         data = traits_df,
                         family = gaussian(),
                         data2 = list(A = A),
                         prior = c(
                           prior(normal(0, .01), "b"),
                           prior(normal(-0.0337147, .5), "Intercept"),
                           prior(normal(0,.1), class = "sd", lb = 0),
                           prior(normal(.04417872,.01), class = "sigma", lb = 0)),
                         control = list(adapt_delta = 0.999,
                                        max_treedepth = 20),
                         iter = mcmc_pars$iter,
                         warmup = mcmc_pars$warmup,
                         save_pars = save_pars(latent = TRUE))
plant_models[[8]] <- brm(cv_effect|mi(sdof_cv_effect) ~ imperfect_trans + (1|gr(plant_label, cov = A)),
                         data = traits_df,
                         family = gaussian(),
                         data2 = list(A = A),
                         prior = c(
                           prior(normal(0, .1), "b"),
                           prior(normal(-0.0337147, .5), "Intercept"),
                           prior(normal(0,.1), class = "sd", lb = 0),
                           prior(normal(.04417872,.01), class = "sigma", lb = 0)),
                         control = list(adapt_delta = 0.999,
                                        max_treedepth = 20),
                         iter = mcmc_pars$iter,
                         warmup = mcmc_pars$warmup,
                         save_pars = save_pars(latent = TRUE))

posterior_pm <- as.array(plant_models[[7]])
# nuts_pm <- nuts_params(plant_models[[1]])
# mcmc_trace(plant_model, pars = "gen_time")
# mcmc_parcoord(posterior_pm, np = nuts_pm, pars = c("bsp_megen_timesd_gen_time", "b_Intercept", "sd_plant_label__Intercept", "sigma", "meanme_megen_time", "sdme_megen_time"))

plot(conditional_effects(plant_models[[7]], method = "fitted"))

posterior_df <- as_tibble(cubelyr::as.tbl_cube(posterior_pm))
ggplot(data = filter(posterior_df, variable == "b_seed_size"))+
  geom_histogram(aes(x = posterior_pm), bins = 100)+
  facet_wrap(~variable, scales = "free_x")





# getting pagel's lambda 
hyp <- "sd_plant_label__Intercept^2 / (sd_plant_label__Intercept^2 + sigma^2) = 0"
(hyp1 <- hypothesis(plant_models[[1]], hyp, class = NULL))
(hyp1 <- hypothesis(plant_models[[2]], hyp, class = NULL))
(hyp1 <- hypothesis(plant_models[[3]], hyp, class = NULL))
(hyp1 <- hypothesis(plant_models[[4]], hyp, class = NULL))


# this is not right, but trying to mess around with this
plant_pic_values$sd_gen_time <- gen_time_summary[1:6,4]
plant_pic_values$sdof_cv_effect <- lambda_cv_diff[1:6,8]

hist(plant_pic_values$gen_time)

plant_model <- brm(
  cv_effect|mi(sdof_cv_effect) ~  me(gen_time, sd_gen_time) -1,
  data = plant_pic_values,
  family = gaussian(),
  prior = c(
    prior(normal(0, 1), "b"),
    # prior(normal(0, .5), "Intercept"),
    prior(normal(0, 1), class = "meanme"),
    prior(student_t(3,0,2.5), class = "sdme"),
    prior(exponential(1), "sigma")),
  control = list(adapt_delta = 0.999),
  iter = mcmc_pars$iter,
  warmup = mcmc_pars$warmup,
  save_pars = save_pars(latent = TRUE)
)

install.packages("performance")
r2_bayes(plant_model)
pp_check(plant_model)
plot(conditional_effects(plant_model, method = "fitted"))
mcmc_hist(plant_model, pars = "bsp_megen_timesd_gen_time")
post_samps <- extract_draws(plant_model)
plant_model$model
####

epichloe_distances <- cophenetic.phylo(epichloe)




# Calculating pairwise distances between endophyte effects
endo <- tibble(lambda_var_diff[1:7,1], lambda_mean_diff[1:7,1])
dist_e <- dist(endo, method = "euclidean", diag = TRUE)

#pairwise distance between gen time
gen_time <- tibble(gen_time_summary[,3], longev_summary[,3], mean_life_expect_summary[,3], R0_summary[,3], seed_size, relatedness)
dist_G <- dist(gen_time, method = "euclidean", diag = TRUE)

dist_G <- dist(gen_time_summary[,3], method = "euclidean", diag = TRUE)
dist_G <- dist(longev_summary[,3], method = "euclidean", diag = TRUE)
dist_G <- dist(mean_life_expect_summary[,3], method = "euclidean", diag = TRUE)
dist_G <- dist(R0_summary[,3], method = "euclidean", diag = TRUE)
dist_G <- dist(seed_size, method = "euclidean", diag = TRUE)
dist_G <- dist( relatedness, method = "euclidean", diag = TRUE)




plot(dist_G, dist_e)


for_pca <- tibble(lambda_sd_diff[1:7,1], lambda_mean_diff[1:7,1],gen_time_summary[,3], longev_summary[,3], mean_life_expect_summary[,3], R0_summary[,3], relatedness , seed_size, subtribes)
pc <- prcomp(for_pca[,1:8],
             center = TRUE,
             scale = FALSE)
pc <- prcomp(for_pca[,1:2],
             center = FALSE,
             scale = FALSE)
pc
summary(pc)
biplot(pc)

plot( lambda_mean_diff[1:7,1], lambda_var_diff[1:7,1], col = as.factor(subtribes))

pc <- prcomp(traits_df[,2:13])

library(vegan)
nmds <- envfit(pc,traits_df[,3:13])

nmds <- envfit(pc,for_pca[,3:6], permu = 999)

biplot(pc)
plot(nmds)
