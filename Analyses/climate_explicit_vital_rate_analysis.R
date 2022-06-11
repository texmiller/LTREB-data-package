## Title: Climate explicit Grass endophyte vital rate models with a Bayesian framework
## Purpose: Assembles data lists and runs vital rate kernels written in STAN with mixed effects, 
## and does visualisation of posterior predictive checks
## Authors: Joshua and Tom
#############################################################

library(tidyverse)
library(rstan)
library(StanHeaders)
# library(shinystan)
library(bayesplot)
library(devtools)
library(moments)
library(gridExtra)

library(RColorBrewer)


invlogit<-function(x){exp(x)/(1+exp(x))}
logit = function(x) { log(x/(1-x)) }
Lkurtosis=function(x) log(kurtosis(x)); 


#############################################################################################
####### Data manipulation to prepare data for Stan models------------------
#############################################################################################

#  data are prepared in the endodemog_data_processing.R file, 
source("Analyses/endodemog_data_processing.R")


#############################################################################################
####### Preparing data lists for vital rate kernels ------------------
#############################################################################################
# I'm gonna start with just one species
## Clean up the main data frame for NA's, other small data entry errors
LTREB_data_forsurv <- LTREB_full %>% 
  # filter(species == "FESU", year_t1 < 2019) %>%
  filter(!is.na(surv_t1)) %>% 
  filter(!is.na(logsize_t)) %>% 
  filter(!is.na(endo_01)) %>%   # There are a few LOAR that don't have a plot level endo assigned
  filter(origin_01 == 1 & year_t != birth | origin_01 == 0) %>%    # filtering out first year germinants (including those that are bigger than 1 tiller)
  filter(!is.na(spei12)) # I think I can try to get more recent climate data, but there are NA's for the 2020 and 2019 AGPE data right now
dim(LTREB_data_forsurv)

LTREB_surv_means <- LTREB_data_forsurv %>% 
  group_by(species, endo_01, year_t1) %>% 
  summarize(mean_surv = mean(surv_t1),
            spei12 = mean(as.numeric(spei12)),
            annual_precip = mean(as.numeric(annual_precip)),
            annual_temp = mean(as.numeric(annual_temp)),
            count = n())
# 
# ggplot(data = LTREB_surv_means)+
#   geom_point(aes(x = year_t1, y = mean_surv, color = as.factor(endo_01), size = count))+ facet_wrap(~species)+
#   theme_classic()
# # 
# ggplot(data = LTREB_surv_means)+
#   geom_point(aes(x = year_t1, y = spei12, color = endo_01, size = count))+ facet_wrap(~species)
# 
# ggplot(data = LTREB_surv_means)+
#   geom_point(aes(x = year_t1, y = annual_precip, color = endo_01, size = count))+ facet_wrap(~species)

# ggplot(data = LTREB_surv_means)+
#   geom_point(aes(x = year_t1, y = annual_temp, color = endo_01, size = count))+ facet_wrap(~species)
# 
# ggplot(data = LTREB_surv_means)+
#   geom_point(aes(x = spei12, y = mean_surv, color = as.factor(endo_01), size = count), alpha = .6)+ facet_wrap(~species) +
#   theme_classic()
# 
# ggplot(data = LTREB_surv_means)+
#   geom_point(aes(x = annual_precip, y = mean_surv, color = as.factor(endo_01), size = count), alpha = .6)+ facet_wrap(~species) +
#   theme_classic()
# 
# ggplot(data = LTREB_surv_means)+
#   geom_point(aes(x = annual_temp, y = mean_surv, color = as.factor(endo_01), size = count), alpha = .6)+ facet_wrap(~species) +
#   theme_classic()
# 
# 
# # 
LTREB_surv_seedling <- LTREB_full %>%
  filter(!is.na(surv_t1)) %>%
  filter(!is.na(logsize_t)) %>%
  filter(!is.na(endo_01)) %>%  # There are a few LOAR that don't have a plot level endo assigned
  filter(origin_01 == 1 & year_t == birth) %>%  #filtering for recruits that are just germinated
  filter(logsize_t == 0) # this is filtering out the plants that are "recruits" but are larger than 1 tiller
dim(LTREB_surv_seedling)
# 
LTREB_seedlingsurv_means <- LTREB_surv_seedling %>%
  group_by(species, endo_01, year_t1) %>%
  summarize(mean_surv = mean(surv_t1),
            spei12 = mean(as.numeric(spei12)),
            annual_precip = mean(as.numeric(annual_precip)),
            annual_temp = mean(as.numeric(annual_temp)),
            count = n())

# ggplot(data = LTREB_seedlingsurv_means)+
#   geom_point(aes(x = (year_t1), y = mean_surv, color = species, size = count))+ facet_wrap(~species) +
#   theme_classic()

# ggplot(data = LTREB_seedlingsurv_means)+
#   geom_point(aes(x = year_t1, y = spei12, color = species, size = count))+ facet_wrap(~species)
# 
# ggplot(data = LTREB_seedlingsurv_means)+
#   geom_point(aes(x = year_t1, y = spei12, color = endo_01, size = count))+ facet_wrap(~species)
# 
# ggplot(data = LTREB_seedlingsurv_means)+
#   geom_point(aes(x = year_t1, y = annual_precip, color = endo_01, size = count))+ facet_wrap(~species)
# 
# ggplot(data = LTREB_seedlingsurv_means)+
#   geom_point(aes(x = year_t1, y = annual_temp, color = endo_01, size = count))+ facet_wrap(~species)
# 
# 
# ggplot(data = LTREB_seedlingsurv_means)+
#   geom_point(aes(x = spei12, y = mean_surv, color = as.factor(endo_01), size = count))+ facet_wrap(~species) +
#   theme_classic()


# ggplot(data = LTREB_seedlingsurv_means)+
#   geom_point(aes(x = annual_precip, y = mean_surv, color = as.factor(endo_01), size = count), alpha = .6)+ facet_wrap(~species) +
#   theme_classic()
# 
# ggplot(data = LTREB_seedlingsurv_means)+
#   geom_point(aes(x = annual_temp, y = mean_surv, color = as.factor(endo_01), size = count), alpha = .6)+ facet_wrap(~species) +
#   theme_classic()


# 
# # I want to look at these "seedlings" that are bigger than 1 tiller, We are dropping these from our seedling model, as they are most likely missed plants from the previous year
# LTREB_surv_big_seedling <- LTREB_full %>% 
#   filter(!is.na(surv_t1)) %>% 
#   filter(!is.na(logsize_t)) %>% 
#   filter(!is.na(endo_01)) %>%  # There are a few LOAR that don't have a plot level endo assigned
#   filter(origin_01 == 1 & year_t == birth) %>%  #filtering for recruits that are just germinated
#   filter(logsize_t != 0) 
# dim(LTREB_surv_big_seedling)
# table(LTREB_surv_big_seedling$year_t,LTREB_surv_big_seedling$size_t, LTREB_surv_big_seedling$species)
# 
LTREB_data_forflw <- LTREB_full %>%
  filter(!is.na(FLW_STAT_T1)) %>%
  filter(!is.na(logsize_t1)) %>%
  filter(!is.na(endo_01))
dim(LTREB_data_forflw)


LTREB_flw_means <- LTREB_data_forflw %>% 
  group_by(species, endo_01, year_t1) %>% 
  summarize(mean_flw = mean(FLW_STAT_T),
            spei12 = mean(as.numeric(spei12)),
            annual_precip = mean(as.numeric(annual_precip)),
            annual_temp = mean(as.numeric(annual_temp)),
            count = n())

# ggplot(data = LTREB_flw_means)+
#   geom_point(aes(x = year_t1, y = mean_flw, color = as.factor(endo_01), size = count))+ facet_wrap(~species)
# 
# ggplot(data = LTREB_flw_means)+
#   geom_point(aes(x = year_t1, y = spei12, color = species, size = count))+ facet_wrap(~species)
# 
# ggplot(data = LTREB_flw_means)+
#   geom_point(aes(x = year_t1, y = spei12, color = endo_01, size = count))+ facet_wrap(~species)
# 
# ggplot(data = LTREB_flw_means)+
#   geom_point(aes(x = year_t1, y = annual_precip, color = endo_01, size = count))+ facet_wrap(~species)
# 
# ggplot(data = LTREB_flw_means)+
#   geom_point(aes(x = year_t1, y = annual_temp, color = endo_01, size = count))+ facet_wrap(~species)
# 
# 
# ggplot(data = LTREB_flw_means)+
#   geom_point(aes(x = spei12, y = mean_flw, color = as.factor(endo_01), size = count))+ facet_wrap(~species) +
#   theme_classic()
# 
# ggplot(data = LTREB_flw_means)+
#   geom_point(aes(x = annual_precip, y = mean_flw, color = as.factor(endo_01), size = count))+ facet_wrap(~species) +
#   theme_classic()
# 
# ggplot(data = LTREB_flw_means)+
#   geom_point(aes(x = annual_temp, y = mean_flw, color = as.factor(endo_01), size = count))+ facet_wrap(~species) +
#   theme_classic()

# 
LTREB_data_forgrow <- LTREB_full %>%
  filter(!is.na(logsize_t)) %>%
  filter(!is.na(size_t1)) %>%
  filter(!is.na(endo_01)) %>%
  filter(origin_01 == 1 & year_t != birth | origin_01 == 0)  # filtering out first year germinants (including those that are bigger than 1 tiller)
dim(LTREB_data_forgrow)

# 
# LTREB_grow_seedling <- LTREB_full %>%
#   filter(!is.na(logsize_t)) %>% 
#   filter(!is.na(size_t1)) %>% 
#   filter(!is.na(endo_01)) %>% 
#   filter(origin_01 == 1 & year_t == birth) %>%  #filtering for recruits that are just germinated
#   filter(logsize_t == 0) # this is filtering out the plants that are "recruits" but are larger than 1 tiller
# 
# dim(LTREB_grow_seedling)
# 
# 
# LTREB_data_forfert <- LTREB_full %>% 
#   filter(!is.na(FLW_COUNT_T)) %>% 
#   filter(FLW_COUNT_T > 0) %>% 
#   filter(!is.na(logsize_t))
# dim(LTREB_data_forfert)
# 
# LTREB_data_forspike <- LTREB_full %>%
#   dplyr::select(-FLW_COUNT_T1, -FLW_STAT_T1, -SPIKE_A_T1, -SPIKE_B_T1, -SPIKE_C_T1, -SPIKE_D_T1, -SPIKE_AGPE_MEAN_T1, -endo_status_from_check, -plot_endo_for_check, -endo_mismatch, -dist_a, -dist_b) %>% 
#   filter(!is.na(FLW_STAT_T)) %>% 
#   filter(FLW_STAT_T>0) %>% 
#   melt(id.var = c("plot_fixed" ,   "plot_index",         "pos"         ,           "id",
#                   "species"       ,         "species_index"  ,        "endo_01",
#                   "endo_index"  ,           "origin_01"       ,       "birth" ,
#                   "year_t1"         ,       "year_t1_index"       ,   "surv_t1" ,
#                   "size_t1"         ,       "logsize_t1"       ,
#                   "year_t",
#                   "year_t_index"     ,      "size_t"           ,      "logsize_t"  ,
#                   "FLW_COUNT_T"      ,      "FLW_STAT_T"),
#        value.name = "spike_count_t") %>% 
#   rename(spikelet_id = variable) %>% 
#   filter(!is.na(spike_count_t), spike_count_t > 0) %>% 
#   mutate(spike_count_t = as.integer(spike_count_t))
# 
# # ggplot(LTREB_data_forspike)+
#   geom_histogram(aes(x=spike_count_t))+
#   facet_grid(year_t~species)
## I don't think there are enough data to fit year variances
## so I am just going to fit fixed effects of size and endo

# rm(LTREB_full)


# Create data lists to be used for the Stan model

surv_data_list <- list(y = LTREB_data_forsurv$surv_t1,
                       spei = as.numeric(LTREB_data_forsurv$spei12),
                       spei_nl = as.numeric(LTREB_data_forsurv$spei12^2),
                       logsize = LTREB_data_forsurv$logsize_t,
                       origin_01 = LTREB_data_forsurv$origin_01,
                       endo_01 = as.integer(LTREB_data_forsurv$endo_01),
                       endo_index = as.integer(LTREB_data_forsurv$endo_index),
                       spp = as.integer(LTREB_data_forsurv$species_index),
                       year_t = as.integer(LTREB_data_forsurv$year_t_index),
                       plot = as.integer(as.factor(LTREB_data_forsurv$plot_index)),
                       N = nrow(LTREB_data_forsurv),
                       nSpp = length(unique(LTREB_data_forsurv$species_index)),
                       nYear = max(unique(LTREB_data_forsurv$year_t_index)),
                       nPlot = length(unique(LTREB_data_forsurv$plot_index)),
                       nEndo =   length(unique(LTREB_data_forsurv$endo_01)))
str(surv_data_list)

seed_surv_data_list <- list(y = LTREB_surv_seedling$surv_t1,
                            spei = as.numeric(LTREB_surv_seedling$spei12),
                            spei_nl = as.numeric(LTREB_surv_seedling$spei12^2),
                            logsize_t = LTREB_surv_seedling$logsize_t,
                            origin_01 = LTREB_surv_seedling$origin_01,
                            endo_01 = as.integer(LTREB_surv_seedling$endo_01),
                            endo_index = as.integer(LTREB_surv_seedling$endo_index),
                            spp = as.integer(LTREB_surv_seedling$species_index),
                            year_t = as.integer(LTREB_surv_seedling$year_t_index),
                            plot = as.integer(LTREB_surv_seedling$plot_index),
                            N = nrow(LTREB_surv_seedling),
                            nSpp = length(unique(LTREB_surv_seedling$species_index)),
                            nYear = max(unique(LTREB_surv_seedling$year_t_index)),
                            nPlot = length(unique(LTREB_data_forsurv$plot_index)),
                            nEndo =   length(unique(LTREB_surv_seedling$endo_01)))
str(seed_surv_data_list)



flw_data_list <- list(y = LTREB_data_forflw$FLW_STAT_T1,
                      spei = as.numeric(LTREB_data_forflw$spei12),
                      spei_nl = as.numeric(LTREB_data_forflw$spei12^2),
                      logsize = LTREB_data_forflw$logsize_t1,
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
# 
grow_data_list <- list(y = as.integer(LTREB_data_forgrow$size_t1),
                       spei = as.numeric(LTREB_data_forgrow$spei12),
                       spei_nl = as.numeric(LTREB_data_forgrow$spei12^2),
                       logsize = LTREB_data_forgrow$logsize_t,
                       origin_01 = as.integer(LTREB_data_forgrow$origin_01),
                       endo_01 = as.integer(LTREB_data_forgrow$endo_01),
                       endo_index = as.integer(LTREB_data_forgrow$endo_index),
                       spp = as.integer(LTREB_data_forgrow$species_index),
                       year_t = as.integer(LTREB_data_forgrow$year_t_index),
                       plot = as.integer(LTREB_data_forgrow$plot_index),
                       N = nrow(LTREB_data_forgrow),
                       nSpp = length(unique(LTREB_data_forgrow$species_index)),
                       nYear = max(unique(LTREB_data_forgrow$year_t_index)),
                       nPlot = length(unique(LTREB_data_forgrow$plot_index)),
                       nEndo =   length(unique(LTREB_data_forgrow$endo_01)))
str(grow_data_list)
# seed_grow_data_list <- list(y = as.integer(LTREB_grow_seedling$size_t1),
#                             logsize_t = LTREB_grow_seedling$logsize_t,
#                             origin_01 = as.integer(LTREB_grow_seedling$origin_01),
#                             endo_01 = as.integer(LTREB_grow_seedling$endo_01),
#                             endo_index = as.integer(LTREB_grow_seedling$endo_index),
#                             spp = as.integer(LTREB_grow_seedling$species_index),
#                             year_t = as.integer(LTREB_grow_seedling$year_t_index),
#                             plot = as.integer(LTREB_grow_seedling$plot_index),
#                             N = nrow(LTREB_grow_seedling),
#                             nSpp = length(unique(LTREB_grow_seedling$species_index)),
#                             nYear = max(unique(LTREB_grow_seedling$year_t_index)),
#                             nPlot = max(unique(LTREB_grow_seedling$plot_index)),
#                             nEndo =   length(unique(LTREB_grow_seedling$endo_01)))
# str(seed_grow_data_list)
# 
# 
# 
# fert_data_list <- list(y = as.integer(LTREB_data_forfert$FLW_COUNT_T),
#                        logsize_t = LTREB_data_forfert$logsize_t,
#                        origin_01 = LTREB_data_forfert$origin_01,
#                        endo_01 = as.integer(LTREB_data_forfert$endo_01),
#                        endo_index = as.integer(LTREB_data_forfert$endo_index),
#                        spp = as.integer(LTREB_data_forfert$species_index),
#                        year_t = as.integer(LTREB_data_forfert$year_t_index),
#                        plot = as.integer(LTREB_data_forfert$plot_index),
#                        N = nrow(LTREB_data_forfert),
#                        nSpp = length(unique(LTREB_data_forfert$species_index)),
#                        nYear = max(unique(LTREB_data_forfert$year_t_index)),
#                        nPlot = max(unique(LTREB_data_forfert$plot_index)),
#                        nEndo =   length(unique(LTREB_data_forfert$endo_01)))
# str(fert_data_list)
# 
# 
# 
# spike_data_list <- list(nYear = max(unique(LTREB_data_forspike$year_t_index)),
#                         nPlot = max(unique(LTREB_data_forspike$plot_index)),
#                         nSpp = length(unique(LTREB_data_forspike$species)),
#                         nEndo=length(unique(LTREB_data_forspike$endo_01)),
#                         N = nrow(LTREB_data_forspike),
#                         year_t = as.integer(LTREB_data_forspike$year_t_index),
#                         plot = as.integer(LTREB_data_forspike$plot_index),
#                         spp = as.integer(as.numeric(as.factor(LTREB_data_forspike$species))),
#                         y = LTREB_data_forspike$spike_count_t,
#                         logsize_t = LTREB_data_forspike$logsize_t,
#                         endo_01 = LTREB_data_forspike$endo_01,
#                         origin_01 = LTREB_data_forspike$origin_01)
# str(spike_data_list)

#########################################################################################################
# Stan model runs------------------------------
#########################################################################################################
## run this code to optimize computer system settings for MCMC
rstan_options(auto_write = FALSE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## MCMC settings
mcmc_pars <- list(
  warmup = 2500, 
  iter = 5000, 
  thin = 1, 
  chains = 3
)

sm_surv <- stan(file = "Analyses/climate_endo_spp_surv_flw.stan", data = surv_data_list,
                iter = mcmc_pars$iter,
                warmup = mcmc_pars$warmup,
                chains = mcmc_pars$chains, 
                thin = mcmc_pars$thin)
# saveRDS(sm_surv, file = "~/Dropbox/EndodemogData/Model_Runs/climate_endo_spp_surv_woseedling.rds")
# The model with squared climate terrms gives max_treedepth warnings but otherwise converges okay but without non-linear terms does not have wraning.

#Running the seedling survival model for 10000 iterations running for only 5000 led to low effective sample size warnings, most likely in sigma0 parameter
# running for 10000 iterations had enough effective sample size, but had 3 divergent transitions. Could try re-running and hope that goes away, or try simplifying, like no non-linear term
# running without non-linear term: model fits with no errors
# running with tighter priors for climate parameters: fits with no errors
sm_seed_surv <- stan(file = "Analyses/climate_seedling_surv.stan", data = seed_surv_data_list,
                     iter = mcmc_pars$iter,
                     warmup = mcmc_pars$warmup,
                     chains = mcmc_pars$chains, 
                     thin = mcmc_pars$thin)
saveRDS(sm_seed_surv, file = "~/Dropbox/EndodemogData/Model_Runs/climate_endo_seedling_surv_withtightpriors.rds")

sm_seed_surv_linear <- readRDS(file = "~/Dropbox/EndodemogData/Model_Runs/climate_endo_seedling_surv_withoutquadraticterm.rds")
sm_seed_surv_nonlinear <- readRDS(file = "~/Dropbox/EndodemogData/Model_Runs/climate_endo_seedling_surv_withtightpriors.rds")

quote_bare <- function( ... ){
  substitute( alist(...) ) %>% 
    eval( ) %>% 
    sapply( deparse )
}
surv_sdlg_parlinear <- rstan::extract(sm_seed_surv_linear, pars =quote_bare(beta0,betaendo,betaspei_endo,
                                                                tau_year, tau_plot))
surv_sdlg_parnonlinear <- rstan::extract(sm_seed_surv_nonlinear, pars =quote_bare(beta0,betaendo,betaspei_endo,betaspei_nl_endo,
                                                                     tau_year, tau_plot))
linear_s_seed_eminus <- matrix(NA, 7,100)
linear_s_seed_eplus <- matrix(NA, 7,100)

nonlinear_s_seed_eminus <- matrix(NA, 7,100)
nonlinear_s_seed_eplus <- matrix(NA, 7,100)

for(s in 1:7){
linear_s_seed_eminus[s,] <- invlogit(mean(surv_sdlg_parlinear$beta0[,s])+mean(surv_sdlg_parlinear$betaendo[,s])*0+mean(surv_sdlg_parlinear$betaspei_endo[,s,1])*seq(-1,2, length.out = 100))#+mean(surv_sdlg_parlinear$betaspei_nl_endo[,s,1]*(seq(-1,2, length.out = 100)^2)))
linear_s_seed_eplus[s,] <- invlogit(mean(surv_sdlg_parlinear$beta0[,s])+mean(surv_sdlg_parlinear$betaendo[,s])*1+mean(surv_sdlg_parlinear$betaspei_endo[,s,2])*seq(-1,2, length.out = 100))#+mean(surv_sdlg_parlinear$betaspei_nl_endo[,s,2]*(seq(-1,2, length.out = 100)^2)))
nonlinear_s_seed_eminus[s,] <- invlogit(mean(surv_sdlg_parnonlinear$beta0[,s])+mean(surv_sdlg_parnonlinear$betaendo[,s])*0+mean(surv_sdlg_parnonlinear$betaspei_endo[,s,1])*seq(-1,2, length.out = 100)+mean(surv_sdlg_parnonlinear$betaspei_nl_endo[,s,1])*(seq(-1,2, length.out = 100)^2))
nonlinear_s_seed_eplus[s,] <- invlogit(mean(surv_sdlg_parnonlinear$beta0[,s])+mean(surv_sdlg_parnonlinear$betaendo[,s])*1+mean(surv_sdlg_parnonlinear$betaspei_endo[,s,2])*seq(-1,2, length.out = 100)+mean(surv_sdlg_parnonlinear$betaspei_nl_endo[,s,2])*(seq(-1,2, length.out = 100)^2))
}


# plotting the predicted linear and nonlinear fits vs the seedling survival data
pred <- as_tibble(t(linear_s_seed_eminus)) %>% 
  rename("lin_Eminus_"= contains("V")) %>% 
  cbind(spei = seq(-1,2,length.out=100)) %>% 
  pivot_longer(cols = contains("lin")) %>% 
  separate(name, c("effect", "endo", "species"))

pred2 <- as_tibble(t(linear_s_seed_eplus)) %>% 
  rename("lin_Eplus_"= contains("V")) %>% 
  cbind(spei = seq(-1,2,length.out=100)) %>% 
  pivot_longer(cols = contains("lin")) %>% 
  separate(name, c("effect", "endo", "species"))

pred3 <- as_tibble(t(nonlinear_s_seed_eminus)) %>% 
  rename("nonlin_Eminus_"= contains("V")) %>% 
  cbind(spei = seq(-1,2,length.out=100)) %>% 
  pivot_longer(cols = contains("nonlin")) %>% 
  separate(name, c("effect", "endo", "species"))

pred4 <- as_tibble(t(nonlinear_s_seed_eplus)) %>% 
  rename("nonlin_Eplus_"= contains("V")) %>% 
  cbind(spei = seq(-1,2,length.out=100)) %>% 
  pivot_longer(cols = contains("nonlin")) %>% 
  separate(name, c("effect", "endo", "species"))

pred0 <- pred %>% 
  full_join(pred2) %>% 
  full_join(pred3) %>% 
  full_join(pred4) %>% 
  mutate(species = case_when(species == 1 ~ "AGPE",
                             species == 2 ~ "ELRI",
                             species == 3 ~ "ELVI",
                             species == 4 ~ "FESU",
                             species == 5 ~ "LOAR",
                             species == 6 ~ "POAL",
                             species == 7 ~ "POSY"))

ggplot(data = LTREB_seedlingsurv_means)+
  geom_point(aes(x = spei12, y = mean_surv, color = species, shape = as.factor(endo_01), size = count))+ 
  geom_line(data = pred0, aes(x = spei, y = value, color = species, lty = endo))+
  facet_wrap(~species+effect) +
  scale_shape_manual(values = c(16,1))+
  theme_classic()

sm_flw <- stan(file = "Analyses/climate_endo_spp_surv_flw.stan", data = flw_data_list,
                iter = mcmc_pars$iter,
                warmup = mcmc_pars$warmup,
                chains = mcmc_pars$chains, 
                thin = mcmc_pars$thin)
saveRDS(sm_flw, file = "~/Dropbox/EndodemogData/Model_Runs/climate_endo_spp_flw.rds")
# gave the same max_treedepth warning, but otherwise ran okay


# Running the growth model with the PIG with climate effects
sm_grow <- stan(file = "Analyses/climate_endo_spp_grow_fert_PIG.stan", data = grow_data_list,
               iter = mcmc_pars$iter,
               warmup = mcmc_pars$warmup,
               chains = mcmc_pars$chains, 
               thin = mcmc_pars$thin)
saveRDS(sm_grow, file = "~/Dropbox/EndodemogData/Model_Runs/climate_endo_grow_fert_PIG.rds")
# gave the same max_treedepth warning, but otherwise ran okay


#########################################################################################################
# Model Diagnostics ------------------------------
#########################################################################################################
# Function for looking at binned size_t fits, particularly important for the growth kernel as this determines the transitions through the matrix model
size_moments_ppc <- function(data,y_name,sim, n_bins, title = NA){
  require(tidyverse)
  data$y_name <- data[[y_name]]
  bins <- data %>%
    ungroup() %>% 
    arrange(logsize_t) %>% 
    mutate(size_bin = cut_number(logsize_t, n_bins)) %>% 
    group_by(size_bin)  %>% 
    dplyr::summarize(mean_t1 = mean(y_name),
                     sd_t1 = sd(y_name),
                     skew_t1 = skewness(y_name),
                     kurt_t1 = Lkurtosis(y_name),
                     bin_mean = mean(logsize_t),
                     bin_n = n())
  sim_moments <- bind_cols(enframe(data$logsize_t), as_tibble(t(sim))) %>%
    rename(logsize_t = value) %>%
    arrange(logsize_t) %>%
    mutate(size_bin = cut_number(logsize_t, n_bins)) %>%
    pivot_longer(., cols = starts_with("V"), names_to = "post_draw", values_to = "sim") %>%
    group_by(size_bin, post_draw) %>%
    summarize( mean_sim = mean((sim)),
               sd_sim = sd((sim)),
               skew_sim = skewness((sim)),
               kurt_sim = Lkurtosis((sim)),
               bin_mean = mean(logsize_t),
               bin_n = n())
  sim_medians <- sim_moments %>%
    group_by(size_bin, bin_mean) %>%
    summarize(median_mean_sim = median(mean_sim),
              median_sd_sim = median(sd_sim),
              median_skew_sim = median(skew_sim),
              median_kurt_sim = median(kurt_sim))
  meanplot <-  ggplot(data = bins)+
    geom_point(data = sim_moments, aes(x = bin_mean, y = mean_sim), color = "gray72") +
    geom_point(data = sim_medians, aes(x = bin_mean, y = median_mean_sim),shape = 1, color = "black") +
    geom_point(aes(x = bin_mean, y = mean_t1), shape = 1, color = "firebrick2") +
    theme_classic()
  sdplot <-  ggplot(data = bins)+
    geom_point(data = sim_moments, aes(x = bin_mean, y = sd_sim), color = "gray72") +
    geom_point(data = sim_medians, aes(x = bin_mean, y = median_sd_sim),shape = 1, color = "black") +
    geom_point(aes(x = bin_mean, y = sd_t1), shape = 1, color = "firebrick2") + theme_classic()
  skewplot <-  ggplot(data = bins)+
    geom_point(data = sim_moments, aes(x = bin_mean, y = skew_sim), color = "gray72") +
    geom_point(data = sim_medians, aes(x = bin_mean, y = median_skew_sim),shape = 1, color = "black") +
    geom_point(aes(x = bin_mean, y = skew_t1), shape = 1, color = "firebrick2") + theme_classic()
  kurtplot <- ggplot(data = bins)+
    geom_point(data = sim_moments, aes(x = bin_mean, y = kurt_sim), color = "gray72") +
    geom_point(data = sim_medians, aes(x = bin_mean, y = median_kurt_sim),shape = 1, color = "black") +
    geom_point(aes(x = bin_mean, y = kurt_t1), shape = 1, color = "firebrick2") + theme_classic()
  size_ppc_plot <- grid.arrange(meanplot, sdplot,skewplot, kurtplot, top = title)
  return(size_ppc_plot)
}


#### survival ppc ####
surv_fit <- sm_surv
# surv_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_spp_surv_woseedling.rds")
predS <- rstan::extract(surv_fit, pars = c("p"))$p # extract the linear predictor
n_post_draws <- 100
post_draws <- sample.int(dim(predS)[1], n_post_draws) # draw samples from the posterior of the linear predictor
y_s_sim <- matrix(NA,n_post_draws,length(surv_data_list$y))
for(i in 1:n_post_draws){
  y_s_sim[i,] <- rbinom(n=length(surv_data_list$y), size=1, prob = invlogit(predS[post_draws[i],]))
}
# ppc_dens_overlay(surv_data_list$y, y_s_sim)
surv_densplot <- ppc_dens_overlay(surv_data_list$y, y_s_sim) + theme_classic() + labs(title = "Adult Survival", x = "Survival status", y = "Density")
surv_densplot
# ggsave(surv_densplot, filename = "surv_densplot.png", width = 4, height = 4)

mean_s_plot <-   ppc_stat(surv_data_list$y, y_s_sim, stat = "mean")
sd_s_plot <- ppc_stat(surv_data_list$y, y_s_sim, stat = "sd")
skew_s_plot <- ppc_stat(surv_data_list$y, y_s_sim, stat = "skewness")
kurt_s_plot <- ppc_stat(surv_data_list$y, y_s_sim, stat = "Lkurtosis")
grid.arrange(mean_s_plot,sd_s_plot,skew_s_plot,kurt_s_plot,  top = "Survival")


# now we want to look at how the the model is fitting across sizes

surv_size_ppc <- size_moments_ppc(data = LTREB_data_forsurv,
                                  y_name = "surv_t1",
                                  sim = y_s_sim, 
                                  n_bins = 4, 
                                  title = "Survival")

