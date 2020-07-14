## Title: Grass endophyte population model with a bayesian framework
## Purpose: Assembles data lists and runs vital rate kernels written in STAN with mixed effects, 
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
####### Data manipulation to prepare data for Stan models------------------
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
  filter(!is.na(endo_01)) %>%   # There are a few LOAR that don't have a plot level endo assigned
  filter(origin_01 == 1 & year_t != birth | origin_01 == 0)  # filtering out first year germinants (including those that are bigger than 1 tiller)
dim(LTREB_data_forsurv)

LTREB_surv_seedling <- LTREB_full %>% 
  filter(!is.na(surv_t1)) %>% 
  filter(!is.na(logsize_t)) %>% 
  filter(!is.na(endo_01)) %>%  # There are a few LOAR that don't have a plot level endo assigned
  filter(origin_01 == 1 & year_t == birth) %>%  #filtering for recruits that are just germinated
  filter(logsize_t == 0) # this is filtering out the plants that are "recruits" but are larger than 1 tiller
dim(LTREB_surv_seedling)

# I want to look at these "seedlings" that are bigger than 1 tiller, We are dropping these from our seedling model, as they are most likely missed plants from the previous year
LTREB_surv_big_seedling <- LTREB_full %>% 
  filter(!is.na(surv_t1)) %>% 
  filter(!is.na(logsize_t)) %>% 
  filter(!is.na(endo_01)) %>%  # There are a few LOAR that don't have a plot level endo assigned
  filter(origin_01 == 1 & year_t == birth) %>%  #filtering for recruits that are just germinated
  filter(logsize_t != 0) 
dim(LTREB_surv_big_seedling)
table(LTREB_surv_big_seedling$year_t,LTREB_surv_big_seedling$size_t, LTREB_surv_big_seedling$species)

LTREB_data_forflw <- LTREB_full %>% 
  filter(!is.na(FLW_STAT_T)) %>% 
  filter(!is.na(logsize_t)) %>% 
  filter(!is.na(endo_01))
dim(LTREB_data_forflw)


LTREB_data_forgrow <- LTREB_full %>%
  filter(!is.na(logsize_t)) %>% 
  filter(!is.na(size_t1)) %>% 
  filter(!is.na(endo_01)) %>% 
  filter(origin_01 == 1 & year_t != birth | origin_01 == 0)  # filtering out first year germinants (including those that are bigger than 1 tiller)
dim(LTREB_data_forgrow)

LTREB_grow_seedling <- LTREB_full %>%
  filter(!is.na(logsize_t)) %>% 
  filter(!is.na(size_t1)) %>% 
  filter(!is.na(endo_01)) %>% 
  filter(origin_01 == 1 & year_t == birth) %>%  #filtering for recruits that are just germinated
  filter(logsize_t == 0) # this is filtering out the plants that are "recruits" but are larger than 1 tiller

dim(LTREB_grow_seedling)


LTREB_data_forfert <- LTREB_full %>% 
  filter(!is.na(FLW_COUNT_T)) %>% 
  filter(FLW_COUNT_T > 0) %>% 
  filter(!is.na(logsize_t))
dim(LTREB_data_forfert)

LTREB_data_forspike <- LTREB_full %>%
  dplyr::select(-FLW_COUNT_T1, -FLW_STAT_T1, -SPIKE_A_T1, -SPIKE_B_T1, -SPIKE_C_T1, -SPIKE_D_T1, -SPIKE_AGPE_MEAN_T1, -endo_status_from_check, -plot_endo_for_check, -endo_mismatch, -dist_a, -dist_b) %>% 
  filter(!is.na(FLW_STAT_T)) %>% 
  filter(FLW_STAT_T>0) %>% 
  melt(id.var = c("plot_fixed" ,   "plot_index",         "pos"         ,           "id",
                  "species"       ,         "species_index"  ,        "endo_01",
                  "endo_index"  ,           "origin_01"       ,       "birth" ,
                  "year_t1"         ,       "year_t1_index"       ,   "surv_t1" ,
                  "size_t1"         ,       "logsize_t1"       ,
                  "year_t",
                  "year_t_index"     ,      "size_t"           ,      "logsize_t"  ,
                  "FLW_COUNT_T"      ,      "FLW_STAT_T"),
       value.name = "spike_count_t") %>% 
  rename(spikelet_id = variable) %>% 
  filter(!is.na(spike_count_t), spike_count_t > 0) %>% 
  mutate(spike_count_t = as.integer(spike_count_t))

# ggplot(LTREB_data_forspike)+
#   geom_histogram(aes(x=spike_count_t))+
#   facet_grid(year_t~species)
## I don't think there are enough data to fit year variances
## so I am just going to fit fixed effects of size and endo

# rm(LTREB_full)


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
                       nPlot = max(unique(LTREB_data_forsurv$plot_index)),
                       nEndo =   length(unique(LTREB_data_forsurv$endo_01)))
str(surv_data_list)

seed_surv_data_list <- list(y = LTREB_surv_seedling$surv_t1,
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
                            nPlot = max(unique(LTREB_surv_seedling$plot_index)),
                            nEndo =   length(unique(LTREB_surv_seedling$endo_01)))
str(seed_surv_data_list)


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

grow_data_list <- list(y = as.integer(LTREB_data_forgrow$size_t1),
                       logsize_t = LTREB_data_forgrow$logsize_t,
                       origin_01 = as.integer(LTREB_data_forgrow$origin_01),
                       endo_01 = as.integer(LTREB_data_forgrow$endo_01),
                       endo_index = as.integer(LTREB_data_forgrow$endo_index),
                       spp = as.integer(LTREB_data_forgrow$species_index),
                       year_t = as.integer(LTREB_data_forgrow$year_t_index),
                       plot = as.integer(LTREB_data_forgrow$plot_index),
                       N = nrow(LTREB_data_forgrow),
                       nSpp = length(unique(LTREB_data_forgrow$species_index)),
                       nYear = max(unique(LTREB_data_forgrow$year_t_index)),
                       nPlot = max(unique(LTREB_data_forgrow$plot_index)),
                       nEndo =   length(unique(LTREB_data_forgrow$endo_01)))
str(grow_data_list)
seed_grow_data_list <- list(y = as.integer(LTREB_grow_seedling$size_t1),
                       logsize_t = LTREB_grow_seedling$logsize_t,
                       origin_01 = as.integer(LTREB_grow_seedling$origin_01),
                       endo_01 = as.integer(LTREB_grow_seedling$endo_01),
                       endo_index = as.integer(LTREB_grow_seedling$endo_index),
                       spp = as.integer(LTREB_grow_seedling$species_index),
                       year_t = as.integer(LTREB_grow_seedling$year_t_index),
                       plot = as.integer(LTREB_grow_seedling$plot_index),
                       N = nrow(LTREB_grow_seedling),
                       nSpp = length(unique(LTREB_grow_seedling$species_index)),
                       nYear = max(unique(LTREB_grow_seedling$year_t_index)),
                       nPlot = max(unique(LTREB_grow_seedling$plot_index)),
                       nEndo =   length(unique(LTREB_grow_seedling$endo_01)))
str(seed_grow_data_list)



fert_data_list <- list(y = as.integer(LTREB_data_forfert$FLW_COUNT_T),
                       logsize_t = LTREB_data_forfert$logsize_t,
                       origin_01 = LTREB_data_forfert$origin_01,
                       endo_01 = as.integer(LTREB_data_forfert$endo_01),
                       endo_index = as.integer(LTREB_data_forfert$endo_index),
                       spp = as.integer(LTREB_data_forfert$species_index),
                       year_t = as.integer(LTREB_data_forfert$year_t_index),
                       plot = as.integer(LTREB_data_forfert$plot_index),
                       N = nrow(LTREB_data_forfert),
                       nSpp = length(unique(LTREB_data_forfert$species_index)),
                       nYear = max(unique(LTREB_data_forfert$year_t_index)),
                       nPlot = max(unique(LTREB_data_forfert$plot_index)),
                       nEndo =   length(unique(LTREB_data_forfert$endo_01)))
str(fert_data_list)



spike_data_list <- list(nYear = max(unique(LTREB_data_forspike$year_t_index)),
                        nPlot = max(unique(LTREB_data_forspike$plot_index)),
                        nSpp = length(unique(LTREB_data_forspike$species)),
                        nEndo=length(unique(LTREB_data_forspike$endo_01)),
                        N = nrow(LTREB_data_forspike),
                        year_t = as.integer(LTREB_data_forspike$year_t_index),
                        plot = as.integer(LTREB_data_forspike$plot_index),
                        spp = as.integer(as.numeric(as.factor(LTREB_data_forspike$species))),
                        y = LTREB_data_forspike$spike_count_t,
                        logsize_t = LTREB_data_forspike$logsize_t,
                        endo_01 = LTREB_data_forspike$endo_01,
                        origin_01 = LTREB_data_forspike$origin_01)
str(spike_data_list)

#########################################################################################################
# Stan model runs------------------------------
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

sm_surv <- stan(file = "Analyses/endo_spp_surv_flw.stan", data = surv_data_list,
                iter = mcmc_pars$iter,
                warmup = mcmc_pars$warmup,
                chains = mcmc_pars$chains, 
                thin = mcmc_pars$thin,
                control = list(adapt_delta = .9))
# saveRDS(sm_surv, file = "~/Dropbox/EndodemogData/Model_Runs/endo_spp_surv_woseedling.rds")

sm_seed_surv <- stan(file = "Analyses/seedling_surv.stan", data = seed_surv_data_list,
                     iter = mcmc_pars$iter,
                     warmup = mcmc_pars$warmup,
                     chains = mcmc_pars$chains, 
                     thin = mcmc_pars$thin,
                     control = list(adapt_delta = .9))
# saveRDS(sm_seed_surv, file = "~/Dropbox/EndodemogData/Model_Runs/endo_seedling_surv.rds")


sm_flw <- stan(file = "Analyses/endo_spp_surv_flw.stan", data = flw_data_list,
               iter = mcmc_pars$iter,
               warmup = mcmc_pars$warmup,
               chains = mcmc_pars$chains, 
               thin = mcmc_pars$thin)
# saveRDS(sm_flw, file = "~/Dropbox/EndodemogData/Model_Runs/endo_spp_flw.rds")

sm_grow <- stan(file = "Analyses/endo_spp_grow_fert.stan", data = grow_data_list,
                iter = mcmc_pars$iter,
                warmup = mcmc_pars$warmup,
                chains = mcmc_pars$chains, 
                thin = mcmc_pars$thin)
# saveRDS(sm_grow, file = "~/Dropbox/EndodemogData/Model_Runs/endo_spp_grow.rds")

sm_grow_pig <- stan(file = "Analyses/endo_spp_grow_fert_PIG.stan", data = grow_data_list,
                    iter = mcmc_pars$iter,
                    warmup = mcmc_pars$warmup,
                    chains = mcmc_pars$chains, 
                    thin = mcmc_pars$thin)
# saveRDS(sm_grow_pig, file = "~/Dropbox/EndodemogData/Model_Runs/endo_spp_grow_PIG.rds")

sm_seed_grow <- stan(file = "Analyses/seedling_grow.stan", data = seed_grow_data_list,
                iter = mcmc_pars$iter,
                warmup = mcmc_pars$warmup,
                chains = mcmc_pars$chains, 
                thin = mcmc_pars$thin,
                control = list(adapt_delta = .9))
# saveRDS(sm_seed_grow, file = "~/Dropbox/EndodemogData/Model_Runs/endo_seedling_grow.rds")



sm_fert <- stan(file = "Analyses/endo_spp_grow_fert.stan", data = fert_data_list,
               iter = mcmc_pars$iter,
               warmup = mcmc_pars$warmup,
               chains = mcmc_pars$chains, 
               thin = mcmc_pars$thin,
               control = list(adapt_delta = .99, max_treedepth = 15))
# saveRDS(sm_fert, file = "~/Dropbox/EndodemogData/Model_Runs/endo_spp_fert.rds")

sm_fert_nc <- stan(file = "Analyses/endo_spp_grow_fert_nc.stan", data = fert_data_list,
                iter = mcmc_pars$iter,
                warmup = mcmc_pars$warmup,
                chains = mcmc_pars$chains, 
                thin = mcmc_pars$thin)
# saveRDS(sm_fert, file = "~/Dropbox/EndodemogData/Model_Runs/endo_spp_fert_nc.rds")

sm_fert <- stan(file = "Analyses/endo_spp_grow_fert_noplot.stan", data = fert_data_list,
                iter = mcmc_pars$iter,
                warmup = mcmc_pars$warmup,
                chains = mcmc_pars$chains, 
                thin = mcmc_pars$thin,
                control = list(adapt_delta = .99)) # run without this, there is a small number of divergent transitions
# saveRDS(sm_fert, file = "~/Dropbox/EndodemogData/Model_Runs/endo_spp_fert_noplot.rds")



sm_spike <- stan(file = "Analyses/endo_spp_spike.stan", data = spike_data_list,
                 iter = mcmc_pars$iter,
                 warmup = mcmc_pars$warmup,
                 chains = mcmc_pars$chains, 
                 thin = mcmc_pars$thin)
saveRDS(sm_spike, file = "~/Dropbox/EndodemogData/Model_Runs/endo_spp_spike.rds")

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
surv_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_surv_woseedling.rds")
predS <- rstan::extract(surv_fit, pars = c("p"))$p
n_post_draws <- 100
post_draws <- sample.int(dim(predS)[1], n_post_draws)
y_s_sim <- matrix(NA,n_post_draws,length(surv_data_list$y))
for(i in 1:n_post_draws){
  y_s_sim[i,] <- rbinom(n=length(surv_data_list$y), size=1, prob = invlogit(predS[post_draws[i],]))
}
ppc_dens_overlay(surv_data_list$y, y_s_sim)

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

#### seedling survival ppc ####
surv_seed_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_seedling_surv.rds")
predS <- rstan::extract(surv_seed_fit, pars = c("p"))$p
n_post_draws <- 100
post_draws <- sample.int(dim(predS)[1], n_post_draws)
y_s_sim <- matrix(NA,n_post_draws,length(seed_surv_data_list$y))
for(i in 1:n_post_draws){
  y_s_sim[i,] <- rbinom(n=length(seed_surv_data_list$y), size=1, prob = invlogit(predS[post_draws[i],]))
}
ppc_dens_overlay(seed_surv_data_list$y, y_s_sim)

mean_s_plot <-   ppc_stat(grow_data_list$y, y_s_sim, stat = "mean")
sd_s_plot <- ppc_stat(grow_data_list$y, y_s_sim, stat = "sd")
skew_s_plot <- ppc_stat(grow_data_list$y, y_s_sim, stat = "skewness")
kurt_s_plot <- ppc_stat(grow_data_list$y, y_s_sim, stat = "Lkurtosis")
grid.arrange(mean_s_plot,sd_s_plot,skew_s_plot,kurt_s_plot,  top = "Seedling Survival")


#### flowering ppc ####
flow_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_flw.rds")
predF <- rstan::extract(flow_fit, pars = c("p"))$p
n_post_draws <- 100
post_draws <- sample.int(dim(predF)[1], n_post_draws)
y_f_sim <- matrix(NA,n_post_draws,length(flw_data_list$y))
for(i in 1:n_post_draws){
  y_f_sim[i,] <- rbinom(n=length(flw_data_list$y), size=1, prob = invlogit(predF[post_draws[i],]))
}
ppc_dens_overlay(flw_data_list$y, y_f_sim)

mean_f_plot <-   ppc_stat(flw_data_list$y, y_f_sim, stat = "mean")
sd_f_plot <- ppc_stat(flw_data_list$y, y_f_sim, stat = "sd")
skew_f_plot <- ppc_stat(flw_data_list$y, y_f_sim, stat = "skewness")
kurt_f_plot <- ppc_stat(flw_data_list$y, y_f_sim, stat = "Lkurtosis")
grid.arrange(mean_f_plot,sd_f_plot,skew_f_plot,kurt_f_plot,  top = "Flowering")

# now we want to look at how the the model is fitting across sizes
flw_size_ppc <- size_moments_ppc(data = LTREB_data_forflw,
                                  y_name = "FLW_STAT_T",
                                  sim = y_f_sim, 
                                  n_bins = 3, 
                                  title = "Flowering")


#### growth ppc ####
grow_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_grow.rds")
# pairs(grow_fit, pars = c("beta0"))
grow_par <- rstan::extract(grow_fit, pars = c("lambda","beta0","betasize","betaendo","betaorigin","tau_year","tau_plot", "phi", "od"))
predG <- grow_par$lambda
phiG <- grow_par$phi
odG <- grow_par$od

n_post_draws <- 100
post_draws <- sample.int(dim(predG)[1], n_post_draws)
y_g_sim <- matrix(NA,n_post_draws,length(grow_data_list$y))
for(i in 1:n_post_draws){
  for(j in 1:length(grow_data_list$y)){
    y_g_sim[i,j] <- sample(x = 1:max(grow_data_list$y), size = 1, replace = T, prob = dnbinom(1:max(grow_data_list$y), mu = exp(predG[post_draws[i],j]), size = odG[post_draws[i],j])/(1-dnbinom(0, mu = exp(predG[post_draws[i],j]), size = odG[post_draws[i],j])))
  }
}
ppc_dens_overlay(grow_data_list$y, y_g_sim)
ppc_dens_overlay(grow_data_list$y, y_g_sim) + xlim(0,60) + ggtitle("growth w/o seedling")
ppc_dens_overlay(grow_data_list$y, y_g_sim) + xlim(50,120)

mean_g_plot <-   ppc_stat(grow_data_list$y, y_g_sim, stat = "mean")
sd_g_plot <- ppc_stat(grow_data_list$y, y_g_sim, stat = "sd")
skew_g_plot <- ppc_stat(grow_data_list$y, y_g_sim, stat = "skewness")
kurt_g_plot <- ppc_stat(grow_data_list$y, y_g_sim, stat = "Lkurtosis")
grid.arrange(mean_g_plot,sd_g_plot,skew_g_plot,kurt_g_plot,  top = "Growth ZTNB")


# now we want to look at how the the growth model is fitting especially across sizes, as this determines the transitions through the matrix model
# data bins

ZTNB_growth_size_ppc <- size_moments_ppc(data = LTREB_data_forgrow,
                                    y_name = "size_t1",
                                    sim = y_g_sim, 
                                    n_bins = 6, 
                                    title = "Growth ZTNB")

#growth PIG distribution
grow_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_grow_PIG.rds")
# pairs(grow_fit, pars = c("beta0"))
grow_par <- rstan::extract(grow_fit, pars = c("predG","beta0","betasize","betaendo","betaorigin","tau_year","tau_plot","theta", "sigma"))
predG <- grow_par$predG
theta <- grow_par$theta
n_post_draws <- 100
post_draws <- sample.int(dim(predG)[1], n_post_draws)
y_g_sim   <-  matrix(NA,n_post_draws,length(grow_data_list$y))

# simulate data
for(i in 1:n_post_draws){
  ## sample growth data (zero-truncated PIG)
  for(j in 1:length(grow_data_list$y)){
    # probability without truncation
    prob_v <- dpois( 1:1000, 
                     lambda = (predG[i,j] * theta[i,j]) )
    # probability for trunctation (denominator)
    prob_t <- (1 - dpois(0, lambda = (predG[i,j] * theta[i,j]) ) )
    
    y_g_sim[i,j] <- sample( x=1:1000, 
                            size = 1, replace = T,
                            prob = prob_v / prob_t )
  }
}

# Posterior predictive check
ppc_dens_overlay(grow_data_list$y, y_g_sim)
ppc_dens_overlay(grow_data_list$y, y_g_sim) + xlim(0,60) + ggtitle("growth w/o seedling PIG")
ppc_dens_overlay(grow_data_list$y, y_g_sim) + xlim(50,120)

mean_g_plot <-   ppc_stat(grow_data_list$y, y_g_sim, stat = "mean")
sd_g_plot <- ppc_stat(grow_data_list$y, y_g_sim, stat = "sd")
skew_g_plot <- ppc_stat(grow_data_list$y, y_g_sim, stat = "skewness")
kurt_g_plot <- ppc_stat(grow_data_list$y, y_g_sim, stat = "Lkurtosis")
grid.arrange(mean_g_plot,sd_g_plot,skew_g_plot,kurt_g_plot, top = "Growth PIG")


# now we want to look at how the the growth model is fitting especially across sizes, as this determines the transitions through the matrix model
PIG_growth_size_ppc <- size_moments_ppc(data = LTREB_data_forgrow,
                                         y_name = "size_t1",
                                         sim = y_g_sim, 
                                         n_bins = 6, 
                                         title = "Growth PIG")

#### seedling growth ppc ####

s_grow_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_seedling_grow.rds")
# pairs(grow_fit, pars = c("beta0"))
s_grow_par <- rstan::extract(s_grow_fit, pars = c("lambda","beta0","betaendo","tau_year","tau_plot", "phi", "od"))
predG <- s_grow_par$lambda
phiG <- s_grow_par$phi
odG <- s_grow_par$od

n_post_draws <- 100
post_draws <- sample.int(dim(predG)[1], n_post_draws)
y_g_sim <- matrix(NA,n_post_draws,length(seed_grow_data_list$y))
for(i in 1:n_post_draws){
  for(j in 1:length(seed_grow_data_list$y)){
    y_g_sim[i,j] <- sample(x = 1:max(seed_grow_data_list$y), size = 1, replace = T, prob = dnbinom(1:max(seed_grow_data_list$y), mu = exp(predG[post_draws[i],j]), size = odG[post_draws[i],j])/(1-dnbinom(0, mu = exp(predG[post_draws[i],j]), size = odG[post_draws[i],j])))
  }
}
ppc_dens_overlay(seed_grow_data_list$y, y_g_sim)
ppc_dens_overlay(seed_grow_data_list$y, y_g_sim) + xlim(0,20) + ggtitle("seedling growth")

mean_seed_g_plot <- ppc_stat(seed_grow_data_list$y, y_g_sim, stat = "mean") + ggtitle("mean")
sd_seed_g_plot <- ppc_stat(seed_grow_data_list$y, y_g_sim, stat = "sd")+ ggtitle("sd")
skew_seed_g_plot <- ppc_stat(seed_grow_data_list$y, y_g_sim, stat = "skewness")+ ggtitle("skew")
kurt_seed_g_plot <- ppc_stat(seed_grow_data_list$y, y_g_sim, stat = "Lkurtosis")+ ggtitle("kurt")
grid.arrange(mean_seed_g_plot,sd_seed_g_plot,skew_seed_g_plot,kurt_seed_g_plot, top = "Seedling Growth")

#### fertility ppc ####
# This fit is maybe not super great, but it's good enough for now probably.
fert_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_fert_noplot.rds")
fert_par <- rstan::extract(fert_fit, pars = c("lambda","beta0","betasize","betaendo","betaorigin","tau_year", "phi", "od"))
predFert <- fert_par$lambda
phiFert <- fert_par$phi
odFert <- fert_par$od

n_post_draws <- 100
post_draws <- sample.int(dim(predFert)[1], n_post_draws)
y_fert_sim <- matrix(NA,n_post_draws,length(fert_data_list$y))
for(i in 1:n_post_draws){
  for(j in 1:length(fert_data_list$y)){
    y_fert_sim[i,j] <- sample(x = 1:max(fert_data_list$y), size = 1, replace = T, prob = dnbinom(1:max(fert_data_list$y), mu = exp(predFert[post_draws[i],j]), size = odFert[post_draws[i],j])/(1-dnbinom(0, mu = exp(predFert[post_draws[i],j]), size = odFert[post_draws[i],j])))
  }
}
ppc_dens_overlay(fert_data_list$y, y_fert_sim)
ppc_dens_overlay(fert_data_list$y, y_fert_sim) + xlim(0,20)
ppc_dens_overlay(fert_data_list$y, y_fert_sim) + xlim(20,80)

mean_fert_plot <-   ppc_stat(fert_data_list$y, y_fert_sim, stat = "mean")
sd_fert_plot <- ppc_stat(fert_data_list$y, y_fert_sim, stat = "sd")
skew_fert_plot <- ppc_stat(fert_data_list$y, y_fert_sim, stat = "skewness")
kurt_fert_plot <- ppc_stat(fert_data_list$y, y_fert_sim, stat = "Lkurtosis")
grid.arrange(mean_fert_plot,sd_fert_plot,skew_fert_plot,kurt_fert_plot)

# Now we can look at the binned size fit for fertility
fert_size_ppc <- size_moments_ppc(data = LTREB_data_forfert,
                                        y_name = "FLW_COUNT_T",
                                        sim = y_fert_sim, 
                                        n_bins = 6, 
                                        title = "Inflorescence Count")

mcmc_trace(fert_fit,pars=c("betaendo[1]","betaendo[2]","betaendo[3]"
                           ,"betaendo[4]","betaendo[5]","betaendo[6]"
                           ,"betaendo[7]"))


#### spikelet ppc ####
spike_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_spike.rds")
predSpike <- rstan::extract(spike_fit, pars = c("p"))$p
n_post_draws <- 100
post_draws <- sample.int(dim(predSpike)[1], n_post_draws)
y_spike_sim <- matrix(NA,n_post_draws,length(spike_data_list$y))
for(i in 1:n_post_draws){
  y_spike_sim[i,] <- rpois(n=length(spike_data_list$y), lambda = exp(predSpike[post_draws[i],]))
}
ppc_dens_overlay(spike_data_list$y, y_spike_sim)

mcmc_trace(spike_fit,pars=c("betaendo[1]","betaendo[2]","betaendo[3]"
                           ,"betaendo[4]","betaendo[5]","betaendo[6]"
                           ,"betaendo[7]"))
