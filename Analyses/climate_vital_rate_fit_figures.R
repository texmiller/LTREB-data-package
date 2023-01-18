## Title: Grass endophyte vital rate fit visualizations and plots for all climate_explicite vital rates and all species
## Authors: Joshua and Tom
#############################################################


library(tidyverse)
library(RColorBrewer)
library(rstan)
library(StanHeaders)
library(bayesplot)
library(devtools)
library(moments)
library(patchwork)
library(actuar) # for PIG distribution 


invlogit<-function(x){exp(x)/(1+exp(x))}
logit = function(x) { log(x/(1-x)) }
Lkurtosis=function(x) log(kurtosis(x)); 

quote_bare <- function( ... ){
  substitute( alist(...) ) %>% 
    eval( ) %>% 
    sapply( deparse )
}
#############################################################################################
####### Reading in the demographic data ------------------
#############################################################################################

#  data are prepared in the endodemog_data_processing.R file, 
# source("Analyses/endodemog_data_processing.R")
species_factor_key <- c("AGPE" = 1, "ELRI" = 2, "ELVI" = 3, 
                        "FESU" = 4, "LOAR" = 5, "POAL" = 6, 
                        "POSY" = 7)
LTREB_full <- read_csv("~/Dropbox/EndodemogData/Fulldataplusmetadata/LTREB_full.csv")
LTREB_repro1 <- read_csv("~/Dropbox/EndodemogData/Fulldataplusmetadata/LTREB_repro1.csv")

# LTREB_full <- read_csv("~/Dropbox/EndodemogData/Fulldataplusmetadata/LTREB_full.csv")
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
  filter(!is.na(FLW_STAT_T1)) %>% 
  filter(!is.na(logsize_t1)) %>% 
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
  filter(!is.na(FLW_COUNT_T1)) %>%
  filter(FLW_COUNT_T1 > 0) %>%
  filter(!is.na(logsize_t1))
dim(LTREB_data_forfert)

LTREB_data_forspike <- LTREB_full %>%
  dplyr::select(-FLW_COUNT_T, -FLW_STAT_T, -SPIKE_A_T, -SPIKE_B_T, -SPIKE_C_T, -SPIKE_D_T, -SPIKE_AGPE_MEAN_T, -census_month, -year, -spei1,-spei24, -annual_temp, -annual_precip, -endo_status_from_check, -plot_endo_for_check, -endo_mismatch, -dist_a, -dist_b) %>% 
  filter(!is.na(FLW_STAT_T1)) %>% 
  filter(FLW_STAT_T1>0) %>% 
  reshape2::melt(id.var = c("plot_fixed" ,   "plot_index",         "pos"         ,           "id",
                  "species"       ,         "species_index"  ,        "endo_01",
                  "endo_index"  ,           "origin_01"       ,       "birth" , "spei12",
                  "year_t1"         ,       "year_t1_index"       ,   "surv_t1" ,
                  "size_t1"         ,       "logsize_t1"       ,
                  "year_t",
                  "year_t_index"     ,      "size_t"           ,      "logsize_t"  ,
                  "FLW_COUNT_T1"      ,      "FLW_STAT_T1"),
       value.name = "spike_count_t1") %>% 
  rename(spikelet_id = variable) %>% 
  filter(!is.na(spike_count_t1), spike_count_t1 > 0) %>% 
  mutate(spike_count_t1 = as.integer(spike_count_t1))


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
seed_grow_data_list <- list(y = as.integer(LTREB_grow_seedling$size_t1),
                            spei = as.numeric(LTREB_grow_seedling$spei12),
                            spei_nl = as.numeric(LTREB_grow_seedling$spei12^2),
                            logsize = LTREB_grow_seedling$logsize_t,
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

fert_data_list <- list(y = as.integer(LTREB_data_forfert$FLW_COUNT_T1),
                       spei = as.numeric(LTREB_data_forfert$spei12),
                       spei_nl = as.numeric(LTREB_data_forfert$spei12^2),
                       logsize = LTREB_data_forfert$logsize_t1,
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
                        spp = as.integer(LTREB_data_forspike$species_index),
                        y = LTREB_data_forspike$spike_count_t1,
                        logsize = LTREB_data_forspike$logsize_t1,
                        endo_01 = LTREB_data_forspike$endo_01,
                        origin_01 = LTREB_data_forspike$origin_01, 
                        spei = as.numeric(LTREB_data_forspike$spei12),
                        spei_nl = as.numeric(LTREB_data_forspike$spei12^2))
str(spike_data_list)

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

s_to_s_data_list <- read_rds("Analyses/s_to_s_data_list.rds")



#############################################################################################
####### Read in matrix population functions ------------------
#############################################################################################
# we can use functions from this to set up our parameters for the plots with mean and variance effects
source("Analyses/MPM_functions.R")



#############################################################################################
####### Read in Stan vital rate model outputs ------------------
#############################################################################################

# The stan objects for each vital rate
# spei_ includes only linear spei term, while spei2_ includes quadratice spei term.
spei_surv_fit_seedling <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_seedling_surv_withoutquadraticterm.rds")
spei2_surv_fit_seedling <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_seedling_surv_withtightpriors.rds")

spei_surv_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_spp_surv_woseedling_linear.rds")
spei2_surv_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_spp_surv_woseedling.rds")


spei_grow_fit_seedling <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_seedling_grow_linear_10000iterations.rds")
spei2_grow_fit_seedling <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_seedling_grow_nonlinear_10000iterations.rds") # not run yet


spei_grow_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_grow_PIG_linear.rds")
spei2_grow_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_grow_PIG_nonlinear.rds") 


spei_flw_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_spp_flw.rds")
spei2_flw_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_spp_flw_nonlinear.rds")


spei_fert_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_fert_PIG_linear.rds")
spei2_fert_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_fert_PIG_nonlinear.rds")


spei_spike_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_spp_spike_linear.rds")
spei2_spike_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_spp_spike_nonlinear.rds")



seedmean_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_seed_mean.rds") # doesn't include spei predictor

spei_stos_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_spp_s_to_s.rds") #only fit with linear term



# Pulling out the actual parameters
#survival
spei_surv_par <- rstan::extract(spei_surv_fit, pars =quote_bare(beta0,betasize,betaendo,betaorigin,
                                                                betaspei_endo,
                                                      tau_year, tau_plot, sigma_year))
spei2_surv_par <- rstan::extract(spei2_surv_fit, pars =quote_bare(beta0,betasize,betaendo,betaorigin,
                                                                  betaspei_endo,betaspei_nl_endo,
                                                      tau_year, tau_plot, sigma_year))
# seedling survival
spei_surv_seedling_par <- rstan::extract(spei_surv_fit_seedling, pars =quote_bare(beta0,betaendo,
                                                                                  betaspei_endo,
                                                                tau_year, tau_plot, sigma_year))
spei2_surv_seedling_par <- rstan::extract(spei2_surv_fit_seedling, pars =quote_bare(beta0,betaendo,
                                                                                    betaspei_endo, betaspei_nl_endo,
                                                                  tau_year, tau_plot, sigma_year))
# growth
spei_grow_par <- rstan::extract(spei_grow_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin,
                                                                 betaspei_endo,
                                                       tau_year, tau_plot,
                                                       sigma, sigma_year))
spei2_grow_par <- rstan::extract(spei2_grow_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin,
                                                                 betaspei_endo, betaspei_nl_endo,
                                                                 tau_year, tau_plot,
                                                                 sigma, sigma_year))

# seedling growth
spei_grow_sdlg_par <- rstan::extract(spei_grow_fit_seedling, pars = quote_bare(beta0,betaendo,
                                                                               betaspei_endo,
                                                                     tau_year, tau_plot,
                                                                     sigma, sigma_year))
spei2_grow_sdlg_par <- rstan::extract(spei2_grow_fit_seedling, pars = quote_bare(beta0,betaendo,
                                                                               betaspei_endo,betaspei_nl_endo,
                                                                               tau_year, tau_plot,
                                                                               sigma, sigma_year))

#flowering
spei_flow_par <- rstan::extract(spei_flw_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin,
                                                      betaspei_endo,
                                                      tau_year, tau_plot, sigma_year))
spei2_flow_par <- rstan::extract(spei2_flw_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin,
                                                           betaspei_endo,  betaspei_nl_endo,
                                                           tau_year, tau_plot, sigma_year))
# number of flower tillers
spei_fert_par <- rstan::extract(spei_fert_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin,
                                                                betaspei_endo,
                                                                tau_year, tau_plot, sigma_year))
spei2_fert_par <- rstan::extract(spei2_fert_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin,
                                                                  betaspei_endo,  betaspei_nl_endo,
                                                                  tau_year, tau_plot, sigma_year))
# spikelets/infl
spei_spike_par <- rstan::extract(spei_spike_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin,
                                                                 betaspei_endo,
                                                                 tau_year, tau_plot, sigma_year))
spei2_spike_par <- rstan::extract(spei2_spike_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin,
                                                                   betaspei_endo,  betaspei_nl_endo,
                                                                   tau_year, tau_plot, sigma_year))

seed_par <- rstan::extract(seedmean_fit, pars = quote_bare(beta0,betaendo)) #no plot or year effect

# seedling recruitment
spei_recruit_par <- rstan::extract(spei_stos_fit, pars = quote_bare(beta0,betaendo,
                                                                    betaspei_endo,
                                                                    tau_year, tau_plot, sigma_year))

# Saved y_rep values from the script "vital_rate_analysis.R", "seed_means.R" and "seed_to_seedling.R" 
# for the linear models
spei_y_s_sim <- readRDS(file = "yrep_climatesurvivalmodel_linear.rds")
spei_y_seed_s_sim <- readRDS(file = "yrep_climate_seedlingsurvivalmodel_linear.rds")
spei_y_f_sim <- readRDS(file = "yrep_climate_floweringmodel_linear.rds")
spei_y_g_sim <- readRDS(file = "yrep_climate_growthPIGmodel_linear.rds")
spei_y_seed_g_sim <- readRDS(file = "yrep_climate_seedlinggrowthPIGmodel_linear.rds")
spei_y_fert_sim <- readRDS(file = "yrep_climate_fertilityPIGmodel_linear.rds")
spei_y_spike_sim <- readRDS(file = "yrep_climate_spikeletNBmodel_linear.rds")

y_seedmean_sim <- readRDS(file = "yrep_seedmeanmodel.rds")
spei_y_recruit_sim <- readRDS(file = "yrep_climate_stosmodel_linear.rds")



#########################################################################################################
# Plots for all species all vital rates and model fits ------------------------------
#########################################################################################################
# Set color scheme based on analine blue
endophyte_color_scheme <- c("#fdedd3","#f3c8a8", "#5a727b", "#4986c7", "#181914",  "#163381")
color_scheme_set(endophyte_color_scheme)
# color_scheme_view()

# And creating a color palette for each year
yearcount = length(unique(LTREB_full$year_t))
yearcolors<- colorRampPalette(brewer.pal(8,"Dark2"))(yearcount)
# scales::show_col(yearcolors)
species_list <- c("AGPE", "ELRI", "ELVI", "FESU", "LOAR", "POAL", "POSY")


## Plot for all models vital rate fits #####
# survival
surv_densplot <- ppc_dens_overlay(surv_data_list$y, spei_y_s_sim) + theme_classic() + labs(title = "Adult Survival", x = "Survival status", y = "Density")
# surv_densplot
# ggsave(surv_densplot, filename = "surv_densplot.png", width = 4, height = 4)

mean_s_plot <-   ppc_stat(surv_data_list$y, spei_y_s_sim, stat = "mean")
sd_s_plot <- ppc_stat(surv_data_list$y, spei_y_s_sim, stat = "sd")
skew_s_plot <- ppc_stat(surv_data_list$y, spei_y_s_sim, stat = "skewness")
kurt_s_plot <- ppc_stat(surv_data_list$y, spei_y_s_sim, stat = "Lkurtosis")
surv_moments <- mean_s_plot+sd_s_plot+skew_s_plot+kurt_s_plot + plot_annotation(title = "Survival")
# surv_moments
# ggsave(surv_moments, filename = "surv_momentplot.png", width = 4, height = 4)

#seedliing survival
seedsurv_densplot <- ppc_dens_overlay(seed_surv_data_list$y, spei_y_seed_s_sim) + theme_classic() + labs(title = "Seedling Survival", x = "Survival status", y = "Density")
# seedsurv_densplot
# ggsave(seedsurv_densplot, filename = "seedsurv_densplot.png", width = 4, height = 4)

mean_s_plot <-   ppc_stat(seed_surv_data_list$y, spei_y_seed_s_sim, stat = "mean")
sd_s_plot <- ppc_stat(seed_surv_data_list$y, spei_y_seed_s_sim, stat = "sd")
skew_s_plot <- ppc_stat(seed_surv_data_list$y, spei_y_seed_s_sim, stat = "skewness")
kurt_s_plot <- ppc_stat(seed_surv_data_list$y, spei_y_seed_s_sim, stat = "Lkurtosis")
seedsurv_moments <- mean_s_plot+sd_s_plot+skew_s_plot+kurt_s_plot + plot_annotation(title = "Seedling Survival")
# seedsurv_moments
# ggsave(seedsurv_moments, filename = "seedsurv_momentsplot.png", width = 4, height = 4)

# Flowering
flw_densplot <- ppc_dens_overlay(flw_data_list$y, spei_y_f_sim) + theme_classic() + labs(title = "Flowering", x = "Flowering status", y = "Density")
# flw_densplot
# ggsave(flw_densplot, filename = "flw_densplot.png", width = 4, height = 4)

mean_f_plot <-   ppc_stat(flw_data_list$y, spei_y_f_sim, stat = "mean")
sd_f_plot <- ppc_stat(flw_data_list$y, spei_y_f_sim, stat = "sd")
skew_f_plot <- ppc_stat(flw_data_list$y, spei_y_f_sim, stat = "skewness")
kurt_f_plot <- ppc_stat(flw_data_list$y, spei_y_f_sim, stat = "Lkurtosis")
flw_moments <- mean_f_plot+sd_f_plot+skew_f_plot+kurt_f_plot +plot_annotation(title = "Flowering")
# flw_moments
# ggsave(flw_moments, filename = "flw_momentsplot.png", width = 4, height = 4)

# Growth (with the poisson inverse gaussian distribution)
grow_densplot <- ppc_dens_overlay(grow_data_list$y, spei_y_g_sim) + xlim(0,60) + theme_classic() + labs(title = "Adult Growth with PIG", x = "No. of Tillers", y = "Density")
# grow_densplot
# ggsave(grow_densplot, filename = "grow_densplot.png", width = 4, height = 4)

mean_g_plot <-   ppc_stat(grow_data_list$y, spei_y_g_sim, stat = "mean")
sd_g_plot <- ppc_stat(grow_data_list$y, spei_y_g_sim, stat = "sd")
skew_g_plot <- ppc_stat(grow_data_list$y, spei_y_g_sim, stat = "skewness")
kurt_g_plot <- ppc_stat(grow_data_list$y, spei_y_g_sim, stat = "Lkurtosis")
grow_moments <- mean_g_plot+sd_g_plot+skew_g_plot+kurt_g_plot+ plot_annotation(title = "Growth PIG")
# grow_moments
# ggsave(grow_moments, filename = "grow_momentsplot.png", width = 4, height = 4)

# seedling growth (with poisson inverse gaussian distribution)
seedgrow_densplot <- ppc_dens_overlay(seed_grow_data_list$y, spei_y_seed_g_sim) + xlim(0,60) + theme_classic() + labs(title = "Seedling Growth with PIG", x = "No. of Tillers", y = "Density")
# seedgrow_densplot
# ggsave(seedgrow_densplot, filename = "seed_grow_densplot.png", width = 4, height = 4)

mean_seed_g_plot <-   ppc_stat(seed_grow_data_list$y, spei_y_seed_g_sim, stat = "mean")
sd_seed_g_plot <- ppc_stat(seed_grow_data_list$y, spei_y_seed_g_sim, stat = "sd")
skew_seed_g_plot <- ppc_stat(seed_grow_data_list$y, spei_y_seed_g_sim, stat = "skewness")
kurt_seed_g_plot <- ppc_stat(seed_grow_data_list$y, spei_y_seed_g_sim, stat = "Lkurtosis")
seedgrow_moments <- mean_seed_g_plot+sd_seed_g_plot+skew_seed_g_plot+kurt_seed_g_plot+ plot_annotation(title = "Seedling Growth PIG")
# seedgrow_moments
# ggsave(seedgrow_moments, filename = "seedgrow_momentsplot.png", width = 4, height = 4)

# Fertility
fert_densplot <- ppc_dens_overlay(fert_data_list$y, spei_y_fert_sim) + xlim(0,40) + ggtitle("Fertility with PIG")
# fert_densplot
# ggsave(fert_densplot, filename = "fert_densplot_withPIG.png", width = 4, height = 4)

# This doesn't fit all the moments quite as well, but it's likely good enough
mean_fert_plot <-   ppc_stat(fert_data_list$y, spei_y_fert_sim, stat = "mean")
sd_fert_plot <- ppc_stat(fert_data_list$y, spei_y_fert_sim, stat = "sd")
skew_fert_plot <- ppc_stat(fert_data_list$y, spei_y_fert_sim, stat = "skewness")
kurt_fert_plot <- ppc_stat(fert_data_list$y, spei_y_fert_sim, stat = "Lkurtosis")
fert_moments <- mean_fert_plot+sd_fert_plot+skew_fert_plot+kurt_fert_plot + plot_annotation(title = "Fertility with PIG")
# fert_moments
# ggsave(fert_moments, filename = "fert_momentsplot_withPIG.png", width = 4, height = 4)

# Spikelets per infl
spike_densplot <- ppc_dens_overlay(spike_data_list$y, spei_y_spike_sim) + xlim(0,250) + ggtitle("Spikelet Count")
# spike_densplot
# ggsave(spike_densplot, filename = "spike_densplot.png", width = 4, height = 4)

mean_spike_plot <-   ppc_stat(spike_data_list$y, spei_y_spike_sim, stat = "mean")
sd_spike_plot <- ppc_stat(spike_data_list$y, spei_y_spike_sim, stat = "sd")
skew_spike_plot <- ppc_stat(spike_data_list$y, spei_y_spike_sim, stat = "skewness")
kurt_spike_plot <- ppc_stat(spike_data_list$y, spei_y_spike_sim, stat = "Lkurtosis")
spike_moments <- mean_spike_plot+sd_spike_plot+skew_spike_plot+kurt_spike_plot+ plot_annotation(title = "Spikelets per Infl. ZTNB")
# spike_moments
# ggsave(spike_moments, filename = "spike_momentsplot.png", width = 4, height = 4)

# Seed mean per spikelet
seedmean_densplot <-ppc_dens_overlay(seed_mean_data_list$seed, y_seedmean_sim)+ theme_classic() + labs(title = "Seed Means", x = "Seeds per Spikelet.", y = "Density") 
# ggsave(seedmean_densplot, filename = "seedmean_densplot.png", width = 4, height = 4)

mean_sm_plot <-   ppc_stat(seed_mean_data_list$seed, y_seedmean_sim, stat = "mean")
sd_sm_plot <- ppc_stat(seed_mean_data_list$seed, y_seedmean_sim, stat = "sd")
skew_sm_plot <- ppc_stat(seed_mean_data_list$seed, y_seedmean_sim, stat = "skewness")
kurt_sm_plot <- ppc_stat(seed_mean_data_list$seed, y_seedmean_sim, stat = "Lkurtosis")
seedmean_moments <- mean_sm_plot + sd_sm_plot + skew_sm_plot + kurt_sm_plot
# seedmean_moments
# ggsave(seedmean_moments, filename = "seedmean_moments.png", width = 4, height = 4)

# Seed to Seedling
stos_densplot <- ppc_dens_overlay(s_to_s_data_list$tot_recruit_t1, spei_y_recruit_sim) +xlim(0,75) + labs(title = "Recruitment", x = "Successful Germination", y = "Density") 
# ggsave(stos_densplot, filename = "stos_densplot.png", width = 4, height = 4)

mean_stos_plot <-   ppc_stat(s_to_s_data_list$tot_recruit_t1, spei_y_recruit_sim, stat = "mean")
sd_stos_plot <- ppc_stat(s_to_s_data_list$tot_recruit_t1, spei_y_recruit_sim, stat = "sd")
skew_stos_plot <- ppc_stat(s_to_s_data_list$tot_recruit_t1, spei_y_recruit_sim, stat = "skewness")
kurt_stos_plot <- ppc_stat(s_to_s_data_list$tot_recruit_t1, spei_y_recruit_sim, stat = "Lkurtosis")
stos_moments <- mean_stos_plot+sd_stos_plot+skew_stos_plot+kurt_stos_plot
# stos_moments
# ggsave(stos_moments, filename = "stos_moments.png", width = 4, height = 4)


# all together
## Plot for all models fits and moments
fitsandmoments_plot <- (surv_densplot + surv_moments)/
  (seedsurv_densplot + seedsurv_moments)/
  (grow_densplot + grow_moments)/ 
  (seedgrow_densplot + seedgrow_moments)/
  (flw_densplot + flw_moments)/
  (fert_densplot + fert_moments)/
  (spike_densplot + spike_moments)/
  (seedmean_densplot + seedmean_moments)/
  (stos_densplot + stos_moments)+
  plot_annotation(title = "Vital rate fits and moments with 500 posterior draws",
                  subtitle = "~SPEI")
ggsave(fitsandmoments_plot, filename = "spei_fitsandmoments_plot.png", width = 18, height = 20)

## Plot for size-specific moments for growth model minus seedlings, which have only one size ####
## could plot these for other vital rates if desired

# Function for looking at binned size_t fits, particularly important for the growth kernel as this determines the transitions through the matrix model
# plots the mean, sd, skew and kertosis of the posteriors (grey) as well as the mean of the posteriors for each moment (black) and the data (red) for size bins
size_moments_ppc <- function(data,y_name,sim, n_bins, title = NA){
  require(tidyverse)
  require(patchwork)
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
  size_ppc_plot <- meanplot+ sdplot+skewplot+ kurtplot+plot_annotation(title = title)
  return(size_ppc_plot)
}

spei_PIG_growth_size_ppc <- size_moments_ppc(data = LTREB_data_forgrow,
                                        y_name = "size_t1",
                                        sim = spei_y_g_sim, 
                                        n_bins = 6, 
                                        title = "Growth PIG")
# ggsave(spei_PIG_growth_size_ppc, filename = "spei_PIG_growth_size_pcc.png", width = 4, height = 4)

spei_size_ppc_plot <- (spei_PIG_growth_size_ppc+ plot_annotation(title = 'Growth'))+
  plot_annotation(title = "Size specific vital rate moments",
                  subtitle = "~SPEI")
# size_ppc_plot
ggsave(spei_size_ppc_plot, filename = "spei_size_ppc_plot.png", width = 6, height = 8)


## Traceplots for select parameters from all vital rates for AGPE
spei_surv_trace <- mcmc_trace(spei_surv_fit, pars = c("beta0[1]", "betaendo[1]", "sigmaendo[1]"))+ggtitle("Survival")
spei_seedsurv_trace <- mcmc_trace(spei_surv_fit_seedling, pars = c("beta0[1]", "betaendo[1]", "sigmaendo[1]"))+ggtitle("Seedling Survival")
spei_flw_trace <- mcmc_trace(spei_flw_fit, pars = c("beta0[1]", "betaendo[1]", "sigmaendo[1]"))+ggtitle("Flowering")
spei_grow_trace <- mcmc_trace(spei_grow_fit, pars = c("beta0[1]", "betaendo[1]", "sigmaendo[1]"))+ggtitle("Growth")
spei_seedgrow_trace <- mcmc_trace(spei_grow_fit_seedling, pars = c("beta0[1]", "betaendo[1]", "sigmaendo[1]"))+ggtitle("Seedling Growth")
spei_fert_trace <- mcmc_trace(spei_fert_fit, pars = c("beta0[1]", "betaendo[1]", "sigmaendo[1]"))+ggtitle("Fertility")
spei_spike_trace <- mcmc_trace(spei_spike_fit, pars = c("beta0[1]", "betaendo[1]", "sigmaendo[1]"))+ggtitle("Spikelets per inflorescence")
seedmean_trace <- mcmc_trace(seedmean_fit, pars = c("beta0[1]", "betaendo[1]"))+ggtitle("Mean Seeds per spikelet")
spei_stos_trace <- mcmc_trace(spei_stos_fit, pars = c("beta0[1]", "betaendo[1]", "sigmaendo[1]"))+ggtitle("Recruitment")

spei_vr_traceplots <- (spei_surv_trace)/
  (spei_seedsurv_trace)/
  (spei_grow_trace)/
  (spei_seedgrow_trace)/
  (spei_flw_trace)/
  (spei_fert_trace)/
  (spei_spike_trace)/
  (seedmean_trace)/
  (spei_stos_trace) + 
  plot_annotation(title = "Traceplots for select parameters from all vital rates (AGPE)",
                  subtitle = "~SPEI")
ggsave(spei_vr_traceplots, filename = "spei_vr_traceplots.png", width = 25, height = 20)



## Plots for mean and year variance endophyte effects on all vital rates, all species ####
mean_species_key <- c("betaendo[1]" = "AGPE", "betaendo[2]" = "ELRI", "betaendo[3]" = "ELVI", "betaendo[4]" = "FESU", "betaendo[5]" = "LOAR", "betaendo[6]" = "POAL", "betaendo[7]" = "POSY")
sd_species_key <- c("sigmaendo[1]" = "AGPE", "sigmaendo[2]" = "ELRI", "sigmaendo[3]" = "ELVI", "sigmaendo[4]" = "FESU", "sigmaendo[5]" = "LOAR", "sigmaendo[6]" = "POAL", "sigmaendo[7]" = "POSY")
spei_species_key <- c("betaspei_endo[1,1]" = "E- AGPE", "betaspei_endo[1,2]" = "E+ AGPE", "betaspei_endo[2,1]" = "E- ELRI", "betaspei_endo[2,2]" = "E+ ELRI", "betaspei_endo[3,1]" = "E- ELVI", "betaspei_endo[3,2]" = "E+ ELVI", "betaspei_endo[4,1]" = "E- FESU", "betaspei_endo[4,2]" = "E+ FESU",  "betaspei_endo[5,1]" = "E- LOAR", "betaspei_endo[5,2]" = "E+ LOAR", "betaspei_endo[6,1]" = "E- POAL", "betaspei_endo[6,2]" = "E+ POAL", "betaspei_endo[7,1]" = "E- POSY", "betaspei_endo[7,2]" = "E+ POSY")

# effect of endophytte on mean (betaendo)
spei_surv_endomean_posteriors <-  mcmc_areas(spei_surv_fit, prob = 0.8, regex_pars = c("betaendo"))+labs(title = "Adult Survival", subtitle = "Endophyte effect on mean with 80% credible intervals") + scale_y_discrete(labels = mean_species_key)
spei_seedsurv_endomean_posteriors <-  mcmc_areas(spei_surv_fit_seedling, prob = 0.8, regex_pars = c("betaendo"))+labs(title = "Seedling Survival", subtitle = "Endophyte effect on mean with 80% credible intervals") + scale_y_discrete(labels = mean_species_key)
spei_grow_endomean_posteriors <- mcmc_areas(spei_grow_fit, prob = 0.8, regex_pars = c("betaendo"))+labs(title = "Adult Growth", subtitle = "Endophyte effect on mean with 80% credible intervals") + scale_y_discrete(labels = mean_species_key)
spei_seedgrow_endomean_posteriors <- mcmc_areas(spei_grow_fit_seedling, prob = 0.8, regex_pars = c("betaendo"))+labs(title = "Seedling Growth", subtitle = "Endophyte effect on mean with 80% credible intervals") + scale_y_discrete(labels = mean_species_key)
spei_flw_endomean_posteriors <- mcmc_areas(spei_flw_fit, prob = 0.8, regex_pars = c("betaendo"))+labs(title = "Flowering", subtitle = "Endophyte effect on mean with 80% credible intervals") + scale_y_discrete(labels = mean_species_key)
spei_fert_endomean_posteriors <- mcmc_areas(spei_fert_fit, prob = 0.8, regex_pars = c("betaendo"))+labs(title = "Fertility", subtitle = "Endophyte effect on mean with 80% credible intervals") + scale_y_discrete(labels = mean_species_key)
spei_spike_endomean_posteriors <- mcmc_areas(spei_spike_fit, prob = 0.8, regex_pars = c("betaendo"))+labs(title = "Spikelets per infl.", subtitle = "Endophyte effect on mean with 80% credible intervals") + scale_y_discrete(labels = mean_species_key)
seedmean_endomean_posteriors <- mcmc_areas(seedmean_fit, prob = 0.8, regex_pars = c("betaendo"))+labs(title = "Mean seeds per spikelet", subtitle = "Endophyte effect on mean with 80% credible intervals") + scale_y_discrete(labels = mean_species_key)
spei_stos_endomean_posteriors <- mcmc_areas(spei_stos_fit, prob = 0.8, regex_pars = c("betaendo"))+labs(title = "Germination", subtitle = "Endophyte effect on mean with 80% credible intervals") + scale_y_discrete(labels = mean_species_key)

# effect of endophyte on year standard deviation (seed mean does not have an endophyte effect on variance)
spei_surv_endosd_posteriors <-  mcmc_areas(spei_surv_fit, prob = 0.8, regex_pars = c("sigmaendo"))+labs(title = "Adult Survival", subtitle = "Endophyte effect on SD with 80% credible intervals") + scale_y_discrete(labels = sd_species_key)
spei_seedsurv_endosd_posteriors <-  mcmc_areas(spei_surv_fit_seedling, prob = 0.8, regex_pars = c("sigmaendo"))+labs(title = "Seedling Survival", subtitle = "Endophyte effect on SD with 80% credible intervals") + scale_y_discrete(labels = sd_species_key)
spei_grow_endosd_posteriors <- mcmc_areas(spei_grow_fit, prob = 0.8, regex_pars = c("sigmaendo"))+labs(title = "Adult Growth", subtitle = "Endophyte effect on SD with 80% credible intervals") + scale_y_discrete(labels = sd_species_key)
spei_seedgrow_endosd_posteriors <- mcmc_areas(spei_grow_fit_seedling, prob = 0.8, regex_pars = c("sigmaendo"))+labs(title = "Seedling Growth", subtitle = "Endophyte effect on SD with 80% credible intervals") + scale_y_discrete(labels = sd_species_key)
spei_flw_endosd_posteriors <- mcmc_areas(spei_flw_fit, prob = 0.8, regex_pars = c("sigmaendo"))+labs(title = "Flowering", subtitle = "Endophyte effect on SD with 80% credible intervals") + scale_y_discrete(labels = sd_species_key)
spei_fert_endosd_posteriors <- mcmc_areas(spei_fert_fit, prob = 0.8, regex_pars = c("sigmaendo"))+labs(title = "Fertility", subtitle = "Endophyte effect on SD with 80% credible intervals") + scale_y_discrete(labels = sd_species_key)
spei_spike_endosd_posteriors <- mcmc_areas(spei_spike_fit, prob = 0.8, regex_pars = c("sigmaendo"))+labs(title = "Spikelets per infl.", subtitle = "Endophyte effect on SD with 80% credible intervals") + scale_y_discrete(labels = sd_species_key)
spei_stos_endosd_posteriors <- mcmc_areas(spei_stos_fit, prob = 0.8, regex_pars = c("sigmaendo"))+labs(title = "Germination", subtitle = "Endophyte effect on SD with 80% credible intervals") + scale_y_discrete(labels = sd_species_key)

# spei effect
spei_surv_betaspei_posteriors <-  mcmc_areas(spei_surv_fit, prob = 0.8, regex_pars = c("betaspei_endo"))+labs(title = "Adult Survival", subtitle = "Effect of SPEI by Endophyte Status with 80% credible interval") + scale_y_discrete(labels = spei_species_key)
spei_seedsurv_betaspei_posteriors <-  mcmc_areas(spei_surv_fit_seedling, prob = 0.8, regex_pars = c("betaspei_endo"))+labs(title = "Seedling Survival", subtitle = "Effect of SPEI by Endophyte Status with 80% credible interval") + scale_y_discrete(labels = spei_species_key)
spei_grow_betaspei_posteriors <- mcmc_areas(spei_grow_fit, prob = 0.8, regex_pars = c("betaspei_endo"))+labs(title = "Adult Growth", subtitle = "Effect of SPEI by Endophyte Status with 80% credible interval") + scale_y_discrete(labels = spei_species_key)
spei_seedgrow_betaspei_posteriors <- mcmc_areas(spei_grow_fit_seedling, prob = 0.8, regex_pars = c("betaspei_endo"))+labs(title = "Seedling Growth", subtitle = "Effect of SPEI by Endophyte Status with 80% credible interval") + scale_y_discrete(labels = spei_species_key)
spei_flw_betaspei_posteriors <- mcmc_areas(spei_flw_fit, prob = 0.8, regex_pars = c("betaspei_endo"))+labs(title = "Flowering", subtitle = "Effect of SPEI by Endophyte Status with 80% credible interval") + scale_y_discrete(labels = spei_species_key)
spei_fert_betaspei_posteriors <- mcmc_areas(spei_fert_fit, prob = 0.8, regex_pars = c("betaspei_endo"))+labs(title = "Fertility", subtitle = "Effect of SPEI by Endophyte Status with 80% credible interval") + scale_y_discrete(labels = spei_species_key)
spei_spike_betaspei_posteriors <- mcmc_areas(spei_spike_fit, prob = 0.8, regex_pars = c("betaspei_endo"))+labs(title = "Spikelets per infl.", subtitle = "Effect of SPEI by Endophyte Status with 80% credible interval") + scale_y_discrete(labels = spei_species_key)
spei_stos_betaspei_posteriors <- mcmc_areas(spei_stos_fit, prob = 0.8, regex_pars = c("betaspei_endo"))+labs(title = "Germination", subtitle = "Effect of SPEI by Endophyte Status with 80% credible interval") + scale_y_discrete(labels = spei_species_key)



spei_endomeanandvar_posteriors <- (spei_surv_endomean_posteriors + spei_surv_endosd_posteriors)/
  (spei_seedsurv_endomean_posteriors + spei_seedsurv_endosd_posteriors)/
  (spei_grow_endomean_posteriors + spei_grow_endosd_posteriors)/
  (spei_seedgrow_endomean_posteriors + spei_seedgrow_endosd_posteriors)/
  (spei_flw_endomean_posteriors + spei_flw_endosd_posteriors)/
  (spei_fert_endomean_posteriors + spei_fert_endosd_posteriors)/
  (spei_spike_endomean_posteriors + spei_spike_endosd_posteriors)/
  (seedmean_endomean_posteriors + plot_spacer())/
  (spei_stos_endomean_posteriors + spei_stos_endosd_posteriors) + plot_annotation(title = "Effect of SPEI by Endophyte Status",
                                                                                  subtitle = "~SPEI")

ggsave(spei_endomeanandvar_posteriors, filename = "spei_endomeanandvar_posteriors.png", width = 12, height = 30)


spei_betaspei_posteriors <- (spei_surv_betaspei_posteriors/
                               spei_seedsurv_betaspei_posteriors/
                               spei_grow_betaspei_posteriors/
                               spei_seedgrow_betaspei_posteriors/
                               spei_flw_betaspei_posteriors/
                               spei_fert_betaspei_posteriors/
                               spei_spike_betaspei_posteriors/
                               spei_stos_betaspei_posteriors )
ggsave(spei_betaspei_posteriors, filename = "spei_betaspei_posteriors.png", width = 6, height = 30)


## making plots with the histogram of E+ and E- variance for each vital rate ####
#surv
dimnames(spei_surv_par$sigma_year) <- list(Draw = paste0("i",1:dim(spei_surv_par$sigma_year)[1]), Species = species_list, Endo = c("E-", "E+"))
spei_surv_sigmayear_cube <- cubelyr::as.tbl_cube(spei_surv_par$sigma_year)
spei_surv_sigmayear_df <- as_tibble(spei_surv_sigmayear_cube)  %>% 
  rename(sigma_year = `spei_surv_par$sigma_year`)

spei_surv_sigmayear_hist <- ggplot(data = spei_surv_sigmayear_df)+
  geom_histogram(aes(x = sigma_year, fill = Endo), alpha = .9, bins = 250, position = "identity") +
  scale_fill_manual(values = c( endophyte_color_scheme[2], endophyte_color_scheme[4]))+
  facet_wrap(~Species, scales = "free", ncol = 1)+
  labs(title = "Survival", x = "Interannual SD")+
  theme_minimal()
# spei_surv_sigmayear_hist

#seedling surv
dimnames(spei_surv_seedling_par$sigma_year) <- list(Draw = paste0("i",1:dim(spei_surv_seedling_par$sigma_year)[1]), Species = species_list, Endo = c("E-", "E+"))
spei_seedsurv_sigmayear_cube <- cubelyr::as.tbl_cube(spei_surv_seedling_par$sigma_year)
spei_seedsurv_sigmayear_df <- as_tibble(spei_seedsurv_sigmayear_cube)  %>% 
  rename(sigma_year = `spei_surv_seedling_par$sigma_year`)

spei_seedsurv_sigmayear_hist <- ggplot(data = spei_seedsurv_sigmayear_df)+
  geom_histogram(aes(x = sigma_year, fill = Endo), alpha = .9, bins = 250, position = "identity") +
  scale_fill_manual(values = c( endophyte_color_scheme[2], endophyte_color_scheme[4]))+
  facet_wrap(~Species, scales = "free", ncol = 1)+
  labs(title = "Seedling Survival", x = "Interannual SD")+
  theme_minimal()
# spei_seedsurv_sigmayear_hist

#grow
dimnames(spei_grow_par$sigma_year) <- list(Draw = paste0("i",1:dim(spei_grow_par$sigma_year)[1]), Species = species_list, Endo = c("E-", "E+"))
spei_grow_sigmayear_cube <- cubelyr::as.tbl_cube(spei_grow_par$sigma_year)
spei_grow_sigmayear_df <- as_tibble(spei_grow_sigmayear_cube)  %>% 
  rename(sigma_year = `spei_grow_par$sigma_year`)

spei_grow_sigmayear_hist <- ggplot(data = spei_grow_sigmayear_df)+
  geom_histogram(aes(x = sigma_year, fill = Endo), alpha = .9, bins = 250, position = "identity") +
  scale_fill_manual(values = c( endophyte_color_scheme[2], endophyte_color_scheme[4]))+
  facet_wrap(~Species, scales = "free", ncol = 1)+
  labs(title = "Growth", x = "Interannual SD")+
  theme_minimal()
# spei_grow_sigmayear_hist

#seedling grow
dimnames(spei_grow_sdlg_par$sigma_year) <- list(Draw = paste0("i",1:dim(spei_grow_sdlg_par$sigma_year)[1]), Species = species_list, Endo = c("E-", "E+"))
spei_seedgrow_sigmayear_cube <- cubelyr::as.tbl_cube(spei_grow_sdlg_par$sigma_year)
spei_seedgrow_sigmayear_df <- as_tibble(spei_seedgrow_sigmayear_cube)  %>% 
  rename(sigma_year = `spei_grow_sdlg_par$sigma_year`)

spei_seedgrow_sigmayear_hist <- ggplot(data = spei_seedgrow_sigmayear_df)+
  geom_histogram(aes(x = sigma_year, fill = Endo), alpha = .9, bins = 250, position = "identity") +
  scale_fill_manual(values = c( endophyte_color_scheme[2], endophyte_color_scheme[4]))+
  facet_wrap(~Species, scales = "free", ncol = 1)+
  labs(title = "Seeding Growth", x = "Interannual SD")+
  theme_minimal()
# spei_seedgrow_sigmayear_hist

#flw
dimnames(spei_flow_par$sigma_year) <- list(Draw = paste0("i",1:dim(spei_flow_par$sigma_year)[1]), Species = species_list, Endo = c("E-", "E+"))
spei_flow_sigmayear_cube <- cubelyr::as.tbl_cube(spei_flow_par$sigma_year)
spei_flow_sigmayear_df <- as_tibble(spei_flow_sigmayear_cube)  %>% 
  rename(sigma_year = `spei_flow_par$sigma_year`)

spei_flow_sigmayear_hist <- ggplot(data = spei_flow_sigmayear_df)+
  geom_histogram(aes(x = sigma_year, fill = Endo), alpha = .9, bins = 250, position = "identity") +
  scale_fill_manual(values = c( endophyte_color_scheme[2], endophyte_color_scheme[4]))+
  facet_wrap(~Species, scales = "free", ncol = 1)+
  labs(title = "Flowering", x = "Interannual SD")+
  theme_minimal()
# spei_flow_sigmayear_hist

#fert
dimnames(spei_fert_par$sigma_year) <- list(Draw = paste0("i",1:dim(spei_fert_par$sigma_year)[1]), Species = species_list, Endo = c("E-", "E+"))
spei_fert_sigmayear_cube <- cubelyr::as.tbl_cube(spei_fert_par$sigma_year)
spei_fert_sigmayear_df <- as_tibble(spei_fert_sigmayear_cube)  %>% 
  rename(sigma_year = `spei_fert_par$sigma_year`)

spei_fert_sigmayear_hist <- ggplot(data = spei_fert_sigmayear_df)+
  geom_histogram(aes(x = sigma_year, fill = Endo), alpha = .9, bins = 250, position = "identity") +
  scale_fill_manual(values = c( endophyte_color_scheme[2], endophyte_color_scheme[4]))+
  facet_wrap(~Species, scales = "free", ncol = 1)+
  labs(title = "# of Flw Tillers", x = "Interannual SD")+
  theme_minimal()
# spei_fert_sigmayear_hist

#spike
dimnames(spei_spike_par$sigma_year) <- list(Draw = paste0("i",1:dim(spei_spike_par$sigma_year)[1]), Species = species_list, Endo = c("E-", "E+"))
spei_spike_sigmayear_cube <- cubelyr::as.tbl_cube(spei_spike_par$sigma_year)
spei_spike_sigmayear_df <- as_tibble(spei_spike_sigmayear_cube)  %>% 
  rename(sigma_year = `spei_spike_par$sigma_year`)

spei_spike_sigmayear_hist <- ggplot(data = spei_spike_sigmayear_df)+
  geom_histogram(aes(x = sigma_year, fill = Endo), alpha = .9, bins = 250, position = "identity") +
  scale_fill_manual(values = c( endophyte_color_scheme[2], endophyte_color_scheme[4]))+
  facet_wrap(~Species, scales = "free", ncol = 1)+
  labs(title = "Spikelets", x = "Interannual SD")+
  theme_minimal()
# spei_spike_sigmayear_hist

#germination
dimnames(spei_recruit_par$sigma_year) <- list(Draw = paste0("i",1:dim(spei_recruit_par$sigma_year)[1]), Species = species_list, Endo = c("E-", "E+"))
spei_recruit_sigmayear_cube <- cubelyr::as.tbl_cube(spei_recruit_par$sigma_year)
spei_recruit_sigmayear_df <- as_tibble(spei_recruit_sigmayear_cube)  %>% 
  rename(sigma_year = `spei_recruit_par$sigma_year`)

spei_recruit_sigmayear_hist <- ggplot(data = spei_recruit_sigmayear_df)+
  geom_histogram(aes(x = sigma_year, fill = Endo), alpha = .9, bins = 250, position = "identity") +
  scale_fill_manual(values = c( endophyte_color_scheme[2], endophyte_color_scheme[4]))+
  facet_wrap(~Species, scales = "free", ncol = 1)+
  labs(title = "Recruitment", x = "Interannual SD")+
  theme_minimal()
# spei_recruit_sigmayear_hist

# posterior histograms of E+ and E- variance

spei_endo_sigmayear_histograms <- spei_surv_sigmayear_hist + spei_seedsurv_sigmayear_hist + spei_grow_sigmayear_hist + spei_seedgrow_sigmayear_hist + spei_flow_sigmayear_hist + spei_fert_sigmayear_hist + spei_spike_sigmayear_hist + spei_recruit_sigmayear_hist+
  plot_layout( nrow  = 1, guides = "collect")+
  plot_annotation(title = "Endophyte effect on interannual SD across vital rates by species",
                  subtitle = "~SPEI")

ggsave(spei_endo_sigmayear_histograms, filename = "spei_endo_sigmayear_histograms.png", width = 30, height = 15)


###############################################################
###### Plots of vital rate fits with data by spei #############
###############################################################

max_size <- LTREB_full %>% 
  dplyr::select(species,species_index, size_t) %>% 
  filter(!is.na(size_t)) %>% 
  group_by(species, species_index) %>% 
  summarise(actual_max_size = max(size_t),
            mean_size = mean(size_t,na.rm = T),
            max_size = quantile(size_t,probs=0.975))

reprod_max_size <- LTREB_data_forfert %>% 
  dplyr::select(species,species_index, size_t) %>% 
  filter(!is.na(size_t)) %>% 
  group_by(species, species_index) %>% 
  summarise(actual_max_size = max(size_t),
            mean_size = mean(size_t,na.rm = T),
            max_size = quantile(size_t,probs=0.975))

spei_minmax <- LTREB_full %>% 
  dplyr::select(species,species_index, spei12) %>% 
  filter(!is.na(spei12)) %>% 
  group_by(species, species_index) %>% 
  summarise(min_spei = min(spei12),
            max_spei = max(spei12))

n_post_draws <- 500
post_draws <- sample.int(7500, n_post_draws)

x_seq_length <- 50
x_seq <- array(dim= c(x_seq_length,7))
for(s in 1:7){
  x_seq[,s] <-  seq(from = 1, to = filter(max_size, species_index == s)$actual_max_size, length.out = x_seq_length)
}
x_mean <- max_size$mean_size # the mean size for each species
repro_x_mean <- reprod_max_size$mean_size # looking at this separately because they fertility model looks like it's underfitting

# stepsize of spei values
spei_steps <- 10
spei_range <- array(dim = c(spei_steps,7))
for(s in 1:7){
  spei_range[,s] <- seq(from = spei_minmax$min_spei[s], to = spei_minmax$max_spei[s], length.out = 10)
}

surv_iter <- grow_iter <- flw_iter <- fert_iter <- spike_iter <-array(dim = c(spei_steps,2,7, n_post_draws))


sx<-function(x,params,spei=0){
  invlogit(params$surv_int + params$surv_spei*spei + params$surv_slope*log(x))
}
gx <- function(x,params,spei=0){
  exp(params$grow_int + params$grow_spei*spei + params$grow_slope*log(x))
}
flwx <- function(x,params,spei=0){
  invlogit(params$flow_int + params$flow_spei*spei + params$flow_slope*log(x))
}
fertx <- function(x,params,spei=0){
  exp(params$fert_int + params$fert_spei*spei + params$fert_slope*log(x))
}

spikex <- function(x,params,spei=0){
  exp(params$spike_int + params$spike_spei*spei + params$spike_slope*log(x))
}


for(i in 1:length(post_draws)){
  for(e in 1:2){
    for(s in 1:7){
      
      surv_iter[,e,s,i] <- sx(make_params(species=s,
                                          endo_mean=(e-1),
                                          endo_var=(e-1),
                                          original = mean(LTREB_full$origin_01), # should be =1 to represent recruit
                                          draw=post_draws[i],
                                          max_size=max_size,
                                          rfx=F,
                                          spei = T,
                                          surv_par=spei_surv_par,
                                          surv_sdlg_par =spei_surv_seedling_par,
                                          grow_par=spei_grow_par,
                                          grow_sdlg_par =spei_grow_sdlg_par,
                                          flow_par=spei_flow_par,
                                          fert_par=spei_fert_par,
                                          spike_par=spei_spike_par,
                                          seed_par=seed_par,
                                          recruit_par=spei_recruit_par), spei =  spei_range[,s], x = x_mean[s])
      grow_iter[,e,s,i] <- gx(make_params(species=s,
                                          endo_mean=(e-1),
                                          endo_var=(e-1),
                                          original = mean(LTREB_full$origin_01), # should be =1 to represent recruit
                                          draw=post_draws[i],
                                          max_size=max_size,
                                          rfx=F,
                                          spei = T,
                                          surv_par=spei_surv_par,
                                          surv_sdlg_par =spei_surv_seedling_par,
                                          grow_par=spei_grow_par,
                                          grow_sdlg_par =spei_grow_sdlg_par,
                                          flow_par=spei_flow_par,
                                          fert_par=spei_fert_par,
                                          spike_par=spei_spike_par,
                                          seed_par=seed_par,
                                          recruit_par=spei_recruit_par), spei =  spei_range[,s], x = x_mean[s])
      flw_iter[,e,s,i] <- flwx(make_params(species=s,
                                           endo_mean=(e-1),
                                           endo_var=(e-1),
                                           original = mean(LTREB_full$origin_01), # should be =1 to represent recruit
                                           draw=post_draws[i],
                                           max_size=max_size,
                                           rfx=F,
                                           spei = T,
                                           surv_par=spei_surv_par,
                                           surv_sdlg_par =spei_surv_seedling_par,
                                           grow_par=spei_grow_par,
                                           grow_sdlg_par =spei_grow_sdlg_par,
                                           flow_par=spei_flow_par,
                                           fert_par=spei_fert_par,
                                           spike_par=spei_spike_par,
                                           seed_par=seed_par,
                                           recruit_par=spei_recruit_par), spei =  spei_range[,s], x = x_mean[s])
      fert_iter[,e,s,i] <- fertx(make_params(species=s,
                                             endo_mean=(e-1),
                                             endo_var=(e-1),
                                             original = mean(LTREB_full$origin_01), # should be =1 to represent recruit
                                             draw=post_draws[i],
                                             max_size=max_size,
                                             rfx=F,
                                             spei = T,
                                             surv_par=spei_surv_par,
                                             surv_sdlg_par =spei_surv_seedling_par,
                                             grow_par=spei_grow_par,
                                             grow_sdlg_par =spei_grow_sdlg_par,
                                             flow_par=spei_flow_par,
                                             fert_par=spei_fert_par,
                                             spike_par=spei_spike_par,
                                             seed_par=seed_par,
                                             recruit_par=spei_recruit_par), spei =  spei_range[,s], x = x_mean[s])
      spike_iter[,e,s,i] <- spikex(make_params(species=s,
                                               endo_mean=(e-1),
                                               endo_var=(e-1),
                                               original = mean(LTREB_full$origin_01), # should be =1 to represent recruit
                                               draw=post_draws[i],
                                               max_size=max_size,
                                               rfx=F,
                                               spei = T,
                                               surv_par=spei_surv_par,
                                               surv_sdlg_par =spei_surv_seedling_par,
                                               grow_par=spei_grow_par,
                                               grow_sdlg_par =spei_grow_sdlg_par,
                                               flow_par=spei_flow_par,
                                               fert_par=spei_fert_par,
                                               spike_par=spei_spike_par,
                                               seed_par=seed_par,
                                               recruit_par=spei_recruit_par), spei =  spei_range[,s], x = x_mean[s])
      
    }
  }
}

#Now I'm gonna make these into  tidy dataframes for plotting

#getting spei values
dimnames(spei_range) <- list(spei = paste0("spei",1:spei_steps),Species = species_list)
spei_range_cube <- cubelyr::as.tbl_cube(spei_range)
spei_range_df <- as_tibble(spei_range_cube) %>% 
  rename(spei_value = spei_range)

dimnames(surv_iter) <-dimnames(grow_iter) <-dimnames(flw_iter) <-dimnames(fert_iter) <- dimnames(spike_iter) <-  list(spei = paste0("spei", 1:spei_steps), Endo = c("Eminus","Eplus"), Species = species_list, Iteration = paste0("iter",1:(n_post_draws)))
surv_iter_cube <- cubelyr::as.tbl_cube(surv_iter)
grow_iter_cube <- cubelyr::as.tbl_cube(grow_iter)
flw_iter_cube <- cubelyr::as.tbl_cube(flw_iter)
fert_iter_cube <- cubelyr::as.tbl_cube(fert_iter)
spike_iter_cube <- cubelyr::as.tbl_cube(spike_iter)

surv_iter_df <- as_tibble(surv_iter_cube) %>% 
  left_join(spei_range_df)
grow_iter_df <- as_tibble(grow_iter_cube) %>% 
  left_join(spei_range_df)
flw_iter_df <- as_tibble(flw_iter_cube) %>% 
  left_join(spei_range_df)
fert_iter_df <- as_tibble(fert_iter_cube) %>% 
  left_join(spei_range_df)
spike_iter_df <- as_tibble(spike_iter_cube) %>% 
  left_join(spei_range_df)

surv_spei_mean <- surv_iter_df %>% 
  group_by(spei_value,Endo,Species) %>% 
  summarize(mean_value = mean(surv_iter))
grow_spei_mean <- grow_iter_df %>% 
  group_by(spei_value,Endo,Species) %>% 
  summarize(mean_value = mean(grow_iter))
flw_spei_mean <- flw_iter_df %>% 
  group_by(spei_value,Endo,Species) %>% 
  summarize(mean_value = mean(flw_iter))
fert_spei_mean <- fert_iter_df %>% 
  group_by(spei_value,Endo,Species) %>% 
  summarize(mean_value = mean(fert_iter))
spike_spei_mean <- spike_iter_df %>% 
  group_by(spei_value,Endo,Species) %>% 
  summarize(mean_value = mean(spike_iter))

bin_by_year_spei <- function(df_raw, vr, nbins){
  require(tidyverse)
  df_raw <- ungroup(df_raw)
  spei_bin_df <- df_raw %>% 
    rename(Endo = endo_01, Year = year_t, Species = species) %>%
    group_by(Year, Endo, Species) %>% 
    summarise(mean_spei = mean((spei12),na.rm=T),
              mean_vr = mean({{vr}},na.rm=T),
              samplesize = n()) %>% 
    mutate(Endo = case_when(Endo == 0 ~ "Eminus", 
                            Endo == 1 ~ "Eplus"))
  
  return(spei_bin_df)
}

surv_spei_binned <- bin_by_year_spei(LTREB_data_forsurv, surv_t1)
grow_spei_binned <- bin_by_year_spei(LTREB_data_forgrow, size_t1)
flw_spei_binned <- bin_by_year_spei(LTREB_data_forflw, FLW_STAT_T1)
fert_spei_binned <- bin_by_year_spei(LTREB_data_forfert, FLW_COUNT_T1)
spike_spei_binned <- bin_by_year_spei(LTREB_data_forspike, spike_count_t1)

#The plots
surv_speiplot <- ggplot()+
  geom_line(data = surv_iter_df, aes(x = spei_value, y = surv_iter, color = Endo, group = interaction(Endo,Iteration)), alpha = .1) +
  geom_line(data = surv_spei_mean, aes(x = spei_value, y = mean_value, color = Endo, group = Endo), lwd = 1)+
  geom_point(data = surv_spei_binned, aes(x = mean_spei, y = mean_vr, color = Endo, lwd = samplesize))+
  scale_color_manual(values =  endophyte_color_scheme[c(2,4)])+
  facet_wrap(~Species, scales = "free",ncol = 1) +
  theme_classic() + theme(strip.background = element_blank()) + labs(title = "Adult Survival", subtitle = "Mean endophyte effect with 500 posteriors draws")
# surv_speiplot

grow_speiplot <- ggplot()+
  geom_line(data = grow_iter_df, aes(x = spei_value, y = grow_iter, color = Endo, group = interaction(Endo,Iteration)), alpha = .05) +
  geom_line(data = grow_spei_mean, aes(x = spei_value, y = mean_value, color = Endo, group = Endo), lwd = 1)+
  geom_point(data = grow_spei_binned, aes(x = mean_spei, y = mean_vr, color = Endo, lwd = samplesize))+
  scale_color_manual(values =  endophyte_color_scheme[c(2,4)])+
  facet_wrap(~Species, scales = "free",ncol = 1) +
  theme_classic() + theme(strip.background = element_blank()) + labs(title = "Adult Growth", subtitle = "Mean endophyte effect with 500 posteriors draws")
# grow_speiplot

flw_speiplot <- ggplot()+
  geom_line(data = flw_iter_df, aes(x = spei_value, y = flw_iter, color = Endo, group = interaction(Endo,Iteration)), alpha = .1) +
  geom_line(data = flw_spei_mean, aes(x = spei_value, y = mean_value, color = Endo, group = Endo), lwd = 1)+
  geom_point(data = flw_spei_binned, aes(x = mean_spei, y = mean_vr, color = Endo, lwd = samplesize))+
  scale_color_manual(values =  endophyte_color_scheme[c(2,4)])+
  facet_wrap(~Species, scales = "free",ncol = 1) +
  theme_classic() + theme(strip.background = element_blank()) + labs(title = "Flowering Probability", subtitle = "Mean endophyte effect with 500 posteriors draws")
# flw_speiplot


fert_speiplot <- ggplot()+
  geom_line(data = fert_iter_df, aes(x = spei_value, y = fert_iter, color = Endo, group = interaction(Endo,Iteration)), alpha = .05) +
  geom_line(data = fert_spei_mean, aes(x = spei_value, y = mean_value, color = Endo, group = Endo), lwd = 1)+
  geom_point(data = fert_spei_binned, aes(x = mean_spei, y = mean_vr, color = Endo, lwd = samplesize))+
  scale_color_manual(values =  endophyte_color_scheme[c(2,4)])+
  facet_wrap(~Species, scales = "free",ncol = 1) +
  theme_classic() + theme(strip.background = element_blank()) + labs(title = "# of Flw Tillers", subtitle = "Mean endophyte effect with 500 posteriors draws")
# fert_speiplot

spike_speiplot <- ggplot()+
  geom_line(data = spike_iter_df, aes(x = spei_value, y = spike_iter, color = Endo, group = interaction(Endo,Iteration)), alpha = .05) +
  geom_line(data = spike_spei_mean, aes(x = spei_value, y = mean_value, color = Endo, group = Endo), lwd = 1)+
  geom_point(data = spike_spei_binned, aes(x = mean_spei, y = mean_vr, color = Endo, lwd = samplesize))+
  scale_color_manual(values =  endophyte_color_scheme[c(2,4)])+
  facet_wrap(~Species, scales = "free",ncol = 1) +
  theme_classic() + theme(strip.background = element_blank()) + labs(title = "# of Spikelets/Infl.", subtitle = "Mean endophyte effect with 500 posteriors draws")
# spike_speiplot


spei_effect_fitplot <- surv_speiplot+grow_speiplot+flw_speiplot+fert_speiplot+spike_speiplot + plot_layout( nrow  = 1)
ggsave(spei_effect_fitplot, filename = "spei_effect_fitplot.png", width = 20, height = 25 )

######################################################################
## Making a plot of the climate time series at the site
######################################################################

LTREB_climateforplot <- LTREB_full %>% 
  group_by(species, year) %>% 
  summarize(spei3 = unique(spei3),
            spei12 = unique(spei12),
            annual_temp = unique(annual_temp),
            annual_precip = unique(annual_precip)) %>%  # converting precip units to mm. from tenths of mm. 
  pivot_longer(cols = c(spei3, spei12, annual_temp, annual_precip))

# Making the plot with just climate data for Elymus census, which is tthe middle of the season, and has nearly the highest sd in spei3 (AGPE has higher, but is the most disimilar census timing)
precip_plot <- ggplot(filter(LTREB_climateforplot, species == "ELVI" & name == "annual_precip"))+
  geom_path(aes(x = year, y = value), col = "blue3")+
  scale_y_continuous(name = "Annual Ppt. \n (mm)")+
  theme_minimal()+
  labs(y = "Annual Ppt. (mm)", x = "Year")+
  theme()
# precip_plot

temp_plot<- ggplot(filter(LTREB_climateforplot, species == "ELVI" & name == "annual_temp"))+
  geom_path(aes(x = year, y = value), col = "firebrick1")+
  theme_minimal()+
  scale_y_continuous(name = "Mean Annual Temp. \n(Celsius)")+
  labs(y = "Mean Annual Temp. (Celsius)", x = "Year")+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())
# temp_plot

spei3_plot<- ggplot(filter(LTREB_climateforplot, species == "ELVI" & name == "spei3"))+
  geom_path(aes(x = year, y = value), col = "thistle3")+
  theme_minimal()+
  labs(y = "3-month SPEI", x = "Year")+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())
# spei3_plot

spei12_plot<- ggplot(filter(LTREB_climateforplot, species == "ELVI" & name == "spei12"))+
  geom_path(aes(x = year, y = value), col = "plum4")+
  theme_minimal()+
  labs(y = "12-month SPEI", x = "Year")+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())
# spei12_plot

climate_plot <- spei3_plot + spei12_plot +  temp_plot + precip_plot +
  plot_layout(ncol = 1) + plot_annotation(tag_levels="A")
climate_plot
ggsave(climate_plot, filename = "climate_plot.png", width =8, height =6.5)
