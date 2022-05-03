## Title: Grass endophyte vital rate fit visualizations and plots for all vital rates and all species
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
source("Analyses/endodemog_data_processing.R")
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
  dplyr::select(-FLW_COUNT_T, -FLW_STAT_T, -SPIKE_A_T, -SPIKE_B_T, -SPIKE_C_T, -SPIKE_D_T, -SPIKE_AGPE_MEAN_T, -census_month, -year, -spei1, -spei12, -spei24, -annual_temp, -annual_precip, -endo_status_from_check, -plot_endo_for_check, -endo_mismatch, -dist_a, -dist_b) %>% 
  filter(!is.na(FLW_STAT_T1)) %>% 
  filter(FLW_STAT_T1>0) %>% 
  melt(id.var = c("plot_fixed" ,   "plot_index",         "pos"         ,           "id",
                  "species"       ,         "species_index"  ,        "endo_01",
                  "endo_index"  ,           "origin_01"       ,       "birth" ,
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
                       logsize = LTREB_data_forsurv$logsize_t,
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
                            logsize = LTREB_surv_seedling$logsize_t,
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


flw_data_list <- list(y = LTREB_data_forflw$FLW_STAT_T1,
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

grow_data_list <- list(y = as.integer(LTREB_data_forgrow$size_t1),
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
                       nPlot = max(unique(LTREB_data_forgrow$plot_index)),
                       nEndo =   length(unique(LTREB_data_forgrow$endo_01)))
str(grow_data_list)
seed_grow_data_list <- list(y = as.integer(LTREB_grow_seedling$size_t1),
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
                        origin_01 = LTREB_data_forspike$origin_01)
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

s_to_s_data_list <- read_rds("s_to_s_data_list.rds")

#############################################################################################
####### Read in matrix population functions ------------------
#############################################################################################
# we can use functions from this to set up our parameters for the plots with mean and variance effects
source("Analyses/MPM_functions.R")


#############################################################################################
####### Read in Stan vital rate model outputs ------------------
#############################################################################################

# The stan objects for each vital rate
surv_fit_seedling <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_seedling_surv.rds")
surv_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_surv_woseedling.rds")
grow_fit_seedling <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_seedling_grow_PIG_10000iterations.rds")
grow_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_grow_PIG.rds")
flw_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_flw.rds")
fert_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_fert_PIG.rds")
spike_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_spike_year_plot_nb.rds")
seedmean_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_seed_mean.rds")
stos_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_s_to_s.rds") 

# Pulling out the actuall parameters
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

# Saved y_rep values from the script "vital_rate_analysis.R", "seed_means.R" and 
y_s_sim <- readRDS(file = "yrep_survivalmodel.rds")
y_seed_s_sim <- readRDS(file = "yrep_seedlingsurvivalmodel.rds")
y_f_sim <- readRDS(file = "yrep_floweringmodel.rds")
y_g_sim <- readRDS(file = "yrep_growthPIGmodel.rds")
y_seed_g_sim <- readRDS(file = "yrep_seedlinggrowthPIGmodel.rds")
y_fert_sim <- readRDS(file = "yrep_fertilityPIGmodel.rds")
y_spike_sim <- readRDS(file = "yrep_spikeletNBmodel.rds")
y_seedmean_sim <- readRDS(file = "yrep_seedmeanmodel.rds")
y_recruit_sim <- readRDS(file = "yrep_stosmodel.rds")

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

## Plot for all models vital rate fits #####
# survival
surv_densplot <- ppc_dens_overlay(surv_data_list$y, y_s_sim) + theme_classic() + labs(title = "Adult Survival", x = "Survival status", y = "Density")
# surv_densplot
# ggsave(surv_densplot, filename = "surv_densplot.png", width = 4, height = 4)

mean_s_plot <-   ppc_stat(surv_data_list$y, y_s_sim, stat = "mean")
sd_s_plot <- ppc_stat(surv_data_list$y, y_s_sim, stat = "sd")
skew_s_plot <- ppc_stat(surv_data_list$y, y_s_sim, stat = "skewness")
kurt_s_plot <- ppc_stat(surv_data_list$y, y_s_sim, stat = "Lkurtosis")
surv_moments <- mean_s_plot+sd_s_plot+skew_s_plot+kurt_s_plot + plot_annotation(title = "Survival")
# surv_moments
# ggsave(surv_moments, filename = "surv_momentplot.png", width = 4, height = 4)

#seedliing survival
seedsurv_densplot <- ppc_dens_overlay(seed_surv_data_list$y, y_seed_s_sim) + theme_classic() + labs(title = "Seedling Survival", x = "Survival status", y = "Density")
# seedsurv_densplot
# ggsave(seedsurv_densplot, filename = "seedsurv_densplot.png", width = 4, height = 4)

mean_s_plot <-   ppc_stat(seed_surv_data_list$y, y_seed_s_sim, stat = "mean")
sd_s_plot <- ppc_stat(seed_surv_data_list$y, y_seed_s_sim, stat = "sd")
skew_s_plot <- ppc_stat(seed_surv_data_list$y, y_seed_s_sim, stat = "skewness")
kurt_s_plot <- ppc_stat(seed_surv_data_list$y, y_seed_s_sim, stat = "Lkurtosis")
seedsurv_moments <- mean_s_plot+sd_s_plot+skew_s_plot+kurt_s_plot + plot_annotation(title = "Seedling Survival")
# seedsurv_moments
# ggsave(seedsurv_moments, filename = "seedsurv_momentsplot.png", width = 4, height = 4)

# Flowering
flw_densplot <- ppc_dens_overlay(flw_data_list$y, y_f_sim) + theme_classic() + labs(title = "Flowering", x = "Flowering status", y = "Density")
# flw_densplot
# ggsave(flw_densplot, filename = "flw_densplot.png", width = 4, height = 4)

mean_f_plot <-   ppc_stat(flw_data_list$y, y_f_sim, stat = "mean")
sd_f_plot <- ppc_stat(flw_data_list$y, y_f_sim, stat = "sd")
skew_f_plot <- ppc_stat(flw_data_list$y, y_f_sim, stat = "skewness")
kurt_f_plot <- ppc_stat(flw_data_list$y, y_f_sim, stat = "Lkurtosis")
flw_moments <- mean_f_plot+sd_f_plot+skew_f_plot+kurt_f_plot +plot_annotation(title = "Flowering")
# flw_moments
# ggsave(flw_moments, filename = "flw_momentsplot.png", width = 4, height = 4)

# Growth (with the poisson inverse gaussian distribution)
grow_densplot <- ppc_dens_overlay(grow_data_list$y, y_g_sim) + xlim(0,60) + theme_classic() + labs(title = "Adult Growth with PIG", x = "No. of Tillers", y = "Density")
# grow_densplot
# ggsave(grow_densplot, filename = "grow_densplot.png", width = 4, height = 4)

mean_g_plot <-   ppc_stat(grow_data_list$y, y_g_sim, stat = "mean")
sd_g_plot <- ppc_stat(grow_data_list$y, y_g_sim, stat = "sd")
skew_g_plot <- ppc_stat(grow_data_list$y, y_g_sim, stat = "skewness")
kurt_g_plot <- ppc_stat(grow_data_list$y, y_g_sim, stat = "Lkurtosis")
grow_moments <- mean_g_plot+sd_g_plot+skew_g_plot+kurt_g_plot+ plot_annotation(title = "Growth PIG")
# grow_moments
# ggsave(grow_moments, filename = "grow_momentsplot.png", width = 4, height = 4)

# seedling growth (with poisson inverse gaussian distribution)
seedgrow_densplot <- ppc_dens_overlay(seed_grow_data_list$y, y_seed_g_sim) + xlim(0,60) + theme_classic() + labs(title = "Seedling Growth with PIG", x = "No. of Tillers", y = "Density")
# seedgrow_densplot
# ggsave(seedgrow_densplot, filename = "seed_grow_densplot.png", width = 4, height = 4)

mean_seed_g_plot <-   ppc_stat(seed_grow_data_list$y, y_seed_g_sim, stat = "mean")
sd_seed_g_plot <- ppc_stat(seed_grow_data_list$y, y_seed_g_sim, stat = "sd")
skew_seed_g_plot <- ppc_stat(seed_grow_data_list$y, y_seed_g_sim, stat = "skewness")
kurt_seed_g_plot <- ppc_stat(seed_grow_data_list$y, y_seed_g_sim, stat = "Lkurtosis")
seedgrow_moments <- mean_seed_g_plot+sd_seed_g_plot+skew_seed_g_plot+kurt_seed_g_plot+ plot_annotation(title = "Seedling Growth PIG")
# seedgrow_moments
# ggsave(seedgrow_moments, filename = "seedgrow_momentsplot.png", width = 4, height = 4)

# Fertility
fert_densplot <- ppc_dens_overlay(fert_data_list$y, y_fert_sim) + xlim(0,40) + ggtitle("Fertility with PIG")
# fert_densplot
# ggsave(fert_densplot, filename = "fert_densplot_withPIG.png", width = 4, height = 4)

# This doesn't fit all the moments quite as well, but it's likely good enough
mean_fert_plot <-   ppc_stat(fert_data_list$y, y_fert_sim, stat = "mean")
sd_fert_plot <- ppc_stat(fert_data_list$y, y_fert_sim, stat = "sd")
skew_fert_plot <- ppc_stat(fert_data_list$y, y_fert_sim, stat = "skewness")
kurt_fert_plot <- ppc_stat(fert_data_list$y, y_fert_sim, stat = "Lkurtosis")
fert_moments <- mean_fert_plot+sd_fert_plot+skew_fert_plot+kurt_fert_plot + plot_annotation(title = "Fertility with PIG")
# fert_moments
# ggsave(fert_moments, filename = "fert_momentsplot_withPIG.png", width = 4, height = 4)

# Spikelets per infl
spike_densplot <- ppc_dens_overlay(spike_data_list$y, y_spike_sim) + xlim(0,250) + ggtitle("Spikelet Count")
# spike_densplot
# ggsave(spike_densplot, filename = "spike_densplot.png", width = 4, height = 4)

mean_spike_plot <-   ppc_stat(spike_data_list$y, y_spike_sim, stat = "mean")
sd_spike_plot <- ppc_stat(spike_data_list$y, y_spike_sim, stat = "sd")
skew_spike_plot <- ppc_stat(spike_data_list$y, y_spike_sim, stat = "skewness")
kurt_spike_plot <- ppc_stat(spike_data_list$y, y_spike_sim, stat = "Lkurtosis")
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
recruit_densplot <- ppc_dens_overlay(s_to_s_data_list$tot_recruit_t1, y_recruit_sim) +xlim(0,75) + labs(title = "Recruitment", x = "Successful Germination", y = "Density") 
# ggsave(recruit_densplot, filename = "recruit_densplot.png", width = 4, height = 4)

mean_stos_plot <-   ppc_stat(s_to_s_data_list$tot_recruit_t1, y_recruit_sim, stat = "mean")
sd_stos_plot <- ppc_stat(s_to_s_data_list$tot_recruit_t1, y_recruit_sim, stat = "sd")
skew_stos_plot <- ppc_stat(s_to_s_data_list$tot_recruit_t1, y_recruit_sim, stat = "skewness")
kurt_stos_plot <- ppc_stat(s_to_s_data_list$tot_recruit_t1, y_recruit_sim, stat = "Lkurtosis")
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
                       plot_annotation(title = "Vital rate fits and moments with 500 posterior draws")
ggsave(fitsandmoments_plot, filename = "fitsandmoments_plot.png", width = 18, height = 20)



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

PIG_growth_size_ppc <- size_moments_ppc(data = LTREB_data_forgrow,
                                        y_name = "size_t1",
                                        sim = y_g_sim, 
                                        n_bins = 6, 
                                        title = "Growth PIG")
# ggsave(PIG_growth_size_ppc, filename = "PIG_growth_size_pcc.png", width = 4, height = 4)

size_ppc_plot <- (PIG_growth_size_ppc+ plot_annotation(title = 'Growth'))+
  plot_annotation(title = "Size specific vital rate moments")
# size_ppc_plot
ggsave(size_ppc_plot, filename = "size_ppc_plot.png", width = 12, height = 20)

## Traceplots for select parameters from all vital rates for AGPE
surv_trace <- mcmc_trace(surv_fit, pars = c("beta0[1]", "betaendo[1]", "sigmaendo[1]"))+ggtitle("Survival")
seedsurv_trace <- mcmc_trace(surv_fit_seedling, pars = c("beta0[1]", "betaendo[1]", "sigmaendo[1]"))+ggtitle("Seedling Survival")
flw_trace <- mcmc_trace(flw_fit, pars = c("beta0[1]", "betaendo[1]", "sigmaendo[1]"))+ggtitle("Flowering")
grow_trace <- mcmc_trace(grow_fit, pars = c("beta0[1]", "betaendo[1]", "sigmaendo[1]"))+ggtitle("Growth")
seedgrow_trace <- mcmc_trace(grow_fit_seedling, pars = c("beta0[1]", "betaendo[1]", "sigmaendo[1]"))+ggtitle("Seedling Growth")
fert_trace <- mcmc_trace(fert_fit, pars = c("beta0[1]", "betaendo[1]", "sigmaendo[1]"))+ggtitle("Fertility")
spike_trace <- mcmc_trace(spike_fit, pars = c("beta0[1]", "betaendo[1]", "sigmaendo[1]"))+ggtitle("Spikelets per inflorescence")
seedmean_trace <- mcmc_trace(seedmean_fit, pars = c("beta0[1]", "betaendo[1]"))+ggtitle("Mean Seeds per spikelet")
stos_trace <- mcmc_trace(stos_fit, pars = c("beta0[1]", "betaendo[1]", "sigmaendo[1]"))+ggtitle("Recruitment")

vr_traceplots <- (surv_trace)/
                 (seedsurv_trace)/
                 (grow_trace)/
                 (seedgrow_trace)/
                 (flw_trace)/
                 (fert_trace)/
                 (spike_trace)/
                 (seedmean_trace)/
                 (stos_trace) + plot_annotation(title = "Traceplots for select parameters from all vital rates (AGPE)")
ggsave(vr_traceplots, filename = "vr_traceplots.png", width = 25, height = 20)


## Plots for mean and year variance endophyte effects on all vital rates, all species ####
mean_species_key <- c("betaendo[1]" = "AGPE", "betaendo[2]" = "ELRI", "betaendo[3]" = "ELVI", "betaendo[4]" = "FESU", "betaendo[5]" = "LOAR", "betaendo[6]" = "POAL", "betaendo[7]" = "POSY")
sd_species_key <- c("sigmaendo[1]" = "AGPE", "sigmaendo[2]" = "ELRI", "sigmaendo[3]" = "ELVI", "sigmaendo[4]" = "FESU", "sigmaendo[5]" = "LOAR", "sigmaendo[6]" = "POAL", "sigmaendo[7]" = "POSY")

# effect of endophytte on mean (betaendo)
surv_endomean_posteriors <-  mcmc_areas(surv_fit, prob = 0.8, regex_pars = c("betaendo"))+labs(title = "Adult Survival", subtitle = "Endophyte effect on mean with 80% credible intervals") + scale_y_discrete(labels = mean_species_key)
seedsurv_endomean_posteriors <-  mcmc_areas(surv_fit_seedling, prob = 0.8, regex_pars = c("betaendo"))+labs(title = "Adult Survival", subtitle = "Endophyte effect on mean with 80% credible intervals") + scale_y_discrete(labels = mean_species_key)
grow_endomean_posteriors <- mcmc_areas(grow_fit, prob = 0.8, regex_pars = c("betaendo"))+labs(title = "Adult Growth", subtitle = "Endophyte effect on mean with 80% credible intervals") + scale_y_discrete(labels = mean_species_key)
seedgrow_endomean_posteriors <- mcmc_areas(grow_fit_seedling, prob = 0.8, regex_pars = c("betaendo"))+labs(title = "Seedling Growth", subtitle = "Endophyte effect on mean with 80% credible intervals") + scale_y_discrete(labels = mean_species_key)
flw_endomean_posteriors <- mcmc_areas(flw_fit, prob = 0.8, regex_pars = c("betaendo"))+labs(title = "Flowering", subtitle = "Endophyte effect on mean with 80% credible intervals") + scale_y_discrete(labels = mean_species_key)
fert_endomean_posteriors <- mcmc_areas(fert_fit, prob = 0.8, regex_pars = c("betaendo"))+labs(title = "Fertility", subtitle = "Endophyte effect on mean with 80% credible intervals") + scale_y_discrete(labels = mean_species_key)
spike_endomean_posteriors <- mcmc_areas(spike_fit, prob = 0.8, regex_pars = c("betaendo"))+labs(title = "Spikelets per infl.", subtitle = "Endophyte effect on mean with 80% credible intervals") + scale_y_discrete(labels = mean_species_key)
seedmean_endomean_posteriors <- mcmc_areas(seedmean_fit, prob = 0.8, regex_pars = c("betaendo"))+labs(title = "Mean seeds per spikelet", subtitle = "Endophyte effect on mean with 80% credible intervals") + scale_y_discrete(labels = mean_species_key)
stos_endomean_posteriors <- mcmc_areas(stos_fit, prob = 0.8, regex_pars = c("betaendo"))+labs(title = "Germination", subtitle = "Endophyte effect on mean with 80% credible intervals") + scale_y_discrete(labels = mean_species_key)

# effect of endophyte on year standard deviation (seed mean does not have an endophyte effect on variance)
surv_endosd_posteriors <-  mcmc_areas(surv_fit, prob = 0.8, regex_pars = c("sigmaendo"))+labs(title = "Adult Survival", subtitle = "Endophyte effect on SD with 80% credible intervals") + scale_y_discrete(labels = sd_species_key)
seedsurv_endosd_posteriors <-  mcmc_areas(surv_fit_seedling, prob = 0.8, regex_pars = c("sigmaendo"))+labs(title = "Adult Survival", subtitle = "Endophyte effect on SD with 80% credible intervals") + scale_y_discrete(labels = sd_species_key)
grow_endosd_posteriors <- mcmc_areas(grow_fit, prob = 0.8, regex_pars = c("sigmaendo"))+labs(title = "Adult Growth", subtitle = "Endophyte effect on SD with 80% credible intervals") + scale_y_discrete(labels = sd_species_key)
seedgrow_endosd_posteriors <- mcmc_areas(grow_fit_seedling, prob = 0.8, regex_pars = c("sigmaendo"))+labs(title = "Seedling Growth", subtitle = "Endophyte effect on SD with 80% credible intervals") + scale_y_discrete(labels = sd_species_key)
flw_endosd_posteriors <- mcmc_areas(flw_fit, prob = 0.8, regex_pars = c("sigmaendo"))+labs(title = "Flowering", subtitle = "Endophyte effect on SD with 80% credible intervals") + scale_y_discrete(labels = sd_species_key)
fert_endosd_posteriors <- mcmc_areas(fert_fit, prob = 0.8, regex_pars = c("sigmaendo"))+labs(title = "Fertility", subtitle = "Endophyte effect on SD with 80% credible intervals") + scale_y_discrete(labels = sd_species_key)
spike_endosd_posteriors <- mcmc_areas(spike_fit, prob = 0.8, regex_pars = c("sigmaendo"))+labs(title = "Spikelets per infl.", subtitle = "Endophyte effect on SD with 80% credible intervals") + scale_y_discrete(labels = sd_species_key)
stos_endosd_posteriors <- mcmc_areas(stos_fit, prob = 0.8, regex_pars = c("sigmaendo"))+labs(title = "Germination", subtitle = "Endophyte effect on SD with 80% credible intervals") + scale_y_discrete(labels = sd_species_key)



endomeanandvar_posteriors <- (surv_endomean_posteriors + surv_endosd_posteriors)/
                       (seedsurv_endomean_posteriors + seedsurv_endosd_posteriors)/
                       (grow_endomean_posteriors + grow_endosd_posteriors)/
                       (seedgrow_endomean_posteriors + seedgrow_endosd_posteriors)/
                       (flw_endomean_posteriors + flw_endosd_posteriors)/
                       (fert_endomean_posteriors + fert_endosd_posteriors)/
                       (spike_endomean_posteriors + spike_endosd_posteriors)/
                       (seedmean_endomean_posteriors + plot_spacer())/
                       (stos_endomean_posteriors + stos_endosd_posteriors) + plot_annotation(title = "Endophyte effect on mean and interannual SD across vital rates by species")
  
ggsave(endomeanandvar_posteriors, filename = "endomeanandvar_posteriors.png", width = 12, height = 30)

## Plots for all vital rates, all species by size with data #####
max_size <- LTREB_full %>% 
  dplyr::select(species,species_index, size_t) %>% 
  filter(!is.na(size_t)) %>% 
  group_by(species, species_index) %>% 
  summarise(actual_max_size = max(size_t),
            max_size = quantile(size_t,probs=0.975))

n_post_draws <- 500
post_draws <- sample.int(7500, n_post_draws)

x_seq_length <- 100
x_seq <- array(dim= c(x_seq_length,7))
for(s in 1:7){
 x_seq[,s] <-  seq(from = 1, to = filter(max_size, species_index == s)$actual_max_size, length.out = 100)
}

surv_iter <- grow_iter <- flw_iter <- fert_iter <-array(dim = c(length(x_seq[,1]),2,7, n_post_draws))
surv_mean <- grow_mean <- flw_mean <- fert_mean <-array(dim = c(length(x_seq[,1]),2,7,3))

survyear_iter <- growyear_iter <- flwyear_iter <- fertyear_iter <-array(dim = c(length(x_seq[,1]),2,7, (length(unique(LTREB_full$year_t_index))),n_post_draws))
survyear_mean <- growyear_mean <- flwyear_mean <- fertyear_mean <-array(dim = c(length(x_seq[,1]),2,7,(length(unique(LTREB_full$year_t_index))),3))


sx<-function(x,params){
  invlogit(params$surv_int + params$surv_slope*log(x))
}
gx <- function(x,params){
  exp(params$grow_int + params$grow_slope*log(x))
}
flwx <- function(x,params){
  invlogit(params$flow_int + params$flow_slope*log(x))
}
fertx <- function(x,params){
  exp(params$fert_int + params$fert_slope*log(x))
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
                      surv_par=surv_par,
                      surv_sdlg_par = surv_sdlg_par,
                      grow_par=grow_par,
                      grow_sdlg_par = grow_sdlg_par,
                      flow_par=flow_par,
                      fert_par=fert_par,
                      spike_par=spike_par,
                      seed_par=seed_par,
                      recruit_par=recruit_par), x = x_seq[,s])
grow_iter[,e,s,i] <- gx(make_params(species=s,
                                    endo_mean=(e-1),
                                    endo_var=(e-1),
                                    original = mean(LTREB_full$origin_01), # should be =1 to represent recruit
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
                                    recruit_par=recruit_par), x = x_seq[,s])
flw_iter[,e,s,i] <- flwx(make_params(species=s,
                                    endo_mean=(e-1),
                                    endo_var=(e-1),
                                    original = mean(LTREB_full$origin_01), # should be =1 to represent recruit
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
                                    recruit_par=recruit_par), x = x_seq[,s])
fert_iter[,e,s,i] <- fertx(make_params(species=s,
                                    endo_mean=(e-1),
                                    endo_var=(e-1),
                                    original = mean(LTREB_full$origin_01), # should be =1 to represent recruit
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
                                    recruit_par=recruit_par), x = x_seq[,s])

    }
  }
}

for(x in 1:length(x_seq[,1])){
  for(e in 1:2){
    for(s in 1:7){
surv_mean[x,e,s,1] <- mean(surv_iter[x,e,s,], na.rm = T)
surv_mean[x,e,s,2:3] <- quantile(surv_iter[x,e,s,], probs = c(.1,.9), na.rm = T)

grow_mean[x,e,s,1] <- mean(grow_iter[x,e,s,], na.rm = T)
grow_mean[x,e,s,2:3] <- quantile(grow_iter[x,e,s,], probs = c(.1,.9), na.rm = T)

flw_mean[x,e,s,1] <- mean(flw_iter[x,e,s,], na.rm = T)
flw_mean[x,e,s,2:3] <- quantile(flw_iter[x,e,s,], probs = c(.1,.9), na.rm = T)

fert_mean[x,e,s,1] <- mean(fert_iter[x,e,s,], na.rm = T)
fert_mean[x,e,s,2:3] <- quantile(fert_iter[x,e,s,], probs = c(.1,.9), na.rm = T)
    }
  }
}

species_list <- c("AGPE", "ELRI", "ELVI", "FESU", "LOAR", "POAL", "POSY")

dimnames(surv_mean) <- list(paste0("size", 1:length(x_seq[,1])), c("Eminus","Eplus"), species_list, c("mean","twenty","eighty"))
dimnames(grow_mean) <- list(paste0("size", 1:length(x_seq[,1])), c("Eminus","Eplus"), species_list, c("mean","twenty","eighty"))
dimnames(flw_mean) <- list(paste0("size", 1:length(x_seq[,1])), c("Eminus","Eplus"), species_list, c("mean","twenty","eighty"))
dimnames(fert_mean) <- list(paste0("size", 1:length(x_seq[,1])), c("Eminus","Eplus"), species_list, c("mean","twenty","eighty"))

# Now the same thing for each year specific vital rate
for(i in 1:length(post_draws)){
  for(e in 1:2){
    for(s in 1:7){
      for(y in 1:(length(unique(LTREB_full$year_t_index)))){
      
      survyear_iter[,e,s,y,i] <- sx(make_params(species=s,
                                          endo_mean=(e-1),
                                          endo_var=(e-1),
                                          original = mean(LTREB_full$origin_01), # should be =1 to represent recruit
                                          draw=post_draws[i],
                                          max_size=max_size,
                                          rfx=T,
                                          year = y,
                                          surv_par=surv_par,
                                          surv_sdlg_par = surv_sdlg_par,
                                          grow_par=grow_par,
                                          grow_sdlg_par = grow_sdlg_par,
                                          flow_par=flow_par,
                                          fert_par=fert_par,
                                          spike_par=spike_par,
                                          seed_par=seed_par,
                                          recruit_par=recruit_par), x = x_seq[,s])
      growyear_iter[,e,s,y,i] <- gx(make_params(species=s,
                                          endo_mean=(e-1),
                                          endo_var=(e-1),
                                          original = mean(LTREB_full$origin_01), # should be =1 to represent recruit
                                          draw=post_draws[i],
                                          max_size=max_size,
                                          rfx=T,
                                          year = y,
                                          surv_par=surv_par,
                                          surv_sdlg_par = surv_sdlg_par,
                                          grow_par=grow_par,
                                          grow_sdlg_par = grow_sdlg_par,
                                          flow_par=flow_par,
                                          fert_par=fert_par,
                                          spike_par=spike_par,
                                          seed_par=seed_par,
                                          recruit_par=recruit_par), x = x_seq[,s])
      flwyear_iter[,e,s,y,i] <- flwx(make_params(species=s,
                                           endo_mean=(e-1),
                                           endo_var=(e-1),
                                           original = mean(LTREB_full$origin_01), # should be =1 to represent recruit
                                           draw=post_draws[i],
                                           max_size=max_size,
                                           rfx=T,
                                           year = y,
                                           surv_par=surv_par,
                                           surv_sdlg_par = surv_sdlg_par,
                                           grow_par=grow_par,
                                           grow_sdlg_par = grow_sdlg_par,
                                           flow_par=flow_par,
                                           fert_par=fert_par,
                                           spike_par=spike_par,
                                           seed_par=seed_par,
                                           recruit_par=recruit_par), x = x_seq[,s])
      fertyear_iter[,e,s,y,i] <- fertx(make_params(species=s,
                                             endo_mean=(e-1),
                                             endo_var=(e-1),
                                             original = mean(LTREB_full$origin_01), # should be =1 to represent recruit
                                             draw=post_draws[i],
                                             max_size=max_size,
                                             rfx=T,
                                             year = y,
                                             surv_par=surv_par,
                                             surv_sdlg_par = surv_sdlg_par,
                                             grow_par=grow_par,
                                             grow_sdlg_par = grow_sdlg_par,
                                             flow_par=flow_par,
                                             fert_par=fert_par,
                                             spike_par=spike_par,
                                             seed_par=seed_par,
                                             recruit_par=recruit_par), x = x_seq[,s])

      }
    }
  }
}
        
for(x in 1:length(x_seq[,1])){
  for(e in 1:2){
    for(s in 1:7){
      for(y in 1:(length(unique(LTREB_full$year_t_index)))){
      survyear_mean[x,e,s,y,1] <- mean(survyear_iter[x,e,s,y,], na.rm = T)
      survyear_mean[x,e,s,y,2:3] <- quantile(survyear_iter[x,e,s,y,], probs = c(.1,.9), na.rm = T)

      growyear_mean[x,e,s,y,1] <- mean(growyear_iter[x,e,s,y,], na.rm = T)
      growyear_mean[x,e,s,y,2:3] <- quantile(growyear_iter[x,e,s,y,], probs = c(.1,.9), na.rm = T)

      flwyear_mean[x,e,s,y,1] <- mean(flwyear_iter[x,e,s,y,], na.rm = T)
      flwyear_mean[x,e,s,y,2:3] <- quantile(flwyear_iter[x,e,s,y,], probs = c(.1,.9), na.rm = T)

      fertyear_mean[x,e,s,y,1] <- mean(fertyear_iter[x,e,s,y,], na.rm = T)
      fertyear_mean[x,e,s,y,2:3] <- quantile(fertyear_iter[x,e,s,y,], probs = c(.1,.9), na.rm = T)
      }
    }
  }
}

dimnames(survyear_mean) <- list(paste0("size", 1:length(x_seq[,1])), c("Eminus","Eplus"), species_list,c(2007:2020), c("mean","twenty","eighty"))
dimnames(growyear_mean) <- list(paste0("size", 1:length(x_seq[,1])), c("Eminus","Eplus"), species_list,c(2007:2020), c("mean","twenty","eighty"))
dimnames(flwyear_mean) <- list(paste0("size", 1:length(x_seq[,1])), c("Eminus","Eplus"), species_list,c(2007:2020), c("mean","twenty","eighty"))
dimnames(fertyear_mean) <- list(paste0("size", 1:length(x_seq[,1])), c("Eminus","Eplus"), species_list,c(2007:2020), c("mean","twenty","eighty"))


#Now I'm gonna make these into  tidy dataframes for plotting

dimnames(x_seq) <- list(paste0("x_size", 1:length(x_seq[,1])), paste0(species_list))

x_seq_df <- as_tibble(x_seq) %>%
  mutate(no_row = row_number()) %>% 
  pivot_longer(cols = -no_row,
               names_to = "Species",
               values_to = "x_seq") %>% 
  mutate(log_x_seq = log(x_seq))


surv_mean_df <- as_tibble(surv_mean) %>%    
  mutate(no_row = row_number()) %>% 
  pivot_longer(cols = starts_with("E"),
               values_to = c("surv") ) %>% 
  separate(name, c("Endo", "Species","quantile")) %>% 
  pivot_wider(id_cols = c("no_row", "Endo", "Species"), names_from = "quantile", values_from = "surv") %>% 
  left_join(x_seq_df)

survyear_mean_df <- as_tibble(survyear_mean) %>%    
  mutate(no_row = row_number()) %>% 
  pivot_longer(cols = -no_row,
               values_to = c("surv")) %>% 
  separate(name, c("Endo", "Species", "Year", "quantile")) %>% 
  pivot_wider(id_cols = c("no_row", "Endo", "Species", "Year"), names_from = "quantile", values_from = "surv") %>% 
  left_join(x_seq_df)

grow_mean_df <- as_tibble(grow_mean) %>%    
  mutate(no_row = row_number()) %>% 
  pivot_longer(cols = starts_with("E"),
               values_to = c("size_t1") ) %>% 
  separate(name, c("Endo", "Species","quantile")) %>% 
  pivot_wider(id_cols = c("no_row", "Endo", "Species"), names_from = "quantile", values_from = "size_t1") %>% 
  left_join(x_seq_df)

growyear_mean_df <- as_tibble(growyear_mean) %>%    
  mutate(no_row = row_number()) %>% 
  pivot_longer(cols = starts_with("E"),
               values_to = c("size_t1") ) %>% 
  separate(name, c("Endo", "Species", "Year", "quantile")) %>% 
  pivot_wider(id_cols = c("no_row", "Endo", "Species", "Year"), names_from = "quantile", values_from = "size_t1") %>% 
  left_join(x_seq_df)

flw_mean_df <- as_tibble(flw_mean) %>%    
  mutate(no_row = row_number()) %>% 
  pivot_longer(cols = starts_with("E"),
               values_to = c("flw_t1") ) %>% 
  separate(name, c("Endo", "Species","quantile")) %>% 
  pivot_wider(id_cols = c("no_row", "Endo", "Species"), names_from = "quantile", values_from = "flw_t1") %>% 
  left_join(x_seq_df)

flwyear_mean_df <- as_tibble(flwyear_mean) %>%    
  mutate(no_row = row_number()) %>% 
  pivot_longer(cols = starts_with("E"),
               values_to = c("flw_t1") ) %>% 
  separate(name, c("Endo", "Species", "Year", "quantile")) %>% 
  pivot_wider(id_cols = c("no_row", "Endo", "Species", "Year"), names_from = "quantile", values_from = "flw_t1") %>% 
  left_join(x_seq_df)

fert_mean_df <- as_tibble(fert_mean) %>%    
  mutate(no_row = row_number()) %>% 
  pivot_longer(cols = starts_with("E"),
               values_to = c("fert_t1") ) %>% 
  separate(name, c("Endo", "Species","quantile")) %>% 
  pivot_wider(id_cols = c("no_row", "Endo", "Species"), names_from = "quantile", values_from = "fert_t1") %>% 
  left_join(x_seq_df)

fertyear_mean_df <- as_tibble(fertyear_mean) %>%    
  mutate(no_row = row_number()) %>% 
  pivot_longer(cols = starts_with("E"),
               values_to = c("fert_t1") ) %>% 
  separate(name, c("Endo", "Species", "Year", "quantile")) %>% 
  pivot_wider(id_cols = c("no_row", "Endo", "Species", "Year"), names_from = "quantile", values_from = "fert_t1") %>% 
  left_join(x_seq_df)

# Bin Data by size and then by year for plotting

bin_by_size_t <- function(df_raw, vr, nbins){
  require(tidyverse)
  df_raw <- ungroup(df_raw)
  size_bin_df <- df_raw %>% 
    rename(Endo = endo_01, Year = year_t, Species = species) %>%
    mutate(size_bin = cut(logsize_t, breaks = nbins)) %>%
    group_by(size_bin, Endo, Species) %>%
    summarise(mean_size = mean((logsize_t),na.rm=T),
              mean_vr = mean({{vr}},na.rm=T),
              samplesize = n()) %>% 
    mutate(Endo = case_when(Endo == 0 ~ "Eminus", 
                           Endo == 1 ~ "Eplus"))
  
  return(size_bin_df)
}
bin_by_size_t1 <- function(df_raw, vr, nbins){
  require(tidyverse)
  df_raw <- ungroup(df_raw)
  size_bin_df <- df_raw %>% 
    rename(Endo = endo_01, Year = year_t, Species = species) %>%
    mutate(size_bin = cut(logsize_t1, breaks = nbins)) %>%
    group_by(size_bin, Endo, Species) %>%
    summarise(mean_size = mean((logsize_t1),na.rm=T),
              mean_vr = mean({{vr}},na.rm=T),
              samplesize = n()) %>% 
    mutate(Endo = case_when(Endo == 0 ~ "Eminus", 
                           Endo == 1 ~ "Eplus"))
  
  return(size_bin_df)
}

bin_by_year_size_t <- function(df_raw, vr, nbins){
  require(tidyverse)
  df_raw <- ungroup(df_raw)
  size_bin_df <- df_raw %>% 
    rename(Endo = endo_01, Year = year_t, Species = species) %>%
    mutate(size_bin = cut(logsize_t, breaks = nbins)) %>%
    group_by(size_bin, Year, Endo, Species) %>%
    summarise(mean_size = mean((logsize_t),na.rm=T),
              mean_vr = mean({{vr}},na.rm=T),
              samplesize = n()) %>% 
    mutate(Endo = case_when(Endo == 0 ~ "Eminus", 
                            Endo == 1 ~ "Eplus"))
  
  return(size_bin_df)
}
bin_by_year_size_t1 <- function(df_raw, vr, nbins){
  require(tidyverse)
  df_raw <- ungroup(df_raw)
  size_bin_df <- df_raw %>% 
    rename(Endo = endo_01, Year = year_t, Species = species) %>%
    mutate(size_bin = cut(logsize_t1, breaks = nbins)) %>%
    group_by(size_bin, Year, Endo, Species) %>%
    summarise(mean_size = mean((logsize_t1),na.rm=T),
              mean_vr = mean({{vr}},na.rm=T),
              samplesize = n()) %>% 
    mutate(Endo = case_when(Endo == 0 ~ "Eminus", 
                            Endo == 1 ~ "Eplus"))
  
  return(size_bin_df)
}


surv_sizebin <- bin_by_size_t(LTREB_data_forsurv,vr = surv_t1, nbins = 20)
seedsurv_sizebin <- bin_by_size_t(LTREB_surv_seedling,vr = surv_t1, nbins = 20)
grow_sizebin <- bin_by_size_t(LTREB_data_forgrow,vr = size_t1, nbins = 20)
seedgrow_sizebin <- bin_by_size_t(LTREB_grow_seedling, vr = size_t1, nbins = 20)
flw_sizebin <- bin_by_size_t1(LTREB_data_forflw, vr = FLW_STAT_T1,nbins = 20)
fert_sizebin <- bin_by_size_t1(LTREB_data_forfert,vr = FLW_COUNT_T1, nbins = 20)
spike_sizebin <- bin_by_size_t1(LTREB_data_forspike, vr = spike_count_t1, nbins = 20)

surv_yearsizebin <- bin_by_year_size_t(LTREB_data_forsurv,vr = surv_t1, nbins = 20)
seedsurv_yearsizebin <- bin_by_year_size_t(LTREB_surv_seedling,vr = surv_t1, nbins = 20)
grow_yearsizebin <- bin_by_year_size_t(LTREB_data_forgrow,vr = size_t1, nbins = 20)
seedgrow_yearsizebin <- bin_by_year_size_t(LTREB_grow_seedling, vr = size_t1, nbins = 20)
flw_yearsizebin <- bin_by_year_size_t1(LTREB_data_forflw, vr = FLW_STAT_T1,nbins = 20)
fert_yearsizebin <- bin_by_year_size_t1(LTREB_data_forfert,vr = FLW_COUNT_T1, nbins = 20)
spike_yearsizebin <- bin_by_year_size_t1(LTREB_data_forspike, vr = spike_count_t1, nbins = 20)


#The plots
surv_meanplot <- ggplot()+
  geom_point(data = surv_sizebin, aes(x = mean_size, y = mean_vr, size = samplesize, shape = Endo), alpha = .5) +
  geom_ribbon(data = surv_mean_df, aes(x = log_x_seq, ymin = twenty, ymax = eighty, fill = Endo), alpha = .3)+
  geom_line(data = surv_mean_df, aes(x = log_x_seq, y = mean, linetype = Endo)) +
  scale_shape_manual(values = c(19,1))+ scale_fill_manual(values = c( endophyte_color_scheme[5], endophyte_color_scheme[3]))+
  facet_wrap(~Species, scales = "free", ncol = 1) + 
  theme_classic() + theme(strip.background = element_blank()) + labs(title = "Adult Survival", subtitle = "Mean endophyte effect with 80% credible intervals")
# surv_meanplot

surv_yearplot <- ggplot()+
  geom_point(data = surv_yearsizebin, aes(x = mean_size, y = mean_vr, size = samplesize, col = as.factor(Year), shape = Endo), alpha = .5) +
  # geom_ribbon(data = survyear_mean_df, aes(x = log_x_seq, ymin = twenty, ymax = eighty, fill = Endo), alpha = .3)+
  geom_line(data = survyear_mean_df, aes(x = log_x_seq, y = mean, linetype = Endo, col = Year)) +
  scale_shape_manual(values = c(19,1))+ 
  scale_color_manual(values = yearcolors)+
  facet_wrap(~Species + Endo, scales = "free", ncol = 2) + 
  theme_classic() + theme(strip.background = element_blank()) + labs(title = "Adult Survival", subtitle = "Endophyte effect on interannual variation")
# surv_yearplot

grow_meanplot <- ggplot()+
  geom_point(data = grow_sizebin, aes(x = mean_size, y = mean_vr, size = samplesize, shape = Endo), alpha = .5) +
  geom_ribbon(data = grow_mean_df, aes(x = log_x_seq, ymin = twenty, ymax = eighty, fill = Endo), alpha = .3)+
  geom_line(data = grow_mean_df, aes(x = log_x_seq, y = mean, linetype = Endo)) +
  scale_shape_manual(values = c(19,1))+ scale_fill_manual(values = c( endophyte_color_scheme[5], endophyte_color_scheme[3]))+
  facet_wrap(~Species, scales = "free", ncol = 1) +
  theme_classic() + theme(strip.background = element_blank()) + labs(title = "Adult Growth", subtitle = "Mean endophyte effect with 80% credible intervals")
# grow_meanplot

grow_yearplot <- ggplot()+
  geom_point(data = grow_yearsizebin, aes(x = mean_size, y = mean_vr, size = samplesize, col = as.factor(Year), shape = Endo), alpha = .5) +
  # geom_ribbon(data = growyear_mean_df, aes(x = log_x_seq, ymin = twenty, ymax = eighty, fill = Endo), alpha = .3)+
  geom_line(data = growyear_mean_df, aes(x = log_x_seq, y = mean, linetype = Endo, col = Year)) +
  scale_shape_manual(values = c(19,1))+ 
  scale_color_manual(values = yearcolors)+
  facet_wrap(~Species + Endo, scales = "free", ncol = 2) + 
  theme_classic() + theme(strip.background = element_blank()) + labs(title = "Adult Growth", subtitle = "Endophyte effect on interannual variation")
# grow_yearplot

flw_meanplot <- ggplot()+
  geom_point(data = flw_sizebin, aes(x = mean_size, y = mean_vr, size = samplesize, shape = Endo), alpha = .5) +
  geom_ribbon(data = flw_mean_df, aes(x = log_x_seq, ymin = twenty, ymax = eighty, fill = Endo), alpha = .3)+
  geom_line(data = flw_mean_df, aes(x = log_x_seq, y = mean, linetype = Endo)) +
  scale_shape_manual(values = c(19,1))+ scale_fill_manual(values = c( endophyte_color_scheme[5], endophyte_color_scheme[3]))+
  facet_wrap(~Species, scales = "free", ncol = 1) + 
  theme_classic() + theme(strip.background = element_blank())+ labs(title = "Flowering", subtitle = "Mean endophyte effect with 80% credible intervals")
# flw_meanplot

flw_yearplot <- ggplot()+
  geom_point(data = flw_yearsizebin, aes(x = mean_size, y = mean_vr, size = samplesize, col = as.factor(Year), shape = Endo), alpha = .5) +
  # geom_ribbon(data = flwyear_mean_df, aes(x = log_x_seq, ymin = twenty, ymax = eighty, fill = Endo), alpha = .3)+
  geom_line(data = flwyear_mean_df, aes(x = log_x_seq, y = mean, linetype = Endo, col = Year)) +
  scale_shape_manual(values = c(19,1))+ 
  scale_color_manual(values = yearcolors)+
  facet_wrap(~Species + Endo, scales = "free", ncol = 2) + 
  theme_classic() + theme(strip.background = element_blank()) + labs(title = "Flowering", subtitle = "Endophyte effect on interannual variation")
# flw_yearplot

fert_meanplot <- ggplot()+
  geom_point(data = fert_sizebin, aes(x = mean_size, y = mean_vr, size = samplesize, shape = Endo), alpha = .5) +
  geom_ribbon(data = fert_mean_df, aes(x = log_x_seq, ymin = twenty, ymax = eighty, fill = Endo), alpha = .3)+
  geom_line(data = fert_mean_df, aes(x = log_x_seq, y = mean, linetype = Endo)) +
  scale_shape_manual(values = c(19,1))+ scale_fill_manual(values = c( endophyte_color_scheme[5], endophyte_color_scheme[3]))+
  facet_wrap(~Species, scales = "free", ncol = 1) + 
  theme_classic() + theme(strip.background = element_blank())+ labs(title = "Fertility", subtitle = "Mean endophyte effect with 80% credible intervals")
# fert_meanplot

fert_yearplot <- ggplot()+
  geom_point(data = fert_yearsizebin, aes(x = mean_size, y = mean_vr, size = samplesize, col = as.factor(Year), shape = Endo), alpha = .5) +
  # geom_ribbon(data = fertyear_mean_df, aes(x = log_x_seq, ymin = twenty, ymax = eighty, fill = Endo), alpha = .3)+
  geom_line(data = fertyear_mean_df, aes(x = log_x_seq, y = mean, linetype = Endo, col = Year)) +
  scale_shape_manual(values = c(19,1))+ 
  scale_color_manual(values = yearcolors)+
  facet_wrap(~Species + Endo, scales = "free", ncol = 2) + 
  theme_classic() + theme(strip.background = element_blank()) + labs(title = "Fertility", subtitle = "Endophyte effect on interannual variation")
# fert_yearplot

meaneffect_fitplot <- surv_meanplot+grow_meanplot+flw_meanplot+fert_meanplot + plot_layout( nrow  = 1)
ggsave(meaneffect_fitplot, filename = "meaneffect_fitplot.png", width = 20, height = 25 )  


vareffect_fitplot <- surv_yearplot+grow_yearplot+flw_yearplot+fert_yearplot + plot_layout( nrow  = 1)
ggsave(vareffect_fitplot, filename = "vareffect_fitplot.png", width = 30, height = 25 )  

meanvareffect_fitplot <-  surv_meanplot + surv_yearplot+grow_meanplot +grow_yearplot +flw_meanplot+flw_yearplot+fert_meanplot +fert_yearplot + plot_layout( nrow  = 1,
                                                                                                                                                            widths = c(1,2,1,2,1,2,1,2),
                                                                                                                                                            guides = "collect",
                                                                                                                                                            )
ggsave(meanvareffect_fitplot, filename = "meanvareffect_fitplot.png", width = 40, height = 25 )  
