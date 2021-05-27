## Create maps for LTREB Indiana EndoDemog data collection
## Authors: Josh 
## Purpose: Create maps for LTREB Indiana EndoDemog data collection using data compiled by endodemog_data_processing.R
## Last Update: May 27, 2021
######################################################
library(tidyverse)
library(reshape2)
library(lubridate)
library(readxl)
library(SPEI)


##############################################################################
####### Source in endodemog_data_processing ------------------------------
##############################################################################

source("Analyses/endodemog_data_processing.R")

LTREB_2020_surviving <- LTREB_full %>% 
  filter(surv_t1 == 1 & year_t1 == 2020) %>% 
  filter(!is.na(dist_a) | !is.na(dist_b)) %>% 
  select(plot_fixed, pos, id, species,birth, year_t1, size_t1, dist_a, dist_b) %>% 
  rename(census_year = year_t1, size_2020 = size_t1) %>% 
  group_by(species, plot_fixed)

write_csv(LTREB_2020_surviving, path = "~/Documents/R_projects/Grass-Endophyte-Stochastic-Demography/LTREB_2020_surviving.csv")

