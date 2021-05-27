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



# Get x and y distances from the a and b corner posts
# we calculate the area based on the perimeter; area = sqrt(s(s-a)(s-b)(s-c)) where s is half the perimeter
# Then we get the height (y value) from the area; h = 2*A/b
# Then we can get the x value as one side of a right triangle with the base and one side (x = sqrt(b^2 - y^2))
dist_between_posts <- 3
LTREB_xy <- LTREB_2020_surviving %>% 
  mutate(triangle_perimeter = dist_a + dist_b + dist_between_posts,
         s = triangle_perimeter*.5,
         triangle_area = sqrt(s*(s-dist_a)*(s-dist_b)*(s-dist_between_posts)),
         y = 2*triangle_area/dist_between_posts,
         x = sqrt(dist_a^2-y^2))


maps_aorigin <- ggplot(data = subset(LTREB_xy, plot_fixed == 1))+
  geom_text(aes(x = x, y = y, label = id)) +
  geom_text(aes(x = 0, y = 0, label = "A"), lwd = 4)+
  geom_text(aes(x = dist_between_posts, y = 0, label = "B"), lwd = 4)+
  lims(x = c(0,dist_between_posts), y = c(0, dist_between_posts)) +
  theme_classic()
maps_aorigin

maps_aright <- ggplot(data = subset(LTREB_xy, plot_fixed == 1))+
  geom_text(aes(x = dist_between_posts-x, y = dist_between_posts-y, label = id)) +
  geom_text(aes(x = dist_between_posts, y = dist_between_posts, label = "A"), lwd = 4)+
  geom_text(aes(x = 0, y = dist_between_posts, label = "B"), lwd = 4)+
  lims(x = c(0,dist_between_posts), y = c(0, dist_between_posts)) +
  theme_classic()
maps_aright

maps_aleft <- ggplot(data = subset(LTREB_xy, plot_fixed == 1))+
  geom_text(aes(x = dist_between_posts-x, y = dist_between_posts-y, label = id)) +
  geom_text(aes(x = dist_between_posts, y = dist_between_posts, label = "A"), lwd = 4)+
  geom_text(aes(x = 0, y = dist_between_posts, label = "B"), lwd = 4)+
  scale_x_reverse(c(dist_between_posts, 0))+
  lims(y = c(0, dist_between_posts)) +
  theme_classic()
maps_aleft



