## Create maps for LTREB Indiana EndoDemog data collection
## Authors: Josh 
## Purpose: Create maps for LTREB Indiana EndoDemog data collection using data compiled by endodemog_data_processing.R
## Last Update: May 21, 2022
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

LTREB_2021_surviving <- LTREB_full %>% 
  filter(surv_t1 == 1 & year_t1 == 2021) %>% 
  filter(!is.na(dist_a) | !is.na(dist_b)) %>% 
  select(plot_fixed, plot_index, pos, id, species,birth, year_t1, size_t1, dist_a, dist_b) %>% 
  rename(census_year = year_t1, size_2020 = size_t1) %>% 
  group_by(species, plot_fixed)

# write_csv(LTREB_2021_surviving, path = "~/Documents/R_projects/Grass-Endophyte-Stochastic-Demography/LTREB_2021_surviving.csv")



# Get x and y distances from the a and b corner posts
# we calculate the area based on the perimeter; area = sqrt(s(s-a)(s-b)(s-c)) where s is half the perimeter
# Then we get the height (y value) from the area; h = 2*A/b
# Then we can get the x value as one side of a right triangle with the base and one side (x = sqrt(b^2 - y^2))

# Here I read in the A to B post distances for each plot
dist_between_posts <- read_excel(path = "~/Dropbox/EndodemogData/Fulldataplusmetadata/IndianaPlotPostsABDistances180608.xlsx", sheet = "AB distances") %>% 
  rename(plot_fixed = plot) %>% 
  select(-Notes, -species)

# merge those distances for each plot with our location data for surviving plants
LTREB_ab_distances <- LTREB_2021_surviving %>% 
  left_join(dist_between_posts, by = c("plot_fixed"))

# then calculate the x and y distances based on the circle perimeter and area
# There are a few NA's, probably where the plot distance doesn't make sense with the measured lazer distances
LTREB_xy <- LTREB_ab_distances %>% 
  mutate(triangle_perimeter = dist_a + dist_b + AB_distance,
         s = triangle_perimeter*.5,
         triangle_area = sqrt(s*(s-dist_a)*(s-dist_b)*(s-AB_distance)),
         y = 2*triangle_area/AB_distance,
         x = sqrt(dist_a^2-y^2))


maps_aorigin <- ggplot(data = subset(LTREB_xy, plot_fixed == 6))+
  geom_text(aes(x = x, y = y, label = id)) +
  geom_text(aes(x = 0, y = 0, label = "A"), lwd = 4)+
  geom_text(aes(x = AB_distance, y = 0, label = "B"), lwd = 4) +
  lims(x = c(0,max(subset(LTREB_xy, plot_fixed == 1)$AB_distance)), y = c(0, max(subset(LTREB_xy, plot_fixed == 1)$AB_distance))) +
  theme_classic()
maps_aorigin

maps_aright <- ggplot(data = subset(LTREB_xy, plot_fixed == 1))+
  geom_text(aes(x = AB_distance-x, y = AB_distance-y, label = id)) +
  geom_text(aes(x = AB_distance, y = AB_distance, label = "A"), lwd = 4)+
  geom_text(aes(x = 0, y = AB_distance, label = "B"), lwd = 4)+
  lims(x = c(0,max(subset(LTREB_xy, plot_fixed == 1)$AB_distance)), y = c(0,max(subset(LTREB_xy, plot_fixed == 1)$AB_distance))) +
  theme_classic()
maps_aright

maps_aleft <- ggplot(data = subset(LTREB_xy, plot_index == 1))+
  geom_text(aes(x = AB_distance-x, y = AB_distance-y, label = id)) +
  geom_text(aes(x = AB_distance, y = AB_distance, label = "A"), lwd = 4)+
  geom_text(aes(x = 0, y = AB_distance, label = "B"), lwd = 4)+
  scale_x_reverse(c(max(subset(LTREB_xy, plot_fixed == 1)$AB_distance), 0))+
  lims(y = c(0, max(subset(LTREB_xy, plot_fixed == 1)$AB_distance))) +
  ggtitle(subset(LTREB_xy, plot_index == 1)$species)+
  theme_classic()
maps_aleft


# Saving the left and right plots to pdf.
#Just printing the Elymus maps for now
pdf("2022_Indiana_Maps.pdf")
for (i in sort(unique(subset(LTREB_xy, species == "ELVI" | species == "ELRI" | species == "FESU")$plot_index))){
  maps_aright <- ggplot(data = subset(LTREB_xy,plot_index == i))+
    geom_text(aes(x = AB_distance-x, y = AB_distance-y, label = id)) +
    geom_text(aes(x = AB_distance, y = AB_distance, label = "A"), lwd = 4)+
    geom_text(aes(x = 0, y = AB_distance, label = "B"), lwd = 4)+
    lims(x = c(0,max(subset(LTREB_xy, plot_index ==i)$AB_distance)), y = c(0,max(subset(LTREB_xy, plot_fixed == 1)$AB_distance))) +
    ggtitle(paste(unique(subset(LTREB_xy, plot_index == i)$species), "plot", unique(subset(LTREB_xy, plot_index == i)$plot_fixed), "right"))+
    theme_classic()
  plot(maps_aright)
  
  maps_aleft <- ggplot(data = subset(LTREB_xy, plot_index == i))+
    geom_text(aes(x = AB_distance-x, y = AB_distance-y, label = id)) +
    geom_text(aes(x = AB_distance, y = AB_distance, label = "A"), lwd = 4)+
    geom_text(aes(x = 0, y = AB_distance, label = "B"), lwd = 4)+
    scale_x_reverse(c(max(subset(LTREB_xy, plot_index == i)$AB_distance), 0))+
    lims(y = c(0, max(subset(LTREB_xy, plot_index == 1)$AB_distance))) +
    ggtitle(paste(unique(subset(LTREB_xy, plot_index == i)$species), "plot", unique(subset(LTREB_xy, plot_index == i)$plot_fixed), "left"))+
    theme_classic()
  plot(maps_aleft)
} 

dev.off()


