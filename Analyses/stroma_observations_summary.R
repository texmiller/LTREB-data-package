## Authors: Josh and Tom	## Grass endophyte stroma observation counts
## Purpose: Pulling observations of stroma formation on plants out of the raw data files and getting some summaries	
## Last Update: Dec. 20, 2023
######################################################
library(tidyverse)
library(reshape2)
library(lubridate)
library(readxl)

########################################################################################################################
###### Stroma observations metadata ---------------------------
########################################################################################################################
## During the annual plot censuses, observations of stroma, indicating sexual reproduction of the Epichloe endophytes, have been recorded, often times as notes in the raw data files
## For the most part, stroma are recorded when they are observed, but it's not a part of the typical demographic data collection. This means that there are some years where observations are recorded in a dedicated column, others where they appear in the notes column and some years where there are no observations recorded in any way. In these cases it is unclear if there were truly no stroma, but we assume that there were none, or at least that they were infrequent enough to not be recorded
## ELRI has the most observations, followed by ELVI.
## POSY has 1 stroma observation, and one census for one plant with the note "FUNGUS!" and it unclear whether this is referring to Epichloe stroma, or an endophyte check of seeds, or another fungus.
## There are no observations recorded for AGPE, FESU, LOAR, POAL
## For ELRI, there are certain years which seem to have many stroma observations and others with few/none.

## The more contemporary data are stored separately for each year
## There are no observations for stromata for any species in the 2019 data
## 2020 has ELRI stromata observations
## 2021 has ELRI stromata observations
## 2022 has no stromata observations for any species, and that census data hasn't been included in the full cleaned dataset yet anyways
########################################################################################################################
###### Reading in raw excel files which have the detailed notes for each species---------------------------
########################################################################################################################

ELRI_2020 <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/Field Data/2020/LTREB_data_2020.xlsx", sheet = "ELRI", .name_repair = "universal") %>% 
  dplyr::select(species, origin,plot,pos,id, observation_year,notes) %>% 
  mutate(stroma_obs = case_when(str_detect(notes, fixed("strom", ignore_case = TRUE)) ~ 1,
                                TRUE ~ 0),
         stroma_count = case_when(stroma_obs == 1 & str_detect(notes, fixed("one", ignore_case = TRUE)) ~ 1,
                          stroma_obs == 1 & str_detect(notes, fixed("two", ignore_case = TRUE)) ~ 2,
                          stroma_obs == 1 ~ as.numeric(str_extract(notes, "\\d+")),
                          TRUE ~ 0))

ELRI_2021 <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/Field Data/2021/LTREB_data_2021.xlsx", sheet = "ELRI", .name_repair = "universal") %>% 
  dplyr::select(species, origin,plot,pos,id, observation_year,notes) %>% 
  mutate(stroma_obs = case_when(str_detect(notes, fixed("strom", ignore_case = TRUE)) ~ 1,
                                TRUE ~ 0),
         stroma_count = case_when(stroma_obs == 1 & str_detect(notes, fixed("one", ignore_case = TRUE)) ~ 1,
                                  stroma_obs == 1 & str_detect(notes, fixed("two", ignore_case = TRUE)) ~ 2,
                                  stroma_obs == 1 ~ as.numeric(str_extract(notes, "\\d+")),
                                  TRUE ~ 0))
