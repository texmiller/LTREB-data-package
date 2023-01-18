## Authors: Josh and Tom	## Grass endophyte stroma observation counts
## Purpose: Pulling observations of stroma formation on plants out of the raw data files and getting some summaries	
## Last Update: Dec. 20, 2023
######################################################
library(tidyverse)
library(reshape2)
library(lubridate)
library(readxl)
library(patchwork)

########################################################################################################################
###### Stroma observations metadata ---------------------------
########################################################################################################################
## During the annual plot censuses, observations of stroma, indicating sexual reproduction of the Epichloe endophytes, have been recorded, often times as notes in the raw data files
## For the most part, stroma are recorded when they are observed, but it's not a part of the typical demographic data collection. This means that there are some years where observations are recorded in a dedicated column, others where they appear in the notes column and some years where there are no observations recorded in any way. In these cases it is unclear if there were truly no stroma, but we assume that there were none, or at least that they were infrequent enough to not be recorded
## ELRI has the most observations, followed by ELVI and POSY.
## POSY has a few stroma observations, and one census for one plant with the note "FUNGUS!" (unclear whether this is referring to Epichloe stroma, or an endophyte check of seeds, or another fungus).
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
# ELRI
ELRI_2020 <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/Field Data/2020/LTREB_data_2020.xlsx", sheet = "ELRI", .name_repair = "universal") %>% 
  dplyr::select(species, origin,plot,pos,id, observation_year,notes) %>% 
  mutate(stroma_obs = case_when(str_detect(notes, fixed("strom", ignore_case = TRUE)) ~ 1,
                                TRUE ~ 0),
         stroma_count = case_when(stroma_obs == 1 & str_detect(notes, fixed("one", ignore_case = TRUE)) ~ 1,
                          stroma_obs == 1 & str_detect(notes, fixed("two", ignore_case = TRUE)) ~ 2,
                          stroma_obs == 1 ~ as.numeric(str_extract(notes, "\\d+")),
                          TRUE ~ 0)) %>% 
  dplyr::select(plot, pos, id, observation_year, stroma_obs, stroma_count)


ELRI_2021 <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/Field Data/2021/LTREB_data_2021.xlsx", sheet = "ELRI", .name_repair = "universal") %>% 
  dplyr::select(species, origin,plot,pos,id, observation_year,notes) %>% 
  mutate(stroma_obs = case_when(str_detect(notes, fixed("strom", ignore_case = TRUE)) ~ 1,
                                TRUE ~ 0),
         stroma_count = case_when(stroma_obs == 1 & str_detect(notes, fixed("one", ignore_case = TRUE)) ~ 1,
                                  stroma_obs == 1 & str_detect(notes, fixed("two", ignore_case = TRUE)) ~ 2,
                                  stroma_obs == 1 ~ as.numeric(str_extract(notes, "\\d+")),
                                  TRUE ~ 0)) %>% 
  dplyr::select(plot, pos, id, observation_year, stroma_obs, stroma_count)


ELRI_pre <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/ELRI data up to 2016.xlsx", sheet = "ELRI", .name_repair = "universal") %>% 
  pivot_longer(cols = c("Stroma_Tillers09", "Stroma_Tillers10", "Stroma_Tilelrs11", "Notes12", "Notes13", "Notes15", "Notes16"),
               names_to = "stroma_source",
               values_to = "stroma_notes", values_transform = list(stroma_notes = as.character)) %>% 
  mutate(observation_year = case_when(stroma_source == "Stroma_Tillers09" ~ 2009,
                                      stroma_source == "Stroma_Tillers10" ~ 2010,
                                      stroma_source == "Stroma_Tilelrs11" ~ 2011,
                                      stroma_source == "Notes12" ~ 2012,
                                      stroma_source == "Notes13" ~ 2013,
                                      stroma_source == "Notes15" ~ 2015,
                                      stroma_source == "Notes16" ~ 2016)) %>% 
  rename(plot = PLOT, pos = POS, id = TAG) %>% 
  mutate(stroma_count = case_when(str_detect(stroma_notes, fixed("S-T")) | str_detect(stroma_notes, fixed("strom", ignore_case = TRUE)) ~ as.numeric(str_extract(stroma_notes, "\\d+")),
                                  str_detect(stroma_notes, fixed("T/C")) ~ 0,
                                  str_detect(stroma_notes, "\\d+")~as.numeric(stroma_notes),
                                  TRUE ~ 0),
         stroma_obs = case_when(stroma_count>0 ~ 1,
                                TRUE ~ 0)) %>% 
  dplyr::select(plot, pos, id, observation_year, stroma_obs, stroma_count) 
  
ELRI_2017 <- read.csv(file = "~/Dropbox/EndodemogData/Fulldataplusmetadata/endo_demog_long.csv") %>% 
    filter(str_detect(notes, fixed("strom", ignore_case = TRUE))) %>% 
    mutate(stroma_obs = 1, stroma_count = 1) %>% 
    rename(observation_year = year_t) %>% 
    dplyr::select(plot, pos, id, observation_year, stroma_obs, stroma_count)

ELRI_stroma <-rbind(ELRI_pre, ELRI_2017, ELRI_2020, ELRI_2021) %>% 
  mutate(species = "ELRI")
# ELVI

ELVI_stroma<- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/ELVI(IN) originals up to 2016.xlsx", sheet = "ELVI", .name_repair = "universal") %>% 
  rename(plot = PLOT, pos = POS, id = TAG) %>% 
  dplyr::select(plot, pos, id, Stroma_Tillers09) %>% 
  rename(stroma_count = Stroma_Tillers09) %>% 
  mutate(stroma_obs = case_when(stroma_count >0 ~ 1,
                                stroma_count == 0 ~ 0),
        observation_year = 2009,
        species = "ELVI") %>% 
  filter(!is.na(stroma_obs))

# POSY


POSY_new_pre <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/POSYnew complete with 2016 data.xlsx", sheet = "POSY", .name_repair = "universal") %>% 
  rename(id = tag) %>% 
  filter(notes3 == "FUNGUS!") %>% 
  dplyr::select(plot, pos, id) %>% 
  mutate(observation_year = 2010, 
         stroma_count = 1,
         stroma_obs = 1)
  

POSY_old_pre <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/POSYold complete with 2016 data.xlsx", sheet = "POSY(Old)", .name_repair = "universal") %>% 
  rename(plot = PLOT, pos = POS, id = TAG) %>% 
  dplyr::select(plot, pos, id, notes6) %>% 
  mutate(observation_year = 2013,
        stroma_count = as.numeric(case_when(str_detect(notes6, fixed("strom", ignore_case = TRUE)) ~ str_extract(notes6, "\\d+"))),
        stroma_obs = case_when(stroma_count>0 ~ 1)) %>% 
  filter(!is.na(stroma_obs)) %>% 
  dplyr::select(-notes6)
  

POSY_old_r_pre <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/POSYold complete with 2016 data.xlsx", sheet = "POSY(Old)recruits", .name_repair = "universal") %>% 
  rename(plot = Plot, pos = RecruitNo, id = Tag) %>% 
  dplyr::select(plot, pos, id, notes.12) %>% 
  mutate(observation_year = 2012,
         stroma_obs = case_when(str_detect(notes.12, fixed("strom", ignore_case = TRUE)) ~ 1),
         stroma_count = as.numeric(case_when(str_detect(notes.12, fixed("strom", ignore_case = TRUE)) ~ str_extract(notes.12, "\\d+")))) %>% 
  filter(!is.na(stroma_obs)) %>% 
  dplyr::select(-notes.12)
         
POSY_stroma <- rbind(POSY_new_pre, POSY_old_pre, POSY_old_r_pre) %>% 
  mutate(species = "POSY")

stroma_df <- rbind(ELRI_stroma, ELVI_stroma, POSY_stroma) %>% 
  rename(year_t1 = observation_year) %>% 
  mutate(plot_fixed = as.numeric(plot)) %>% 
  filter(stroma_obs>0)

stroma_unique_ids <- stroma_df %>% 
  group_by(species, plot,pos,id) %>% 
  summarize(n(),
            years = toString(year_t1))

stroma_counts_byspecies <- stroma_unique_ids %>% 
  group_by(species) %>% 
  summarize(n())
########################################################################################################################
###### Merging these observations with the compiled dataframe of all the census data---------------------------
########################################################################################################################
LTREB_full <- read_csv(file = "~/Dropbox/EndodemogData/Fulldataplusmetadata/LTREB_full.csv")

stroma_full <- LTREB_full %>% 
  left_join(stroma_df) %>% 
  mutate(stroma_obs = case_when(surv_t1 == 1 & !is.na(stroma_obs) ~ stroma_obs,
                                 surv_t1 == 1 & is.na(stroma_obs) ~ 0),
         stroma_count = case_when(surv_t1 == 1 & !is.na(stroma_count) ~ stroma_count,
                                   surv_t1 == 1 & is.na(stroma_count) ~ 0))

# summarizing how many plants out of the whole experiment have stroma. This is a bit off because some individuals have stroma in multiple years.   
stroma_summary_byyear <- stroma_full %>% 
  group_by(species, year_t1) %>% 
  summarize(plants_obs_with_stromata = sum(stroma_obs, na.rm = T),
            percent_with_stromata_perspeciesyear = sum(stroma_obs, na.rm = T)/length(unique(id))*100,
            percent_with_stromata = sum(stroma_obs, na.rm = T)/length(unique(LTREB_full$id))*100)
            

stroma_summary <- stroma_full %>% 
  group_by(species) %>% 
  summarize(no_obs_with_stromata = sum(stroma_obs, na.rm = T),
            percent_with_stromata_perspecies= sum(stroma_obs, na.rm = T)/length(unique(id))*100,
            percent_with_stromata_overall = sum(stroma_obs, na.rm = T)/length(unique(LTREB_full$id))*100)

stroma_summary <- stroma_full %>% 
  group_by() %>% 
  summarize(no_obs_with_stromata = sum(stroma_obs, na.rm = T),
            percent_with_stromata_perspecies= sum(stroma_obs, na.rm = T)/length(unique(id))*100,
            percent_with_stromata_overall = sum(stroma_obs, na.rm = T)/length(unique(LTREB_full$id))*100)





########################################################################################################################
###### Making a plot of the endophyte checks to quantify the level of plot endophyte status changes ---------------------------
########################################################################################################################

LTREB_status_changes_species <- LTREB_full %>% 
  distinct(species, plot_fixed, id, endo_01, endo_mismatch) %>% 
  group_by(species, endo_01) %>% 
  summarize(Same = sum(endo_mismatch == 0, na.rm = TRUE),
            `Endophyte Lost` = sum(endo_mismatch > 0, na.rm = TRUE),
            `Endophyte Gain` = sum(endo_mismatch <0, na.rm = TRUE)) %>% 
  # mutate(percent_gain = (gain_endo/same)*100, percent_lose = (lose_endo/same)*100) %>% 
  pivot_longer(cols = c(Same, `Endophyte Lost`, `Endophyte Gain`))

endo_check_Eminus_plot <- ggplot(filter(LTREB_status_changes_species, endo_01 == 0 & name !="Endophyte Lost"))+
  geom_bar(aes(x = name, y = value), stat = "identity")+
  facet_wrap(~species, ncol = 1) +
  labs(title = "E- Plots", y = "# of plants scored", x = "")+
  scale_y_continuous(n.breaks = 4)+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+
  theme_classic()+
  theme(axis.text.x = element_text(face = "bold"))
# endo_check_Eminus_plot

endo_check_Eplus_plot <- ggplot(filter(LTREB_status_changes_species, endo_01 == 1 & name !="Endophyte Gain"))+
  geom_bar(aes(x = name, y = value), stat = "identity")+
  facet_wrap(~species, ncol = 1) +
  labs(title = "E+ Plots",y = "# of plants scored", x = "")+
  scale_y_continuous(n.breaks = 4)+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+
  theme_classic()+
  theme(axis.text.x = element_text(face = "bold"))
# endo_check_Eplus_plot

endo_check_plot <- endo_check_Eminus_plot+endo_check_Eplus_plot+
  plot_layout(nrow = 1) + plot_annotation(title = "Endophyte Status Checks", tag_levels = "A")
# endo_check_plot
ggsave(endo_check_plot, filename = "endo_check_plot.png", width = 5, height = 7)
