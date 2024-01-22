library(tidyverse)
## read in NOAA data from nearest weathr stations to LDW and SFAEF
indiana<-read.csv("C:/Users/tm9/Dropbox/EndodemogData/climate/IN_weather_station_NOAA_20072023.csv")
indiana$site="LDW"
texas<-read.csv("C:/Users/tm9/Dropbox/EndodemogData/climate/TX_weather_station_NOAA_20072023.csv")
texas$site="SFA"

summary(texas$PRCP)

bind_rows(indiana,texas) %>% 
  mutate(year = substr(DATE,1,4),
         month = substr(DATE,6,7),
         day = substr(DATE,9,10)) %>% 
  group_by(site,year) %>% 
  summarise(total_precip = sum(PRCP,na.rm=T),
            precip_na = sum(is.na(PRCP)),
            mean_temp = mean(TOBS,na.rm=T),
            temp_na = sum(is.na(TOBS))) -> weather

weather %>% 
  ggplot()+
  geom_point(aes(x=year,y=total_precip))+
  facet_grid(~site)

weather %>% 
  ggplot()+
  geom_point(aes(x=year,y=mean_temp))+
  facet_grid(~site)

str(indiana)
indiana$DATE[1]
substr(indiana$DATE[1],1,4)
substr(indiana$DATE[1],6,7)
substr(indiana$DATE[1],9,10)
