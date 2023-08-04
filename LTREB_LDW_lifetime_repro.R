library(tidyverse)
library(ggExtra)
library(patchwork)
library(RColorBrewer)

## read in LTREB data
LTREB<-read.csv("LDW_LTREB_20072022.csv") %>% 
  ## if there is size data but infs is NA then set infs to zero
  mutate(infs_t=ifelse(!is.na(size_t) & is.na(flw_count_t),0,flw_count_t),
         infs_t1=ifelse(!is.na(size_t1) & is.na(flw_count_t1),0,flw_count_t1)) %>% 
  filter(origin_01==1) %>% ## filter out recruits only (true known ages)
  mutate(age = year_t-birth) %>% #calculate age 
  mutate(age_lump = ifelse(age<=4,age,5)) %>% 
  mutate(avg_spike_t= mean(c_across(spike_a_t:spike_d_t), na.rm = T),
         avg_spike_t1= mean(c_across(spike_a_t1:spike_d_t1), na.rm = T)) %>% 
  select(endo_01, id, plot, age, age_lump, surv_t1,
         year_t, year_t1, flw_count_t, flw_count_t1,avg_spike_t,avg_spike_t1,
         species) %>% 
  filter(age>=0)

## where spike data is missing, fill with species mean
LTREB %>% 
  group_by(species) %>% 
  summarise(mean_spikes_t=mean(avg_spike_t,na.rm=T),
            mean_spikes_t1=mean(avg_spike_t1,na.rm=T)) -> mean_spikelets

LTREB %>% 
  ##filter out anything still alive in 2022 (bc we dont have lifetime reproduction)
  mutate(alive22 = ifelse(year_t1==2022 & surv_t1==1,1,0)) %>% 
  group_by(species,id) %>% 
  summarise(alive22=sum(alive22),
            lifespan = max(age),
            total_inf = sum(flw_count_t,na.rm=T)) %>% 
  filter(alive22!=1) -> lifetimes

spp_list<-unique(lifetimes$species)
spp_cols<-brewer.pal(n=7,name="Dark2")
spp_names<-c("Agrostis perennans","Elymus villosus","Elymus virginicus",
            "Festuca subverticillata","Lolium aryndinaceum",
            "Poa alsodes","Poa sylvestris")
par(mfrow=c(2,4))
for(i in 1:length(spp_list)){
  dat<-lifetimes %>% filter(species==spp_list[i])
  plot(jitter(dat$lifespan),jitter(dat$total_inf),main=spp_names[i],
       xlab="Lifespan",ylab="Lifetime inflorescences")
}

plots<-list()
for(i in 1:length(spp_list)){
  dat<-lifetimes %>% filter(species==spp_list[i])
  holdit<-ggplot(dat) +
    geom_jitter(aes(x = lifespan, y = total_inf), alpha = 0.25, shape = 16, color = spp_cols[i]) +  
    labs(x = "Lifespan (years)", y = "Lifetime inflorescences")+
    theme_classic()+ggtitle(spp_names[i])+theme(plot.title = element_text(size = 9, face = "bold"))
  plots[[i]]<-ggMarginal(holdit, type="histogram",fill=spp_cols[i],binwidth=1,
                         xparams = list(size=0), yparams = list(size=0))
}

patchwork::wrap_elements(plots[[1]]) + patchwork::wrap_elements(plots[[2]]) +
  patchwork::wrap_elements(plots[[3]]) + patchwork::wrap_elements(plots[[4]]) +
  patchwork::wrap_elements(plots[[5]]) + patchwork::wrap_elements(plots[[6]]) +
  patchwork::wrap_elements(plots[[7]]) + patchwork::plot_layout(nrow=2)


ggplot(lifetimes) +
  geom_point(aes(x = lifespan, y = total_inf), alpha = 0.6, shape = 16) +  
  labs(x = "Lifespan (years)", y = "Lifetime inflorescences") +
  geom_density()

## calculate Pr(Lifespan|LRO)
