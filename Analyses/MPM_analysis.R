## Title: Grass endophyte population model with a bayesian framework
## Purpose: Builds matrix population model from vital rate estimates
## and runs stochastic LTRE
## Authors: Joshua and Tom
#############################################################

library(tidyverse)
library(scales)
library(bayesplot)
library(popbio)
library(countreg)
library(actuar)
library(rstan)

quote_bare <- function( ... ){
  substitute( alist(...) ) %>% 
    eval( ) %>% 
    sapply( deparse )
}

#############################################################################################
####### Read in Data and creating size bins------------------
#############################################################################################

source("Analyses/endodemog_data_processing.R")

max_size <- LTREB_full %>% 
  dplyr::select(species,size_t) %>% 
  filter(!is.na(size_t)) %>% 
  group_by(species) %>% 
  summarise(max_size = quantile(size_t,probs=0.975))


#############################################################################################
####### Read in matrix population functions ------------------
#############################################################################################

source("Analyses/MPM_functions.R")


#############################################################################################
####### Read in Stan vital rate model outputs ------------------
#############################################################################################

surv_fit_seedling <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_seedling_surv.rds")
surv_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_surv_woseedling.rds")
grow_fit_seedling <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_seedling_grow.rds")
grow_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_grow_PIG.rds")
flw_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_flw.rds")
fert_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_fert_noplot.rds")
spike_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_spike_year_plot.rds")
seedmean_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/seed_mean.rds")
stos_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_s_to_s.rds") 

surv_par <- rstan::extract(surv_fit, pars =quote_bare(beta0,betasize,betaendo,betaorigin,
                                                                   tau_year, tau_plot))
surv_sdlg_par <- rstan::extract(surv_fit_seedling, pars =quote_bare(beta0,betaendo,
                                                           tau_year, tau_plot))
grow_par <- rstan::extract(grow_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin,
                                                                    tau_year, tau_plot,
                                                       sigma))
grow_sdlg_par <- rstan::extract(grow_fit_seedling, pars = quote_bare(beta0,betaendo,
                                                       tau_year, tau_plot,
                                                       phi))
flow_par <- rstan::extract(flw_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin,
                                                       tau_year, tau_plot))
fert_par <- rstan::extract(fert_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin,
                                                      tau_year)) # plot effects didn't converge
spike_par <- rstan::extract(spike_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin,
                                                                         tau_year, tau_plot))
seed_par <- rstan::extract(seedmean_fit, pars = quote_bare(beta0,betaendo)) #no plot or year effect
recruit_par <- rstan::extract(stos_fit, pars = quote_bare(beta0,betaendo,
                                                          tau_year, tau_plot))

#############################################################################################
####### Run the MPM ------------------
#############################################################################################

# make the list of parameters and calculate mean lambdas
n_draws <- 100
post_draws <- sample.int(1300,size=n_draws) # this is smaller because of the flowering model that was run for a shorter number of iterations
lambda_mean <- array(dim = c(8,2,n_draws))

for(i in 1:length(post_draws)){
  for(e in 1:2){
    for(s in 1:7){
      lambda_mean[s,e,i] <- lambda(bigmatrix(make_params(species=s,
                                                         endo_mean=(e-1),
                                                         endo_var=(e-1),
                                                         original = 1, # should be =1 to represent recruit
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
                                                         recruit_par=recruit_par))$MPMmat)
    }
    lambda_mean[8,e,i] <- mean(lambda_mean[1:7,e,i])
  }
}

saveRDS(lambda_mean, file = "~/Documents/lambda_mean.rds")
# lambda_mean <- read_rds(path = "~/Documents/lambda_mean.rds")

# Mean endophyte difference and quantiles
lambda_mean_diff <- matrix(NA,8,7)
for(s in 1:8){
  lambda_mean_diff[s,1] = mean(lambda_mean[s,2,] - lambda_mean[s,1,])
  lambda_mean_diff[s,2:7] = quantile(lambda_mean[s,2,] - lambda_mean[s,1,],probs=c(0.05,0.125,0.25,0.75,0.875,0.95))
}

## now do variance in lambda 

lambda_hold <- array(dim = c(11,7,2,n_draws))
lambda_var <- array(dim = c(8,2,n_draws))
for(i in 1:length(post_draws)){
  for(e in 1:2){
    for(s in 1:7){
      for(y in 1:11){
        
        lambda_hold[y,s,e,i] <- lambda(bigmatrix(make_params(species=s,
                                                             endo_mean=(e-1),
                                                             endo_var=(e-1),
                                                             original = 1, # should be =1 to represent recruit
                                                             draw=post_draws[i],
                                                             max_size=max_size,
                                                             rfx=T,
                                                             year=y,
                                                             surv_par=surv_par,
                                                             surv_sdlg_par = surv_sdlg_par,
                                                             grow_par=grow_par,
                                                             grow_sdlg_par = grow_sdlg_par,
                                                             flow_par=flow_par,
                                                             fert_par=fert_par,
                                                             spike_par=spike_par,
                                                             seed_par=seed_par,
                                                             recruit_par=recruit_par))$MPMmat)
      }
      lambda_var[s,e,i] <- sd(lambda_hold[,s,e,i])
    }
    lambda_var[8,e,i] <- mean(lambda_var[1:7,e,i])
  }
}
saveRDS(lambda_var, file = "~/Documents/lambda_var.rds")
# lambda_var <- read_rds(path = "~/Documents/lambda_var.rds")

lambda_var_diff <- matrix(NA,8,7)
for(s in 1:8){
  lambda_var_diff[s,1] = mean(lambda_var[s,2,] - lambda_var[s,1,])
  lambda_var_diff[s,2:7] = quantile(lambda_var[s,2,] - lambda_var[s,1,],probs=c(0.05,0.125,0.25,0.75,0.875,0.95))
}

################################################################
##### Plot of mean and variance endo effect on lambda
################################################################

lambda_mean_diff_df <- as_tibble(lambda_mean_diff) %>%  
  rename( "mean" = V1, fifth = V2, twelfthpointfive = V3, twentyfifth = V4, seventyfifth = V5, eightyseventhpointfive = V6, ninetyfifth = V7) %>% 
  mutate(rownames = row.names(.)) %>% 
  mutate(species = case_when(rownames == 1 ~ "Agrostis perennans",
                             rownames == 2 ~ "Elymus virginicus",
                             rownames == 3 ~ "Elymus villosus",
                             rownames == 4 ~ "Festuca subverticillata",
                             rownames == 5 ~ "Lolium arundinaceum",
                             rownames == 6 ~ "Poa alsodes",
                             rownames == 7 ~ "Poa sylvestris",
                             rownames == 8 ~ "Species Mean"))

  
ggplot(data = lambda_mean_diff_df) +
  geom_hline(yintercept = 0, col = "white") + 
  geom_linerange(aes(y = mean, x = species, ymin = 0, ymax = mean), color = "white", lwd =4)+
  geom_point(aes(y = mean, x = species, color = species), lwd = 4) +
  geom_linerange(aes(y = mean, x = species, ymin = twentyfifth, ymax = seventyfifth, color = species), lwd = 2) +
  geom_linerange(aes(y = mean, x = species, ymin = twelfthpointfive, ymax = eightyseventhpointfive, color = species), lwd = 1) +
  geom_linerange(aes(y = mean, x = species, ymin = fifth, ymax = ninetyfifth, color = species)) +
  scale_color_manual(values = c("#ffff99", "#b8e3a0", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#0c2c84", "#A9A9A9")) +
  coord_flip()+
  theme(panel.background = element_rect(fill = "#dbd7c3"),
        panel.grid = element_line(color = NA))
  
# Version with raw posterior draws
dimnames(lambda_mean) <- list(Species = paste0("s",1:8), Endo = paste0("e",1:2), Iteration= paste0("i",1:n_draws))
lambda_mean_cube <- as.tbl_cube(lambda_mean)
lambda_mean_df <- as_tibble(lambda_mean_cube) %>% 
  pivot_wider(names_from = Endo, values_from = lambda_mean) %>% 
  mutate(lambda_diff = e2-e1) %>% 
  mutate(species = case_when(Species == "s1" ~ "Agrostis perennans",
                           Species == "s2" ~ "Elymus virginicus",
                           Species == "s3" ~ "Elymus villosus",
                           Species == "s4" ~ "Festuca subverticillata",
                           Species == "s5" ~ "Lolium arundinaceum",
                           Species == "s6" ~ "Poa alsodes",
                           Species == "s7" ~ "Poa sylvestris",
                           Species == "s8" ~ "Species Mean"))

ggplot(data = lambda_mean_df) +
  geom_hline(yintercept = 0, col = "black") + 
  # geom_linerange(data = lambda_mean_diff_df, aes(x = species, y = mean, ymin = 0, ymax = mean, color = species)) 
  geom_jitter( aes(y = lambda_diff, x = species, color = species), width = .2, alpha = .2) +
  stat_summary(aes(y = lambda_diff, x = species), fun = mean,geom = "point", size = 5) +
  geom_point(data = lambda_mean_diff_df, aes(y = mean, x = species, color = species), lwd = 4) +
  scale_color_manual(values = c("#ffff99", "#b8e3a0", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#0c2c84", "#A9A9A9")) +
  coord_flip() +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = NA))

#############################################################################################
####### Stochastic lambda simulations ------------------
#############################################################################################

# Make a list of year specific transition matrices
lambdaS_out <- array(dim = c(8,4,n_draws))
for(i in 1:length(post_draws)){
  for(s in 1:7){
    eminus_list <- eplus_list <- eplus__mean_only_list <- eplus__var_only_list <- list()
    for(y in 1:11){
      eminus_list[[y]] <- bigmatrix(make_params(species=s,
                                                endo_mean=0,
                                                endo_var=0,
                                                original = 1, # should be =1 to represent recruit
                                                draw=post_draws[i],
                                                max_size=max_size,
                                                rfx=T,
                                                year=y,
                                                surv_par=surv_par,
                                                surv_sdlg_par = surv_sdlg_par,
                                                grow_par=grow_par,
                                                grow_sdlg_par = grow_sdlg_par,
                                                flow_par=flow_par,
                                                fert_par=fert_par,
                                                spike_par=spike_par,
                                                seed_par=seed_par,
                                                recruit_par=recruit_par))$MPMmat
      eplus_list[[y]] <- bigmatrix(make_params(species=s,
                                               endo_mean=1,
                                               endo_var=1,
                                               original = 1, # should be =1 to represent recruit
                                               draw=post_draws[i],
                                               max_size=max_size,
                                               rfx=T,
                                               year=y,surv_par=surv_par,
                                               surv_sdlg_par = surv_sdlg_par,
                                               grow_par=grow_par,
                                               grow_sdlg_par = grow_sdlg_par,
                                               flow_par=flow_par,
                                               fert_par=fert_par,
                                               spike_par=spike_par,
                                               seed_par=seed_par,
                                               recruit_par=recruit_par))$MPMmat
      eplus__mean_only_list[[y]] <- bigmatrix(make_params(species=s,
                                                          endo_mean=1,
                                                          endo_var=0,
                                                          original = 1, # should be =1 to represent recruit
                                                          draw=post_draws[i],
                                                          max_size=max_size,
                                                          rfx=T,
                                                          year=y,surv_par=surv_par,
                                                          surv_sdlg_par = surv_sdlg_par,
                                                          grow_par=grow_par,
                                                          grow_sdlg_par = grow_sdlg_par,
                                                          flow_par=flow_par,
                                                          fert_par=fert_par,
                                                          spike_par=spike_par,
                                                          seed_par=seed_par,
                                                          recruit_par=recruit_par))$MPMmat
      eplus__var_only_list[[y]] <- bigmatrix(make_params(species=s,
                                                         endo_mean=0,
                                                         endo_var=1,
                                                         original = 1, # should be =1 to represent recruit
                                                         draw=post_draws[i],
                                                         max_size=max_size,
                                                         rfx=T,
                                                         year=y,surv_par=surv_par,
                                                         surv_sdlg_par = surv_sdlg_par,
                                                         grow_par=grow_par,
                                                         grow_sdlg_par = grow_sdlg_par,
                                                         flow_par=flow_par,
                                                         fert_par=fert_par,
                                                         spike_par=spike_par,
                                                         seed_par=seed_par,
                                                         recruit_par=recruit_par))$MPMmat
    }
    lambdaS_out[s,1,i] <- lambdaSim(eminus_list)$lambdaS
    lambdaS_out[s,2,i] <- lambdaSim(eplus__mean_only_list)$lambdaS
    lambdaS_out[s,3,i] <- lambdaSim(eplus__var_only_list)$lambdaS
    lambdaS_out[s,4,i] <- lambdaSim(eplus_list)$lambdaS
  }
  lambdaS_out[8,1,i] <- mean(lambdaS_out[1:7,1,i])
  lambdaS_out[8,2,i] <- mean(lambdaS_out[1:7,2,i])
  lambdaS_out[8,3,i] <- mean(lambdaS_out[1:7,3,i])
  lambdaS_out[8,4,i] <- mean(lambdaS_out[1:7,4,i])
}

saveRDS(lambdaS_out, file = "~/Documents/lambdaS_out.rds")
# lambdaS_out <- read_rds(path = "~/Documents/lambdaS_out.rds")


lambdaS_diff <- lambdaS_diff_mean_only <- lambdaS_diff_var_only <- matrix(NA,8,7)
for(s in 1:8){
  lambdaS_diff[s,1] = mean(lambdaS_out[s,4,] - lambdaS_out[s,1,], na.rm = T)
  lambdaS_diff[s,2:7] = quantile(lambdaS_out[s,4,] - lambdaS_out[s,1,],probs=c(0.05,0.125,0.25,0.75,0.875,0.95), na.rm = T)
  
  lambdaS_diff_mean_only[s,1] = mean(lambdaS_out[s,2,] - lambdaS_out[s,1,], na.rm = T)
  lambdaS_diff_mean_only[s,2:7] = quantile(lambdaS_out[s,2,] - lambdaS_out[s,1,],probs=c(0.05,0.125,0.25,0.75,0.875,0.95), na.rm = T)
  
  lambdaS_diff_var_only[s,1] = mean(lambdaS_out[s,3,] - lambdaS_out[s,1,], na.rm = T)
  lambdaS_diff_var_only[s,2:7] = quantile(lambdaS_out[s,3,] - lambdaS_out[s,1,],probs=c(0.05,0.125,0.25,0.75,0.875,0.95), na.rm = T)
}

lambdaS_effect <- rowMeans(lambdaS_out[8,2:4,]-lambdaS_out[8,1,], na.rm = T)
barplot(lambdaS_effect,col=rainbow(3),
        names.arg = c("Mean effect","Variance effect","Total effect"))
################################################################
##### Plot of stochastic lambda contributions
################################################################
# I need to check the indexing of the variance endophyte effect because it seems like it is backwards... (total effect is less than mean + var)
# I am relabeling them here for ESA, but I need to fix this in the vital rate models.
lambdaS_diff_df <- as_tibble(lambdaS_diff) %>%  
  rename( "mean" = V1, fifth = V2, twelfthpointfive = V3, twentyfifth = V4, seventyfifth = V5, eightyseventhpointfive = V6, ninetyfifth = V7) %>% 
  mutate(rownames = row.names(.)) %>% 
  mutate(species = case_when(rownames == 1 ~ "Agrostis perennans",
                             rownames == 2 ~ "Elymus virginicus",
                             rownames == 3 ~ "Elymus villosus",
                             rownames == 4 ~ "Festuca subverticillata",
                             rownames == 5 ~ "Lolium arundinaceum",
                             rownames == 6 ~ "Poa alsodes",
                             rownames == 7 ~ "Poa sylvestris",
                             rownames == 8 ~ "Species Mean")) %>% 
  mutate(contribution = "Mean")
lambdaS_diff_mean_only_df <- as_tibble(lambdaS_diff_mean_only) %>%  
  rename( "mean" = V1, fifth = V2, twelfthpointfive = V3, twentyfifth = V4, seventyfifth = V5, eightyseventhpointfive = V6, ninetyfifth = V7) %>% 
  mutate(rownames = row.names(.)) %>% 
  mutate(species = case_when(rownames == 1 ~ "Agrostis perennans",
                             rownames == 2 ~ "Elymus virginicus",
                             rownames == 3 ~ "Elymus villosus",
                             rownames == 4 ~ "Festuca subverticillata",
                             rownames == 5 ~ "Lolium arundinaceum",
                             rownames == 6 ~ "Poa alsodes",
                             rownames == 7 ~ "Poa sylvestris",
                             rownames == 8 ~ "Species Mean")) %>% 
  mutate(contribution = "Total")
lambdaS_diff_var_only_df <- as_tibble(lambdaS_diff_var_only) %>%  
  rename( "mean" = V1, fifth = V2, twelfthpointfive = V3, twentyfifth = V4, seventyfifth = V5, eightyseventhpointfive = V6, ninetyfifth = V7) %>% 
  mutate(rownames = row.names(.)) %>% 
  mutate(species = case_when(rownames == 1 ~ "Agrostis perennans",
                             rownames == 2 ~ "Elymus virginicus",
                             rownames == 3 ~ "Elymus villosus",
                             rownames == 4 ~ "Festuca subverticillata",
                             rownames == 5 ~ "Lolium arundinaceum",
                             rownames == 6 ~ "Poa alsodes",
                             rownames == 7 ~ "Poa sylvestris",
                             rownames == 8 ~ "Species Mean")) %>% 
  mutate(contribution = "Var")

lambdaS_diff_all_df <- bind_rows(lambdaS_diff_df, lambdaS_diff_mean_only_df, lambdaS_diff_var_only_df)
lambdaS_diff_all_df$contribution <- factor(lambdaS_diff_all_df$contribution, levels = c("Var", "Mean", "Total"))
x_labels <- c(expression(paste("Var(", lambda, ")")), expression(bar(lambda)), "Total")
ggplot(data = lambdaS_diff_all_df) +
  geom_hline(yintercept = 0, col = "black") + 
  geom_linerange(aes(y = mean, x = contribution, ymin = twentyfifth, ymax = seventyfifth, color = species), lwd = 2) +
  geom_linerange(aes(y = mean, x = contribution, ymin = twelfthpointfive, ymax = eightyseventhpointfive, color = species), lwd = 1) +
  geom_linerange(aes(y = mean, x = contribution, ymin = fifth, ymax = ninetyfifth, color = species)) +
  geom_point(aes(y = mean, x = contribution, fill = species, pch = contribution), lwd = 4) + 
  facet_wrap(~species, nrow = 2, scales = "free") + coord_flip() + 
  scale_shape_manual(values = c(21, 22, 23)) +
  scale_fill_manual(values = c("#59A1BC", "#4E816D", "#EFAD3A", "#E3AF7B", "#9E78A1", "#E04D55", "#9D5251", "#5C5C5D")) +
  scale_color_manual(values = c("#59A1BC", "#4E816D", "#EFAD3A", "#E3AF7B", "#9E78A1", "#E04D55", "#9D5251", "#5C5C5D")) +
  xlab("") + scale_x_discrete(labels = x_labels) +
  ylab(expression(paste("Endophyte effect on", " ", lambda["s"]))) +
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.x = element_line(size = .5, colour = "black"),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = rel(1.5), face = "bold"),
        axis.text.x = element_text(size = rel(1.5), face = "bold"),
        axis.title = element_text(size = rel(1.5)),
        strip.background = element_blank(),
        legend.position="none")

# Just the mean effect
sppmean_stochplot <- ggplot(data = subset(lambdaS_diff_all_df, species == "Species Mean")) +
  geom_hline(yintercept = 0, col = "black") + 
  geom_linerange(aes(y = mean, x = contribution, ymin = twentyfifth, ymax = seventyfifth, color = species), lwd = 2) +
  geom_linerange(aes(y = mean, x = contribution, ymin = twelfthpointfive, ymax = eightyseventhpointfive, color = species), lwd = 1) +
  geom_linerange(aes(y = mean, x = contribution, ymin = fifth, ymax = ninetyfifth, color = species)) +
  geom_point(aes(y = mean, x = contribution, fill = species, pch = contribution), lwd = 4) + 
  facet_wrap(~species, nrow = 2, scales = "free") + coord_flip() + 
  scale_shape_manual(values = c(21, 22, 23)) +
  scale_fill_manual(values = c( "#5C5C5D")) +
  scale_color_manual(values = c("#5C5C5D")) +
  xlab("") + scale_x_discrete(labels = x_labels) +
  ylab(expression(paste("Endophyte effect on", " ", lambda["s"]))) +
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.x = element_line(size = .5, colour = "black"),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = rel(2), face = "bold"),
        axis.text.x = element_text(size = rel(2), face = "bold"),
        axis.title = element_text(size = rel(3)),
        strip.text = element_blank(),
        strip.background = element_blank(),
        legend.position="none")
sppmean_stochplot
ggsave(sppmean_stochplot, filename = "sppmean_stochplot.tiff",  width = 6, height = 4, bg = "transparent")


################################################################
# Testing out the parameter vectors for mean and variance
for(s in 1:7){
  eminus_listtest <- eplus_listtest <- eplus__meanonly_listtest <- eplus__varonly_listtest <- list()
  for(y in 1:11){
    eminus_listtest[[y]] <- make_params(species=s,
                                       endo_mean=0,
                                       endo_var=0,
                                       original = 1, # should be =1 to represent recruit
                                       draw=post_draws[100],
                                       max_size=max_size,
                                       rfx=T,
                                       year=y,surv_par=surv_par,
                                       surv_sdlg_par = surv_sdlg_par,
                                       grow_par=grow_par,
                                       grow_sdlg_par = grow_sdlg_par,
                                       flow_par=flow_par,
                                       fert_par=fert_par,
                                       spike_par=spike_par,
                                       seed_par=seed_par,
                                       recruit_par=recruit_par)
    eplus_listtest[[y]] <- make_params(species=s,
            endo_mean=1,
            endo_var=1,
            original = 1, # should be =1 to represent recruit
            draw=post_draws[100],
            max_size=max_size,
            rfx=T,
            year=y,surv_par=surv_par,
            surv_sdlg_par = surv_sdlg_par,
            grow_par=grow_par,
            grow_sdlg_par = grow_sdlg_par,
            flow_par=flow_par,
            fert_par=fert_par,
            spike_par=spike_par,
            seed_par=seed_par,
            recruit_par=recruit_par)
    eplus__meanonly_listtest[[y]] <- make_params(species=s,
                                       endo_mean=1,
                                       endo_var=0,
                                       original = 1, # should be =1 to represent recruit
                                       draw=post_draws[100],
                                       max_size=max_size,
                                       rfx=T,
                                       year=y,surv_par=surv_par,
                                       surv_sdlg_par = surv_sdlg_par,
                                       grow_par=grow_par,
                                       grow_sdlg_par = grow_sdlg_par,
                                       flow_par=flow_par,
                                       fert_par=fert_par,
                                       spike_par=spike_par,
                                       seed_par=seed_par,
                                       recruit_par=recruit_par)
    eplus__varonly_listtest[[y]] <- make_params(species=s,
                                                 endo_mean=0,
                                                 endo_var=1,
                                                 original = 1, # should be =1 to represent recruit
                                                 draw=post_draws[100],
                                                 max_size=max_size,
                                                 rfx=T,
                                                 year=y,surv_par=surv_par,
                                                 surv_sdlg_par = surv_sdlg_par,
                                                 grow_par=grow_par,
                                                 grow_sdlg_par = grow_sdlg_par,
                                                 flow_par=flow_par,
                                                 fert_par=fert_par,
                                                 spike_par=spike_par,
                                                 seed_par=seed_par,
                                                 recruit_par=recruit_par)

}}

mean_surv_coef <- lapply(rstan::extract(surv_fit, pars =quote_bare(beta0,betasize,betaendo,betaorigin,
                                                                    tau_year, tau_plot))
                         , colMeans)

mean_surv_seedling_coef <- lapply(rstan::extract(surv_fit_seedling, pars = quote_bare(beta0,betaendo,
                                                                             tau_year, tau_plot))
                                  , colMeans)

mean_grow_coef <- lapply(rstan::extract(grow_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin,
                                                                    tau_year, tau_plot))
                         , colMeans)

mean_grow_seedling_coef <- lapply(rstan::extract(grow_fit_seedling, pars = quote_bare(beta0,betaendo,
                                                                             tau_year, tau_plot))
                                  , colMeans)


mean_flw_coef <- lapply(rstan::extract(flw_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin,
                                                                    tau_year, tau_plot))
                         , colMeans)

mean_fert_coef <- lapply(rstan::extract(fert_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin,
                                                                    tau_year)) # No plot effect
                         , colMeans)

mean_spikelet_coef <- lapply(rstan::extract(spike_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin,
                                                                        tau_year, tau_plot))
                             , colMeans)

mean_s_to_s_coef <- lapply(rstan::extract(stos_fit, pars = quote_bare(beta0,betaendo,
                                                                      tau_year, tau_plot))
                           , colMeans)

mean_seedmean_coef <- lapply(rstan::extract(seedmean_fit, pars = quote_bare(beta0,betaendo)) #no plot or year effect
                              , colMeans)





