# Endo_Stoch_Demo
This repository houses analysis and manuscript for the Stochastic  Demography of Fungal Endophyte - Grass Populations
### Repository Authors: 
Josh Fowler and Tom Miller
### Project Co-Authors: 
Jenn Rudgers, Ken Whitney, and Shaun Ziegler

## README file last updated: 
Jan 11, 2021

### Project Overview:
Fungal endophytes are widespread symbionts of grasses that have been shown to provide a variety of context dependent benefits under environmental stress such as drought or salinity tolerance. Context-dependence may make interactions seem unpredictable as environmental conditions vary between years, but it provides a distinct mechanism by which symbionts may act as mutualists beyond influencing mean population growth rates. We are quantifying the relative importance of mean and variance effects of microbial symbionts using long-term demographic data from experimental plots in Indiana at Lilly-Dickey Woods. 

The experiment, started in 2007, comprises 10-18 plots for each of 7 species of grass hosts. Half the plots were planted with 20 endophyte-infected plants, and half with 20 endophyte-free plants. The plots are censused annually for plant survival, plant size (measured as number of tillers) and reproduction (measured as number of flowering tillers and measurements of seed and spikelets).

We estimate individual level vital rates from this experimental data. From these vital rate estimates, we build stochastic Integral Projection models to make population projections which allow us to assess the contributions of endophyte partnership on both the mean and variance components of host fitness. 

## Repository Folder Description:
This repository is set up with three folders:
### Analyses 
holds scripts to cobble data and run analyses
File Name  | Description
------------- | -------------
MPM_analysis.R  | Script to analyse matrix population model with and without endo effects on mean and variance. This script also performs stochastic simulations for life table response experiment.  Also includes a few plots at the end visualizing growth rates.
MPM_functions.R | Includes functions used in MPM_analysis.R
endo_spp_grow_fert.stan | Stan model for growth and no. of inflorescences vital rate model for multiple species with negative binomial distribution with year and plot random effects
endo_spp_grow_fert_PIG.stan | Stan model for growth and no. of inflorescences vital rate model for multiple species with Poisson Inverse Gaussian distribution with year and plot random effects
endo_spp_grow_fert_noplot.stan | Stan model for growth and no. of inflorescences vital rate model for multiple species with negative binomial distribution with only year random effects
endo_spp_s_to_s.stan | Stan model for seed to seedling transition  (germination) vital rate model for multiple species with Binomial distribution with year and plot random effects
endo_spp_seed_mean.stan | Stan model for mean seed per spikelet estimates for multiple species with normal distribution and no random effects
endo_spp_spike.stan | Stan model for spikelets per inflorescence vital rate model for multiple species with Poisson distribution with year and plot random effects (This model has no endophyte effect on variance)
endo_spp_surv_flw.stan | Stan model for survival and flowering status vital rate models for multiple species with Bernoulli distribution with year and plot random effects.
endodemog_data_processing.R | Script that cleans legacy experimental data (2007-2018) and merges this legacy with ongoing field data (2019-2020). Data is stored in the Dropbox folder "EndodemogData". Legacy data manipulation involves pulling reproductive data out of the spreadsheets and merging this with size and survival data. For recent data, these are cleaned and also stored in the Dropbox sub-folder "Field Data", which include data as collected in the field, as well as cleaned data that correct some missing fields and assign tag not found status to recruits or previously found plants if possible.
seed_means.R | Script to run seed means model and visualize model diagnostics
seed_to_seedling.R | Script to run germination model and visualize model diagnostics
seedling_grow.stan | Stan model for first year seedling growth vital rate model for multiple species with negative binomial distribution with year and plot random effects. Our matrix model assumes a reproductive delay where first year plants of 1 tiller size do not reproduce.
seedling_surv.stan | Stan model for first year seedling survival vital rate model for multiple species with Bernoulli distribution with year and plot random effects. Our matrix model assumes a reproductive delay where first year plants of 1 tiller size do not reproduce.
stochastic_lambda_analysis.R | Test script to run stochastic population growth simulations and life table response experiment, replaced by MPM_analysis.R 
vital_rate_analysis.R | Script to run survival, growth, and fertility vital rate models, and visualize model diagnostics


### Manuscript 
holds drafts and figures of manuscript

File Name  | Description
------------- | -------------
endo_stoch_demo_MS.Rnw | Sweave document for writing the manuscript
endo_stoch_demo_MS.pdf | PDF output of the manuscript sweave document
endo_stoch_demo_MS.bbl | bibtex file for citations
Other misc. files | Various style and log files for latex, based around formatting for PNAS


### Meeting Notes 
holds weekly-ish Josh and Tom's meeting notes, including weekly goals and discussions beyond this project

File Name  | Description
------------- | -------------
MeetingUpdates.tex | Latex document for recording weekly notes and goals
MeetingUpdates.pdf | PDF document output of MeetingUpdates.tex
Other misc. files | Various style and log files for latex



