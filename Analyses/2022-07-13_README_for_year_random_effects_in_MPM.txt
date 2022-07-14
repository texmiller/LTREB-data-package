Lining up the year random effects of vital rate models within the population model\

Meeting Notes - July 13, 2022\
Tom and Josh had a discussion to make sure that we are using the year random effects in the population model correctly for our simulations. \

The overall data from the study that we are using in our analyses includes 14 transition years starting with 2007-2008, and ending with 2020-2021. \

The vital rate models:\

In the vital rate models, we are fitting each predicted value as the value in year t1:\
Surv_t1, size_t1, recruit_t1 are all a function of the relevant year_t values for for plant size or seeds\

Flw_t1, fert_t1, spike_t1 are a function of the relevant year_t1 values for plant size\

We are fitting this by including a year random effect, which is indexed years 1-14. This is a bit confusing because the year index is technically taken from the year_t values, so year_t_index = 1 & year_t = 2007, but this is reflecting the first transition year which is the 2007-2008 transition year. Similarly, the last year index is 14, which is year_t_index = 14 & year = 2020, reflecting the transition year 2020-2021. Our data frame also includes a year_t1_index (2-15), which spans the years 2008-2021. For the purposed of the vital rate models, we use the year_t_index as our index for the year random effects.\

Josh originally set this up while fitting the vital rate models for 2 reasons: (1) it lines up with the main vital rates which are year_t1 values predicted based on the year_t size in the linear predictor, and (2) this was a convenient way to have the year random effects be indexed, rather than doing something like taking the year_t1_index-1.\

It is also worth noting that there isn\'92t any reproductive data for 2007, and there are no/few recruits from before 2009.\

The population model:\

Our population model is made of transition matrices which should represent a transition year including the survival, growth and recruitment from a given transition year, and the fertility (seed production made up from flowering, fertility, and spikelet vital rates) from the preceding year, because reproduction in a given year doesn\'92t contribute to recruitment in that year, but in the next.\

This means that we use the year random effect for a given year for survival etc., and the year random effect for that year minus one for the flowering, fertility and spikelets. Because we don\'92t have flowering data for 2007, and because for sure etc we need we need the preceding year of data the first year in our population model is 2008-2009 transition year. \

This leaves us with 13 transition years in our population model, from 2008-2009 to 2020-2021. This is also probably a reasonable choice because the first transition year was entirely transplanted plants with what seems like unusually high survival, and no recruitment yet.\

We double checked the MPM_functions.R code and it is pulling the year random effects in this way.\

Other Notes:\
- betaendo is the effect of having endophytes, being E+\
- betaorigin its the effect of being a recruit. \
