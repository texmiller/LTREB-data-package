library(rstan)
library(bayesplot)
library(MASS)
library(tidyverse)

invlogit <- function(x){exp(x)/(1+exp(x))}
quote_bare <- function( ... ){
  substitute( alist(...) ) %>% 
    eval( ) %>% 
    sapply( deparse )
}

## one random normal deviate and one random bernoulli deviate
reps <- 100
n_blocks <- 15
y1_mean <- 86.2
y1_block_sigma <- 25
y2_mean <- 0
y2_block_sigma <- 5
y1y2_cor <- -0.25
block_means <- mvrnorm(n=n_blocks,
                       mu=c(y1_mean,y2_mean),
                       Sigma=matrix(c(y1_block_sigma^2,y1_block_sigma*y2_block_sigma*y1y2_cor,
                                      y1_block_sigma*y2_block_sigma*y1y2_cor,y2_block_sigma^2),2,2))
## now treat y1 and y2 as being sampled independently
resid_sigma1 <- 6.6
df1 <- data.frame(block=rep(1:n_blocks,each=reps),y1=rnorm(n_blocks*reps,rep(block_means[,1],each=reps),resid_sigma1))
df2 <- data.frame(block=rep(1:n_blocks,each=reps),y2=rbinom(n_blocks*reps,size=1,prob=invlogit(rep(block_means[,2],each=reps))))

dat <- list(
  y1=df1$y1,
  N1=nrow(df1),
  n_blocks1=length(unique(df1$block)),
  block1=df1$block,
  y2=df2$y2,
  N2=nrow(df2),
  n_blocks2=length(unique(df2$block)),
  block2=df2$block
)

combo_model = stan_model("Analyses/combo_model.stan")
fit_combo <- sampling(combo_model,data = dat,
                      chains = 2, iter = 5000, warmup = 1000)
block_samples <- rstan::extract(fit_combo, pars =quote_bare(eps1,eps2))

block_corr<-c()
for(i in sample(1:nrow(block_samples$eps1),size=1000)){
  block_corr[i]<-cor(block_samples$eps1[i,],block_samples$eps2[i,])
}
hist(block_corr);abline(v=y1y2_cor,col="red");abline(v=cor(block_means[,1],block_means[,2]))
