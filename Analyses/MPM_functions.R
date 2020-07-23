## Title: Grass endophyte population model with a bayesian framework
## Purpose: functions for building matrix population model from vital rate estimates 
## Authors: Joshua and Tom
#############################################################

invlogit<-function(x){exp(x)/(1+exp(x))}

# Parameter assembly function ---------------------------------------------
make_params <- function(species,endo_mean,endo_var,draw,original=0,rfx=F,year=NULL,max_size,
                        surv_par,surv_sdlg_par,grow_par,grow_sdlg_par,flow_par,fert_par,spike_par,seed_par,recruit_par){
  
  if(rfx==F){rfx_surv <- rfx_surv_sdlg <- rfx_grow <- rfx_grow_sdlg <- rfx_flow <- rfx_fert <- rfx_spike <- rfx_rct <-  0}
  if(rfx==T){
    ## timing and survival and growth (size_t / y_t1) is meant to line up with reproduction (size_t1 / y_t1)
    rfx_surv <- surv_par$tau_year[draw,species,(endo_var+1),(year+1)]; 
    rfx_surv_sdlg <-surv_sdlg_par$tau_year[draw,species,(endo_var+1),(year+1)];
    rfx_grow <- grow_par$tau_year[draw,species,(endo_var+1),(year+1)];
    rfx_grow_sdlg <- grow_sdlg_par$tau_year[draw,species,(endo_var+1),(year+1)];
    rfx_flow <- flow_par$tau_year[draw,species,(endo_var+1),year];
    rfx_fert <- fert_par$tau_year[draw,species,(endo_var+1),year]; 
    rfx_spike <- spike_par$tau_year[draw,year]; #no endo effects or species differences
    rfx_rct <- recruit_par$tau_year[draw,species, (endo_var+1),year];
  }
  
  params <- c()
  #survival
  params$surv_int <- surv_par$beta0[draw,species] + 
    endo_mean * surv_par$betaendo[draw,species] + 
    original * surv_par$betaorigin[draw,species] + rfx_surv
  params$surv_slope <- surv_par$betasize[draw,species]
  # seedling survival
  params$surv_sdlg_int <- surv_sdlg_par$beta0[draw,species] + 
    endo_mean * surv_sdlg_par$betaendo[draw,species] + 
    rfx_surv_sdlg
  #growth
  params$grow_int <- grow_par$beta0[draw,species] + 
    endo_mean * grow_par$betaendo[draw,species] + 
    original * grow_par$betaorigin[draw,species] + rfx_grow
  params$grow_slope <- grow_par$betasize[draw,species]  
  params$grow_sigma <- grow_par$sigma[draw] 
  # seedling growth
  params$grow_sdlg_int <- grow_sdlg_par$beta0[draw,species] + 
    endo_mean * grow_sdlg_par$betaendo[draw,species] + 
    rfx_grow_sdlg
  params$grow_sdlg_phi <- grow_sdlg_par$phi[draw] 
  
  #flowering
  params$flow_int <- flow_par$beta0[draw,species] + 
    endo_mean * flow_par$betaendo[draw,species] + 
    original * flow_par$betaorigin[draw,species] + rfx_flow
  params$flow_slope <- flow_par$betasize[draw,species]  
  #fertility
  params$fert_int <- fert_par$beta0[draw,species] +
   endo_mean * fert_par$betaendo[draw,species] +
   original * fert_par$betaorigin[draw,species] + rfx_fert
  params$fert_slope <- fert_par$betasize[draw,species]

  #spikelets
  params$spike_int <- spike_par$beta0[draw,species]  +
    endo_mean * spike_par$betaendo[draw,species] +
    original * spike_par$betaorigin[draw,species] + rfx_spike
  params$spike_slope <- spike_par$betasize[draw,species]  
  #seeds per spikelet
  params$seeds_per_spike <- seed_par$beta0[draw,species] + 
    endo_mean * seed_par$betaendo[draw,species]
  #recruits per seed
  params$recruits_per_seed <- recruit_par$beta0[draw,species] + 
    endo_mean * recruit_par$betaendo[draw,species] + rfx_rct
  #tack on max size
  params$max_size <- max_size$max_size[species]
  
  return(params)
}
# Vital rate functions ----------------------------------------------------
sx<-function(x,params){
  invlogit(params$surv_int + params$surv_slope*log(x))
}
sx_sdlg <- function(params){
  invlogit(params$surv_sdlg_int)
}

gxy <- function(x,y,params){
  grow_mean <- params$grow_int + params$grow_slope*log(x)
  
  grow<-dpoisinvgauss(x=y,mean=exp(grow_mean),shape=(exp(grow_mean)*params$grow_sigma))
  grow<-ifelse(is.nan(grow) | is.infinite(grow),0,grow)
  
  truncLower<-dpoisinvgauss(x=0,mean=exp(grow_mean), shape=(exp(grow_mean)*params$grow_sigma))
  # truncLower<-sum(ifelse(is.nan(truncLower) | is.infinite(truncLower),0,truncLower))
  
  truncUpper<-sum(dpoisinvgauss(x=params$max_size:10000,mean=exp(grow_mean),shape=(exp(grow_mean)*params$grow_sigma)))
  # truncUpper<-sum(ifelse(is.nan(truncUpper) | is.infinite(truncUpper),0,truncUpper))
  return(grow/(1-(truncLower+truncUpper)))
}

gxy_sdlg <- function(x,y,params){
  grow_mean <- params$grow_sdlg_int
  grow <- dnbinom(x=y,mu=exp(grow_mean),size=exp(params$grow_sdlg_phi),log=F)
  truncLower<-dnbinom(x=0,mu=exp(grow_mean),size=exp(params$grow_sdlg_phi),log=F)
  truncUpper<-sum(dnbinom(x=params$max_size:10000,mu=exp(grow_mean),size=exp(params$grow_sdlg_phi),log=F))
  return(grow/(1-(truncLower+truncUpper)))
}

pxy<-function(x,y,params){
  sx(x,params) * gxy(x,y,params)
}

fx<-function(x, params){
  flw <- invlogit(params$flow_int + params$flow_slope*log(x))
  fert <- exp(params$fert_int + params$fert_slope*log(x))
  spike <- exp(params$spike_int + params$spike_slope*log(x))
  seeds_per_spike <- exp(params$seeds_per_spike)
  recruits_per_seed <- invlogit(params$recruits_per_seed)
  seedlings <- flw * fert * spike * seeds_per_spike * recruits_per_seed
  return(seedlings)
}

# Bigmatrix function ------------------------------------------------------
# This includes a reproductive delay till the first tiller
bigmatrix<-function(params){   
  matdim<-params$max_size
  y <- 1:matdim 
  #fertilitty transition
  Fmat <- matrix(0,matdim+1,matdim+1)
  Fmat[1,2:(matdim+1)]<-fx(x = y, params)
  
  #growth/survival transition
  Tmat <-matrix(0,matdim+1,matdim+1)
  Tmat[2:(matdim+1),2:(matdim+1)] <- t(outer(y,y,pxy,params))
  # surviving seedlings emerge into continuous population
  Tmat[2:(matdim+1),1] <- gxy_sdlg(x=1,y=y, params = params)*sx_sdlg(params = params)
  MPMmat<-Tmat + Fmat
  return(list(MPMmat = MPMmat, Tmat = Tmat, Fmat = Fmat))
}

# lambda2 function ------------------------------------------------------
# skips matrices with missing data. I'm not sure what the deal I've gotten an error when trying to run lambda for our posterior draws.
lambda2 <- function(mat) {
  # check mat for missing values: if TRUE return NA, else return lambda(mat)
  ifelse(any(is.na(mat)), NA, lambda(mat))     
}

