
max_size <- data.frame(max_size=c(18,16,9,13,40,28,21))

toy_params <- make_params(species=1,
                          endo_mean=0,
                          endo_var=0,
                          original = 1, # should be =1 to represent recruit
                          draw=1,
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
                          recruit_par=recruit_par)

gxy <- function(x,y,params){
  grow_mean <- params$grow_int + params$grow_slope*log(x)
  grow<-dpois(x=y,lambda=exp(grow_mean))
  truncZero<-dpois(x=0,lambda=exp(grow_mean))
  truncUpper<-sum(dpois(x=(params$max_size+1):50,lambda=exp(grow_mean)))
  return(grow/(1-(truncZero+truncUpper)))
}

g_mat <- t(outer(1:(toy_params$max_size),1:(toy_params$max_size),gxy,toy_params))
image(g_mat)
colSums(g_mat)

size_t=18
plot(g_mat[,size_t],type="l",main=paste(sum(g_mat[,size_t])))
abline(v=size_t)

dpoisinvgauss(x=0,mean=5, shape=2)
              