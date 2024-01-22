library(RColorBrewer)

rm(list=ls(all=TRUE)); graphics.off();

d = read.csv ("LDW_LTREB_20072022.csv")

## if there is size data but flw_count_t (# flowering tillers at the
## start of the interval) is NA then set infs ("inflorescences") to
## zero
d$infs.t = ifelse(!is.na(d$size_t) & is.na(d$flw_count_t), 0, d$flw_count_t)
d$infs.t1 = ifelse(!is.na(d$size_t1) & is.na(d$flw_count_t1), 0,
                   d$flw_count_t1)

## Calculate age
d$age = d$year_t - d$birth

## Only keep rows for which the individual was recruited into the
## population, not greenhouse-reared and transplanted.
d = d[d$origin_01 == 1,]

## Only keep rows with non-negative ages
d = d[d$age >= 0,]

##If we want to quantify LRO in terms of seed production rather than
##inflorescence production then we would want to use the spike
##info. Most of the variation in reproductive output comes in the form
##of variable numbers of inflorescences and less so from variable
##numbers of spikelets within an inflorescence. So in this first pass
##I went with the coarser measure of reproduction.
if (FALSE) { 
  ## Calculate average number of flower spikes
  aveSpike.t = apply(cbind(d$spike_a_t, d$spike_b_t, d$spike_c_t,
                           d$spike_d_t), 1, mean, na.rm=TRUE)
  aveSpike.t1 = apply(cbind(d$spike_a_t1, d$spike_b_t1, d$spike_c_t1,
                            d$spike_d_t1), 1, mean, na.rm=TRUE)


  ## where spike data is missing, fill with species mean
  meanSpikes.t = tapply (aveSpike.t, d$species, mean, na.rm=TRUE)
  meanSpikes.t1 = tapply (aveSpike.t1, d$species, mean, na.rm=TRUE)
  meanSpikelets = data.frame(meanSpikes.t=meanSpikes.t,
                             meanSpikes.t1=meanSpikes.t1)
}

## remove anything still alive in 2022 (our last year of data
## collection) because it's not done living and reproducing
alive22 = ifelse (d$year_t1==2022 & d$surv_t1==1, 1, 0)
d = d[alive22 != 1,]

## calculate lifespan for each individual in each species
##lifespan = tapply (d$age, list(d$species, d$id), max)
lifespan = tapply (d$age, d$id, max)

## calculate total inflorescences
##totInf = tapply (d$flw_count_t, list(d$species, d$id), sum,
##na.rm=TRUE)
R = tapply (d$infs.t, d$id, sum, na.rm=TRUE)

speciesForID = tapply (d$species, d$id, unique)

dd = data.frame (lifespan=lifespan, R=R, species=speciesForID)

spp_list<-unique(dd$species)
spp_cols<-brewer.pal(n=7,name="Dark2")
spp_names<-c("Agrostis perennans","Elymus villosus","Elymus virginicus",
            "Festuca subverticillata","Lolium aryndinaceum",
            "Poa alsodes","Poa sylvestris")

par(mfrow=c(2,4))
for(i in 1:length(spp_list)) {
  plot (jitter (dd$lifespan), jitter(dd$R), main=spp_names[i],
        xlab="Lifespan", ylab="LRO")
}

##maxLifespan = tapply (dd$lifespan, dd$species, max)
##maxR = tapply (dd$R, dd$species, max)

speciesCodes = unique (dd$species)
dev.new()
par(mfrow=c(2,4))
for (ii in 1:7) {
  subset = dd[dd$species==speciesCodes[ii],]
  maxLifespan = max(subset$lifespan)
  maxR = max(subset$R)
  jointProbLifespanR = probLifespanCondR =
    matrix (0, nrow=maxR, ncol=maxLifespan)
  probR = rep(0, maxR)
  for (r in 1:maxR) {
    probR[r] = sum(subset$R==r) / length(subset$R)
    for (l in 1:maxLifespan) {
      jointProbLifespanR[r,l] =
        sum(subset$lifespan==l & subset$R == r) /
        length(subset$lifespan)
    }
    probLifespanCondR[r,] = jointProbLifespanR[r,] /
      probR[r]
  }

  image (0:maxR, 0:maxLifespan, probLifespanCondR, xlab="LRO",
         ylab="Lifespan", main=spp_names[ii])
}

library(qgam)
plot(jitter(dd$R),jitter(dd$lifespan))
gam05<-qgam(lifespan ~ s(R), qu=0.05, data=dd)
gam95<-qgam(lifespan ~ s(R), qu=0.05, data=dd)

library(MASS)
attach(geyser)
plot(jitter(dd$R),jitter(dd$lifespan))
f1 <- kde2d(dd$R, waiting, n = 50, lims = c(0.5, 6, 40, 100))
image(f1, zlim = c(0, 0.05))
f2 <- kde2d(duration, waiting, n = 50, lims = c(0.5, 6, 40, 100),
            h = c(width.SJ(duration), width.SJ(waiting)) )
image(f2, zlim = c(0, 0.05))
persp(f2, phi = 30, theta = 20, d = 5)
