
lambda<-rnorm(15,1,0.15)
toptwo <- lambda%in%sort(lambda)[c(1,15)]
topfour <- lambda%in%sort(lambda)[c(1:2,14:15)]
topsix <- lambda%in%sort(lambda)[c(1:3,13:15)]

par(mfrow=c(4,1),mar=c(2,4,0,0))
hist(lambda,xlim=c(0.7,1.3),breaks=10,main=" ",ylim=c(0,4),col="darkgray")
hist(lambda[topsix],xlim=c(0.7,1.3),breaks=10,main=" ",ylim=c(0,4),col="forestgreen")
hist(lambda[topfour],xlim=c(0.7,1.3),breaks=10,main=" ",ylim=c(0,4),col="cornflowerblue")
hist(lambda[toptwo],xlim=c(0.7,1.3),breaks=10,main=" ",ylim=c(0,4),col="tomato")

sd(lambda);sd(lambda[topsix]);sd(lambda[topfour]);sd(lambda[toptwo])
mean(lambda);mean(lambda[topsix]);mean(lambda[topfour]);mean(lambda[toptwo])
