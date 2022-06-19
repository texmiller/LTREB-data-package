library(rstan)

## one random normal deviate and one random bernoulli deviate
y1 <- rnorm(100,50,10)
y2 <- rnorm(150,10,10)

dat <- list(
  y1=y1,
  N1=length(y1)
)

fit_combo <- stan(file = "combo_model.stan",data = dat,
                  chains = 2, iter = 1000, warmup = 100)

