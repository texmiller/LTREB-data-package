//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N1;
  vector[N1] y1;
  int<lower=0> n_blocks1;
  int<lower=0> block1[N1];
  int<lower=0> N2;
  int y2[N2];
  int<lower=0> n_blocks2;
  int<lower=0> block2[N2];
  
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real mu1;
  real eps1[n_blocks1];
  real<lower=0> sigma;
  real<lower=0,upper=1> block1_sigma;  
  real mu2;
  real eps2[n_blocks2];
  real<lower=0,upper=1> block2_sigma;  
  }
  
transformed parameters {
  real lp1[N1];
  real lp2[N2];
  
  for(i in 1:N1){
    lp1[i] = mu1 + eps1[block1[i]];
  }
  for(i in 1:N2){
    lp2[i] = mu2 + eps2[block2[i]];
  }
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  eps1~normal(0,block1_sigma);
  eps2~normal(0,block2_sigma);
  
  y1 ~ normal(lp1, sigma);
  y2 ~ bernoulli_logit(lp2);
}

