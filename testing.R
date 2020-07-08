rm(list=ls())
source('irls.R')
source('emvs.R')
library(tidyverse)
set.seed(1)

# generate data
n = 500
p = 10
X = cbind(1, matrix( rnorm(n * p), nrow = n) )
beta.true = c(1, 1, 1, 1, 1.5, .2, 0, 0, 0, 0, 0)
eta = X %*% beta.true
Y = rbinom(n, 1, prob = plogis(eta))

#summary(glm(Y~X-1, family = binomial))
#binomial.irls(Y~X-1)$coef

# initialize parameters

nu_0 = 0.5
nu_1 = 10
theta_0 = 0.5
beta_0 = rep(0, p+1)

j = EMVS(beta_0, theta_0, nu_0 = .5, nu_1 = 100)




