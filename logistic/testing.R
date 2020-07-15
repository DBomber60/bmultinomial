rm(list=ls())
source('irls.R')
source('logistic/emvs.R')
source('plotting.R')
source('logistic/mMALAVS.R')
library(tidyverse)
library(mvtnorm)
set.seed(1)

# to do on logistic model (before moving onto multinomial)
# 1. adjust EMVS so that intercept is always in model
# 2. adjust plotting for EMVS to choose the right hyperparameters - CHECK
# 3. get full gibbs sampler to work on mMALA case
# 4. make diagnostic plots from MCMC

# does MMALA scale better with more predictors ??

# correlated predictors
# make a covariance matrix
p = 10
n = 500

sig = matrix(nrow = p, ncol = p)

for(i in 1:p) {
  for(j in 1:p) {
    sig[i,j] = abs(i-j)
  }
}

sig.5 = (exp(sig * log(.6)))
X = cbind(rmvnorm(n, sigma = sig.5))
#X = rmvnorm(n, mean = rep(0,p))

beta.true = c(1, 1, 1, 1.5, 3, rep(0, p-5))
eta = X %*% beta.true
Y = rbinom(n, 1, prob = plogis(eta))

summary(glm(Y~X-1, family = binomial))
#binomial.irls(Y~X-1)$coef
#nr.beta( rep(0,p), phat.init = rep(.5, n), dstar = rep(0, p) )

# initialize parameters
beta_0 = rep(0, p)

# send this to plotter for a visually informed choice of hyperparameters
#plot_emvs(p, 5, beta_0, theta.init = .5, nu0_low = 0.01, nu0_high = 1, increment = .01, nu_1=1000)


init.vals = EMVS(beta_0, theta.init = .5, nu_0 = .5, nu_1 = 100)

#Call the MCMC Sampler
epsilon = .05;
Niter = 500;
burnin = floor(0.2 * Niter)
Res.mM = M2MALA_logRegr(Niter,Y,X,epsilon, beta.init = init.vals[[1]])

# in the MCMC, beta is a list that contains the beta vector, gradient, hessian, and terms

plotter = function(output, whichbeta, burnin, Niter, type) {
  vec = output[burnin:Niter,whichbeta]
  plot(vec,type='l',col = rgb(red = 1, green = 0, blue = 0, alpha = 0.5), main = type)
  acf(vec,lag.max=100,col='blue', main = "")
  hist(vec,nclass=50,prob=T,col='blue', main = "")
}

par(mfrow = c(1,3))
plotter(Res.mM, whichbeta = 1, burnin, Niter, type = "MMALA")

