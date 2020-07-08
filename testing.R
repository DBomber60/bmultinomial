rm(list=ls())
source('irls.R')
source('emvs.R')
library(tidyverse)
set.seed(1)

# generate data
n = 500
p = 10
X = cbind(1, matrix( rnorm(n * p), nrow = n) )
beta.true = c(.5, .5, 1, 1, 1.5, .2, 0, 0, 0, 0, 0)
eta = X %*% beta.true
Y = rbinom(n, 1, prob = plogis(eta))

#summary(glm(Y~X-1, family = binomial))
#binomial.irls(Y~X-1)$coef

# initialize parameters

nu_0 = 0.1
nu_1 = 10
theta_0 = 0.5
beta_0 = rep(0, p+1)

EMVS(beta_0, theta_0, nu_0 = 5, nu_1)

# spike and slab regularization plot
nu_0_seq = seq(from = 0.01, to = 2.0, by = .25)
res = matrix(nrow = p+1, ncol = length(nu_0_seq))
colnames(res) = paste("nu0",nu_0_seq)
a = data.frame(res)
a$var = paste("beta",1:(p+1))

# now, fill in the data

for (m in 1:length(nu_0_seq)) {
  em.iter = EMVS(beta_0, theta_0, nu_0_seq[m], nu_1)
  a[,m] = em.iter[[1]]
}


long = pivot_longer(a, cols = starts_with("nu0"), values_to = "beta")
long$nu = rep(nu_0_seq, p+1)

ggplot(data = long, aes(x=nu, y=beta[,1])) + geom_line(aes(colour=var)) + 
  theme_bw()

