rm(list=ls())
source('irls.R')
source('emvs_RG.R')
library(tidyverse)
library(mvtnorm)
library(latex2exp)
set.seed(1)

# generate data
n = 100
p = 100

# correlated predictors
# make a covariance matrix

sig = matrix(nrow = p, ncol = p)

for(i in 1:p) {
  for(j in 1:p) {
    sig[i,j] = abs(i-j)
  }
}

sig.5 = (exp(sig * log(.6)))
X = rmvnorm(n, sigma = sig.5)
beta.true = c(3,2,1, rep(0, p-3))
Y = rnorm(n, mean = X %*% beta.true, sd = sqrt(3))

beta_0 = rep(1,p)
theta_0 = 0.5
j = EMVS2(beta.init = beta_0, theta.init = theta_0, nu_0 = .1, nu_1 = 1000)
j[[3]]

plot(j[[1]])

# spike and slab regularization plot
theta_0 = 0.5
nu_0_seq = seq(from = 0.005, to = .5, by = .01)
res = matrix(nrow = p, ncol = length(nu_0_seq))
colnames(res) = paste("nu0",nu_0_seq)
a = data.frame(res)
a$var = factor(1:p) #paste("beta",1:(p+1))

# now, fill in the data

for (m in 1:length(nu_0_seq)) {
  em.iter = EMVS2(beta_0, theta_0, nu_0_seq[m], nu_1=1000)
  a[,m] = em.iter[[1]]
}

long = pivot_longer(a, cols = starts_with("nu0"), values_to = "beta")
long$nu = rep(nu_0_seq, p)

ggplot(data = long, aes(x=nu, y=beta[,1])) +
  geom_line(aes(colour=var)) +
  scale_color_manual(values=c( rep('#E69F00', 3), rep('#999999', p-3)) ) +
  theme_bw() + ylab(c(-5,3)) + 
  theme(legend.position = "none") + ylab(TeX("$\\hat{\\beta}$")) +
  xlab(TeX("$\\nu_0$"))


thresh = function(v_0) {
  c = sqrt(v_1/v_0)
}

