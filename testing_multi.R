library(nnet)
source('bayesmulti.R')
set.seed(1)

# 1. does it work for proportions? CHECK
# 2. add weights to maximizer CHECK
# 3. test weights CHECK
# 4. tweak it to make it bayesian! CHECK
# 5. substitute into current EM setup

# maximize beta given phat.current/ beta.1.nr

nr.beta.mult = function(beta.init, dstar, X, Y, K, w = rep(1, nrow(Y))) {
  
  beta.current= beta.init
  dnr.ll = 1
  
  while (dnr.ll > 1e-8) { #  
      beta.old = beta.current
      beta.current.m = matrix(beta.current, nrow = p, ncol = K-1)
      e_eta = exp( X %*% beta.current.m )
      phat_m = e_eta / (1 + rowSums(e_eta))
      # Y = unw
      # w = pr
      U = grad_multi(B = beta.current, phat_m, X, Y, K, w, dstar)
      H.curr = H(phat_m, X, K, w, dstar)
      beta.current = beta.current + H.curr %*% U
      dnr.ll = sqrt(crossprod(beta.current - beta.old))
    }
  
  return(beta.current)
  
}

# 1. does maximizer work?

# simulate multinomial data
n = 500
K = 3 # 3 classes
beta.true = matrix(c(c(0,0),
                     c(0,.5),
                     c(0,3)), nrow = 2)

#beta.test = matrix(c(c(0,0), c(0,.5)), nrow = 2)
p = 2 # total covariates (including intercept)
X = cbind(rep(1,n), matrix(rnorm(n * (p-1)), ncol = (p-1)))

e_eta = exp(X %*% beta.true)
probs = e_eta / rowSums(e_eta)
Y = t(apply(probs, 1, function(x) rmultinom(1, 1, prob = x)))
yind = apply(Y, 1, function(x) which(x==1)) # factor response

summary(nnet::multinom(Y ~ X -1))
beta.init = rep(1, (k-1)*p)

#a = (nr.beta(beta.init, dstar = 0, X, Y, K))

library(gtools)

gammas = rdirichlet(500, alpha = c(1,2,3))

# ok, now let's try with weights
pr = rep(.5, 500) #rbeta(500, shape1 = 1, shape2 = 2)

unw = gammas * pr

summary(nnet::multinom(unw ~ X - 1, weights = pr))
b = (nr.beta.mult(beta.init, dstar = rep(0, (K - 1) * p), X, Y = unw, K, w = pr))






