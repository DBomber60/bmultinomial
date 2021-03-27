library(nnet)
source('bayesmulti.R')
set.seed(1)

# 1. make separate EM
# 2. substitute into current larger EM setup

# maximize beta given phat.current/ beta.1.nr

nr.beta.mult = function(beta.init, dstar, X, Y, K, w = rep(1, nrow(Y))) {
  
  beta.current= beta.init
  p = dim(X)[2]
  dnr.ll = 1
  
  while (dnr.ll > 1e-6) { #  
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

EM.iter = function(beta.current, theta.current, nu_0, nu_1, X, Y, K, w) {
  
  # E-Step: get pstar and dstar vectors
  
  # when gamma_i=1, beta_i ~ N(0, v_1); beta_i is a scale mixture of normals (spike + slab)
  p1 = dnorm(beta.current, mean = 0, sd = sqrt(nu_1)) * theta.current
  p0 = dnorm(beta.current, mean = 0, sd = sqrt(nu_0)) * (1-theta.current)
  
  pstar = p1/(p0 + p1)
  dstar = pstar/nu_1 + (1-pstar)/nu_0

  
  # M-Step
  
  # maximize beta using Newtwon Rhapson
  
  beta.new = nr.beta.mult(beta.current, dstar, X, Y, K, w)
  
  # now, update theta
  theta.new = sum(pstar)/((K-1)*p) # assumes a=b=1
  return(list(beta=beta.new, theta = theta.new, pstar = pstar))
}


EMVS <- function(beta.init, theta.init, nu_0, nu_1,  X, Y, K, w) {
  
  beta.current <- beta.init
  theta.current = theta.init
  
  delta.ll <- 1
  while(delta.ll > 1e-6)  { # while(delta.ll > 1e-6) for (i in 1:1)
    beta.old = beta.current
    it = EM.iter(beta.current, theta.current, nu_0, nu_1, X, Y, K, w)
    beta.current = it$beta
    theta.current = it$theta
    delta.ll = sqrt(crossprod(beta.current - beta.old))
    print(delta.ll)
    pstar = it$pstar
  }
  return(list(beta.current, theta.current, pstar))
}

# simulate multinomial data
n = 500
K = 3 # 3 classes
p = 10 # total covariates (including intercept)

beta.true = matrix(  c(1, 2, 3, rep(0,7),
                       2, 3, 4, rep(0,7)) , nrow = p)

X = cbind(rep(1,n), matrix(rnorm(n * (p-1)), ncol = (p-1)))

e_eta = exp(X %*% cbind(0, beta.true) )
probs = e_eta / rowSums(e_eta)
Y = t(apply(probs, 1, function(x) rmultinom(1, 1, prob = x)))
yind = apply(Y, 1, function(x) which(x==1)) # factor response

#summary(nnet::multinom(Y ~ X -1))
beta.init = rep(0, (K-1)*p)
#j = nr.beta.mult(beta.init, dstar = rep(0, 20), X, Y, K)
#print(matrix(j, nrow = p))

#a = EMVS(beta.init, theta.init = .5, nu_0 = .1, nu_1 = 10, X, Y, K, w = rep(1, nrow(Y)))

