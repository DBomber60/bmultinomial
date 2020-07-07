source('irls.R')
set.seed(1)

# first problem to solve: marginalizing out the intercept in logistic regression
set.seed(1)
n = 500
p = 5
X = cbind(1, matrix( rnorm(n * p), nrow = n) )
beta.true = c(.5, .5, 1, 1, 1.5, .2)
eta = X %*% beta.true
Y = rbinom(n, 1, prob = plogis(eta))


summary(glm(Y~X-1, family = binomial))
binomial.irls(Y~X-1)$coef

# centered Y

# initialize parameters

nu_0 = 0.1
nu_1 = 10
theta_0 = 0.5
beta_0 = rep(0, 6)

beta.current = beta_0
theta.current = theta_0


EM.iter = function(beta.current, theta.current) {

  # E-Step: get pstar and dstar vectors
  
  # when gamma_i=1, beta_i ~ N(0, v_1)
  p1 = dnorm(beta.current, mean = 0, sd = nu_1) * theta.current
  p0 = dnorm(beta.current, mean = 0, sd = nu_0) * (1-theta.current)
  pstar = p1/ (p0 + p1)
  print(pstar)
  dstar = pstar/nu_1 + (1-pstar)/nu_0
  
  # M-Step
  # one step of Newtwon Rhapson
  eta.current = X %*% beta.current
  phat.current = plogis(eta.current)
  U = t(X) %*% (Y - phat.current) - dstar * beta.current
  H = solve(diag(dstar) + t(as.numeric(phat.current * (1-phat.current)) * X ) %*% X)
  beta.new = beta.current + H %*% U
  theta.new = sum(pstar)/(p+1) # assumes a=b=1
  
  return(list(beta=beta.new, theta = theta.new))
  
}


EMVS <- function(beta.init, theta.init, nu_0, nu_1) {
  
  beta.current <- beta.init
  theta.current = theta.init
  
  delta.ll <- 1
  
  while(delta.ll > 1e-6) {
    beta.old = beta.current
    it = EM.iter(beta.current, theta.current)
    beta.current = it$beta
    theta.current = it$theta
    delta.ll = sqrt(crossprod(beta.current - beta.old))
  }
  
  return(list(beta.current, theta.current))
}

EMVS(beta_0, theta_0)

# 0.04220186