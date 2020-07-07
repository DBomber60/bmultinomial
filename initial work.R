source('irls.R')

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

# E-Step: get pstar and dstar vectors

# when gamma_i=1, beta_i ~ N(0, v_1)
p1 = dnorm(beta_0, mean = 0, sd = nu_1) * theta_0
p0 = dnorm(beta_0, mean = 0, sd = nu_0) * (1-theta_0)
pstar = p1/ (p0 + p1)
dstar = pstar/nu_1 + (1-pstar)/nu_0

# M-Step





