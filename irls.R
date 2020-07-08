# [ Iteratively reweighted least squares for canonical links ]

# `f` is formula; `b`, `b1`, `b1inv` and `b2` are the cumulant function, its
# first derivative, inverse of b', and second derivative; `init`
# is a function that provides an initial 'mu' as a function of 'y'.
# For simplicity, assumes unit dispersion and canonical link
irls <- function (f, b, b1, b1inv, b2, init, tol = 1e-6) {
  # initialize
  y <- model.response(model.frame(f))
  X <- model.matrix(f) # design matrix
  mu <- init(y)
  eta <- b1inv(mu) # g(mu)
  
  lhood <- sum(y * eta - b(eta))

  # iterate
  repeat {
    W <- sapply(eta, b2) # b''(theta) = V(mu)
    z <- eta + (y - mu) / W # working response
    wlm <- lm(z ~ X - 1, weights = W) # weighted least squares
    beta <- coef(wlm)
    eta <- fitted(wlm) # X %*% beta
    lhood.new <- sum(y * eta - b(eta))
    if (abs((lhood.new - lhood) / lhood) < tol) break # converged?
    lhood <- lhood.new
    mu <- sapply(eta, b1) # b'(theta)
  }

  # report
  list(coef = beta, var = summary(wlm)$cov.unscaled)
}

# e.g. Poisson with canonical link:
poisson.irls <- function (f, tol = 1e-6)
  irls(f, exp, exp, log, exp, function (y) y + 0.5, tol)

# e.g. binomial with canonical link:
b = function(theta) log(1+exp(theta))
b2 = function(theta) plogis(theta)*(1-plogis(theta))
logit = function(x) log(x/(1-x))

binomial.irls <- function (f, tol = 1e-6)
  irls(f, b=b, b1 = plogis, b1inv = logit, b2 = b2, function (y) ifelse(y>.5,.7,.2), tol)


