# implementation of the EMVS procedure of Rockava and George, 2014

EM.iter = function(beta.current, theta.current, nu_0, nu_1) {

  # E-Step: get pstar and dstar vectors
  
  # when gamma_i=1, beta_i ~ N(0, v_1); beta_i is a scale mixture of normals (spike + slab)
  p1 = dnorm(beta.current, mean = 0, sd = nu_1) * theta.current
  p0 = dnorm(beta.current, mean = 0, sd = nu_0) * (1-theta.current)
  pstar = p1/ (p0 + p1)
  dstar = pstar/nu_1 + (1-pstar)/nu_0
  
  
  # M-Step
  
  # one step of Newtwon Rhapson to update beta vector
  eta.current = X %*% beta.current
  phat.current = plogis(eta.current)
  U = t(X) %*% (Y - phat.current) - dstar * beta.current
  H = solve(diag(dstar) + t(as.numeric(phat.current * (1-phat.current)) * X ) %*% X)
  beta.new = beta.current + H %*% U
  
  # now, update theta
  theta.new = sum(pstar)/(p+1) # assumes a=b=1
  
  return(list(beta=beta.new, theta = theta.new))
  
}


EMVS <- function(beta.init, theta.init, nu_0, nu_1) {
  
  beta.current <- beta.init
  theta.current = theta.init
  
  delta.ll <- 1
  
  while(delta.ll > 1e-6) {
    beta.old = beta.current
    it = EM.iter(beta.current, theta.current, nu_0, nu_1)
    beta.current = it$beta
    theta.current = it$theta
    delta.ll = sqrt(crossprod(beta.current - beta.old))
  }
  
  return(list(beta.current, theta.current))
}
