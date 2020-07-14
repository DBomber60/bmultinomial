# replicating the rockava george EM for sanity - normal outcome

EM.iter = function(beta.current, theta.current, nu_0, nu_1, sigma.current = 1) {
  
  # E-Step: get pstar and dstar vectors
  
  # when gamma_i=1, beta_i ~ N(0, v_1); beta_i is a scale mixture of normals (spike + slab)
  p1 = dnorm(beta.current, mean = 0, sd = sqrt(nu_1*3)) * theta.current
  p0 = dnorm(beta.current, mean = 0, sd = sqrt(nu_0*3)) * (1-theta.current)
  
  pstar = p1/(p0 + p1)
  dstar = pstar/nu_1 + (1-pstar)/nu_0
  
  # M-Step
  
  beta.new = solve(t(X) %*% X + diag(as.numeric(dstar)) ) %*% t(X) %*% Y
  theta.new = sum(pstar)/(p+1) # assumes a=b=1
  
  #print(theta.new)
  
  return(list(beta=beta.new, theta = theta.new, pstar = pstar))
}

EMVS2 <- function(beta.init, theta.init, nu_0, nu_1) {
  beta.current <- beta.init
  theta.current = theta.init
  
  delta.ll <- 1
  while(delta.ll > 1e-8)  {
    beta.old = beta.current
    it = EM.iter(beta.current, theta.current, nu_0, nu_1)
    beta.current = it$beta
    theta.current = it$theta
    delta.ll = sqrt(crossprod(beta.current - beta.old))
    dstar = it$pstar/nu_1 + (1-it$pstar)/nu_0
  }
  
  return(list(beta.current, theta.current, dstar))
}
