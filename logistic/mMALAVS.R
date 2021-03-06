set.seed(3)
library(mvtnorm)

# compute the log-posterior dist. -- assumes N(0,I) prior on beta vector
logpi = function(beta,y,X, dstar){
  eta = X%*%beta;
  return(sum(y*eta -log(1+exp(eta))) - .5 * sum( dstar * beta^2 ) )
}

# compute gradient based on a N(0,I) prior on beta
grad = function(beta, y, X, dstar) {
  eta = X%*%beta;
  return(-dstar * beta + t(X) %*% (y - (1+exp(-eta))^-1))
}

# compute inverse hessian based on a N(0,I) prior on beta
# in calderman parlance, this is G^-1
hessian = function(beta,y,X, dstar) {
  eta = X%*%beta;
  phat = exp(eta) / (1 + exp(eta)) 
  W = as.numeric(phat* (1-phat)) # Var(y)
  C = solve(diag(as.numeric(dstar)) + t(as.numeric(phat* (1-phat)) * X ) %*% X)
  return(C)
}

# mMALA beta proposal
beta_prop_mMALA = function(beta_old, epsilon, y, X, dstar) {
  p = dim(X)[2]
  C = hessian(beta_old, y, X, dstar)
  rootC = chol(C)
  return( beta_old + .5 * epsilon^2 * C %*% grad(beta_old, y, X, dstar) + epsilon * rootC %*% rnorm(p) )
}

# MMALA acceptance probability
acc_mMALA = function(beta_old, beta_prop, y, X, epsilon, dstar) {
  a = exp(logpi(beta_prop,y,X, dstar) - logpi(beta_old,y,X, dstar) )

  Ginv_old = hessian(beta_old, y, X, dstar)
  grad_old = grad(beta_old, y, X, dstar)
  
  Ginv_new = hessian(beta_prop, y, X, dstar)
  grad_new = grad(beta_prop, y, X, dstar)
  
  new_given_old = dmvnorm(as.numeric(beta_prop), mean= beta_old + .5 * epsilon^2 * Ginv_old %*% grad_old,
                          sigma = epsilon^2 * Ginv_old, log = T)
  
  old_given_new = dmvnorm(as.numeric(beta_old), mean= beta_prop + .5 * epsilon^2 * Ginv_new %*% grad_new,
                          sigma = epsilon^2 * Ginv_new, log = T)
  
  b = exp(old_given_new - new_given_old)
  return(a * b)
}


# functions for riemann manifold

dGdbetaj = function(beta, y, X, j) {
  n = nrow(X)
  eta = X%*%beta;
  phat = plogis(eta)
  W = as.numeric(phat* (1-phat)) # Var(y)
  Vj = (1-2*phat) * X[,j]
  return(t(X) %*% diag(W, nrow = n) %*% diag(as.numeric(Vj), nrow = n) %*% X)
}

# input: beta,y,X
# output: multidimensional array (p x p matrix for each beta)

dG = function(beta, y, X) {
  p = dim(X)[2]
  dGb = array(0, dim = c(p,p,p))
  for(m in 1:p) {
    dGb[,,m] = dGdbetaj(beta, y, X, m)
  }
  return(dGb)
}


# let's get first element
element1 = function(Ginv, beta, y, X) {
  
  dG_mda = dG(beta, y, X)
  
  p = dim(Ginv)[1]
  output1 = array(0,dim = c(p,1))
  output2 = array(0,dim = c(p,1))
  
  for (i in 1:p) {
    GinvdG = Ginv %*% dG_mda[,,i]  #dGdbetaj(beta,y,X,i) # p x p matrix
    output1 = output1 + (GinvdG %*% Ginv)[,i]
    output2 = output2 + Ginv[,i] * sum(diag(GinvdG))
  }
  return(list(output1, output2))
}

# propose a new beta based on old gradient/ hessian/ hot

beta_prop_M2MALA = function(beta_old, epsilon, y, X, dstar) {
  p = dim(X)[2]
  gradient = beta_old$grad
  C = beta_old$hess
  #dG = dG(beta_old,y,X)
  rootC = chol(C)
  elements = beta_old$elements #element1(C, dG, beta_old, y, X)
  mu_prop = beta_old$beta + .5 * epsilon^2 * C %*% gradient - epsilon^2 * elements[[1]] + 
    .5 * epsilon^2 * elements[[2]]
  beta_prop = mu_prop + epsilon * rootC %*% rnorm(p)
  
  # now get moments associated with the beta proposal 
  
  gnew = grad(beta_prop, y, X, dstar)
  hnew = hessian(beta_prop, y, X, dstar)
  elements = element1(hnew, beta_prop, y, X)
  
  # return proposed beta along with its moments
  return(list(beta = beta_prop, gradient = gnew, hess = hnew, elements = elements))
}

# MMALA acceptance probability
# input: beta_prop: list containing the beta proposal as well as the other terms to get the mean

acc_M2MALA = function(beta_old, beta_prop, y, X, epsilon, dstar) {
  a = exp(logpi(beta_prop$beta, y, X, dstar) - logpi(beta_old$beta,y,X, dstar) )
  
  mu_old = beta_old$beta - epsilon^2 * beta_old$elements[[1]] + .5 * epsilon^2 * beta_old$elements[[2]] + 
             .5 * epsilon^2 * beta_old$hess %*% beta_old$gradient
  
  mu_new = beta_prop$beta - epsilon^2 * beta_prop$elements[[1]] + .5 * epsilon^2 * beta_prop$elements[[2]] + 
           .5 * epsilon^2 * beta_prop$hess %*% beta_prop$gradient
  
  new_given_old = dmvnorm(as.numeric(beta_prop$beta), mean= mu_old,
                          sigma = epsilon^2 * beta_old$hess, log = T)
  
  old_given_new = dmvnorm(as.numeric(beta_old$beta), mean= mu_new,
                          sigma = epsilon^2 * beta_prop$hess, log = T)
  
  b = exp(old_given_new - new_given_old)
  return(a * b)
}


M2MALA_logRegr = function(Niter, y, X, epsilon, beta.init){
  
  p=length(X[1,]); n=length(y);
  Res = matrix(NA,ncol=p,nrow = Niter);
  beta = beta.init; #initialize with zeros
  dstar = rep(1, p)
  
  # gradient
  gradient = grad(beta, y, X, dstar)
  # hessian
  hess = hessian(beta, y, X, dstar)
  # terms
  elements = element1(hess, beta, y, X)
  
  # beta is a list that contains beta as well as the 'terms'
  beta = list(beta=beta, gradient = gradient, hess = hess, elements = elements)
  
  lpi = logpi(beta$beta, y, X, dstar);
  
  for (jj in 1:Niter) {
    print(jj)
    
    # beta prop should return gradient, hessian, elements
    beta_prop =  beta_prop_M2MALA(beta, epsilon, y, X, dstar) 
    alpha = acc_M2MALA(beta_old = beta, beta_prop, y, X, epsilon, dstar)
    Acc = min(1,alpha)
    if (runif(1)<=Acc){
      beta = beta_prop;
    }
    epsilon = epsilon + (1/jj^0.7)*(Acc - 0.5);
    Res[jj,] = beta$beta
  }
  return(Res)
}

# dstar contains diagonal elements of D

mMALA_logRegr = function(Niter, y, X, epsilon, beta.init, theta.init, pstar, v0, v1) { #, theta.init, gamma.init){
  p=length(X[1,]); n=length(y);
  Res = matrix(NA,ncol=p * 2, nrow = Niter);
  dimnames(Res)[[2]] = c( paste("beta",1:p), paste("gamma",1:p) )
  
  # initialize parameters
  beta = beta.init 
  gam = ifelse(pstar > .2, 1, 0)
  dstar = ifelse(gam==1, 1/v1, 1/v0)
  theta = theta.init
  
  lpi = logpi(beta,y,X, dstar);
  
  # iterate
  for ( jj in 1:Niter){
    print(jj)
    
    # sample beta
    beta_prop =  beta_prop_mMALA(beta, epsilon, y, X, dstar) 

    alpha = acc_mMALA(beta_old = beta, beta_prop, y, X, epsilon, dstar)
    Acc = min(1,alpha)
    if (runif(1)<=Acc){
      beta= beta_prop;
    }
    
    # sample theta
    theta = rbeta(1, 1 + sum(gam), 1 + p - sum(gam))
    
    # sample gamma
    p1 = dnorm(beta, sd = sqrt(v1)) * theta
    p0 = dnorm(beta, sd = sqrt(v0)) * (1-theta)
    gam = rbinom(p, 1, prob = p1/(p1 + p0))
    dstar = ifelse(gam==1, 1/v1, 1/v0)
    
    epsilon = epsilon + (1/jj^0.7)*(Acc - 0.5);
    
    Res[jj,] = c(beta, gam)
    
  }
  return(Res)
}
