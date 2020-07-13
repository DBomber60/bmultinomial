set.seed(3)
source('bayesLRsampler.R')
library(mvtnorm)

# compute the log-posterior dist. -- assumes N(0,I) prior on beta vector
logpi = function(beta,y,X){
  eta = X%*%beta;
  return(sum(y*eta -log(1+exp(eta))) - (0.5)*crossprod(beta))
}

# compute gradient based on a N(0,I) prior on beta
grad = function(beta,y,X) {
  eta = X%*%beta;
  return(-beta + t(X) %*% (y - (1+exp(-eta))^-1))
}

# compute inverse hessian based on a N(0,I) prior on beta
# in calderman parlance, this is G^-1
hessian = function(beta,y,X) {
  eta = X%*%beta;
  phat = exp(eta) / (1 + exp(eta)) 
  W = as.numeric(phat* (1-phat)) # Var(y)
  C = solve(diag(ncol(X)) + t(as.numeric(phat* (1-phat)) * X ) %*% X)
  return(C)
}

# functions for riemann manifold
dGdbetaj = function(beta, y, X, j) {
  n = nrow(X)
  eta = X%*%beta;
  phat = plogis(eta)
  W = as.numeric(phat* (1-phat)) # Var(y)
  Vj = (1-2*phat) * X[,j]
  return(t(X) %*% diag(W) %*% diag(as.numeric(Vj), nrow = n) %*% X)
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


element1 = function(Ginv, dG, beta, y, X) {
  p = dim(Ginv)[1]
  output1 = array(0,dim = c(p,1))
  output2 = array(0,dim = c(p,1))
  
  for (i in 1:p) {
    GinvdG = Ginv %*% dGdbetaj(beta,y,X,i) # p x p matrix
    output1 = output1 + (GinvdG %*% Ginv)[,i]
    output2 = output2 + Ginv[,i] * sum(diag(GinvdG))
  }
  return(list(output1, output2))
}

#element1(Ginv, dG, beta,y,X)


beta_prop_M2MALA = function(beta_old, epsilon, y, X) {
  p = dim(X)[2]
  C = hessian(beta_old, y, X)
  dG = dG(beta_old,y,X)
  rootC = chol(C)
  terms = element1(C, dG, beta_old, y, X)
  return(beta_old + .5 * epsilon^2 * C %*% grad(beta_old, y, X) + epsilon * rootC %*% rnorm(p) -
           .5 * epsilon^2 * terms[[1]] )
}



# different beta proposals based on sampling method

# MALA beta proposal
beta_prop_MALA = function(beta_old, epsilon, y, X) {
  return( beta_old + as.numeric(.5 * epsilon^2 * grad(beta_old, y, X) + epsilon * rnorm(p)) )
}

# mMALA beta proposal
beta_prop_MMALA = function(beta_old, epsilon, y, X) {
  p = dim(X)[2]
  C = hessian(beta_old, y, X)
  rootC = chol(C)
  return( beta_old + .5 * epsilon^2 * C %*% grad(beta_old, y, X) + epsilon * rootC %*% rnorm(p) )
}

# MMALA beta proposal

# MALA acceptance probability
acc_MALA = function(beta_old, beta_prop, y, X, epsilon) {
  a = exp(logpi(beta_prop,y,X) - logpi(beta_old,y,X) )

  old_given_new = crossprod(beta_old - beta_prop - .5 * epsilon^2 * grad(beta_prop, y, X))
  new_given_old = crossprod(beta_prop - beta_old - .5 * epsilon^2 * grad(beta_old, y, X))
  b = exp( -1/(2*epsilon^2) * (old_given_new - new_given_old) )

  return(a * b)
}

# MMALA acceptance probability
acc_MMALA = function(beta_old, beta_prop, y, X, epsilon) {
  a = exp(logpi(beta_prop,y,X) - logpi(beta_old,y,X) )
  
  Ginv_old = hessian(beta_old, y, X)
  grad_old = grad(beta_old, y, X)
  
  Ginv_new = hessian(beta_prop, y, X)
  grad_new = grad(beta_prop, y, X)
  
  new_given_old = dmvnorm(as.numeric(beta_prop), mean= beta_old + .5 * epsilon^2 * Ginv_old %*% grad_old,
                          sigma = epsilon^2 * Ginv_old, log = T)
  
  old_given_new = dmvnorm(as.numeric(beta_old), mean= beta_prop + .5 * epsilon^2 * Ginv_new %*% grad_new,
                          sigma = epsilon^2 * Ginv_new, log = T)
  
  b = exp(old_given_new - new_given_old)
  return(a * b)
}

# MMALA acceptance probability
acc_M2ALA = function(beta_old, beta_prop, y, X, epsilon) {
  a = exp(logpi(beta_prop,y,X) - logpi(beta_old,y,X) )
  
  Ginv_old = hessian(beta_old, y, X)
  grad_old = grad(beta_old, y, X)
  
  Ginv_new = hessian(beta_prop, y, X)
  grad_new = grad(beta_prop, y, X)
  
  new_given_old = dmvnorm(as.numeric(beta_prop), mean= beta_old + .5 * epsilon^2 * Ginv_old %*% grad_old,
                          sigma = epsilon^2 * Ginv_old, log = T)
  
  old_given_new = dmvnorm(as.numeric(beta_old), mean= beta_prop + .5 * epsilon^2 * Ginv_new %*% grad_new,
                          sigma = epsilon^2 * Ginv_new, log = T)
  
  b = exp(old_given_new - new_given_old)
  return(a * b)
}


# random walk metropolis
RWM_logRegr = function(Niter, y, X, epsilon){
  p=length(X[1,]); n=length(y);
  Res = matrix(NA,ncol=p,nrow = Niter);
  beta = double(p); #initialize with zeros
  lpi = logpi(beta,y,X);
  for ( jj in 1:Niter){
    beta_prop =  beta + epsilon*rnorm(p); #new beta
    lpi_prop = logpi(beta_prop,y,X);
    Acc = min(1,exp(lpi_prop-lpi))
    if (runif(1)<=Acc){
      beta= beta_prop;
      lpi=lpi_prop;
    }
    epsilon = epsilon + (1/jj^0.7)*(Acc - 0.4);
    Res[jj,] = beta
  }
  return(Res)
}


# MALA as described here: 
# https://en.wikipedia.org/wiki/Metropolis-adjusted_Langevin_algorithm

MALA_logRegr = function(Niter, y, X, epsilon){
  p=length(X[1,]); n=length(y);
  Res = matrix(NA,ncol=p,nrow = Niter);
  beta = double(p); #initialize with zeros
  lpi = logpi(beta,y,X);
  for ( jj in 1:Niter){
    beta_prop =  beta_prop_MALA(beta, epsilon, y, X) 
    alpha = acc_MALA(beta_old = beta, beta_prop, y, X, epsilon)
    Acc = min(1,alpha)
    if (runif(1)<=Acc){
      beta= beta_prop;
    }
    epsilon = epsilon + (1/jj^0.7)*(Acc - 0.5);
    Res[jj,] = beta
  }
  return(Res)
}


mMALA_logRegr = function(Niter, y, X, epsilon){
  p=length(X[1,]); n=length(y);
  Res = matrix(NA,ncol=p,nrow = Niter);
  beta = double(p); #initialize with zeros
  lpi = logpi(beta,y,X);
  for ( jj in 1:Niter){
    beta_prop =  beta_prop_MMALA(beta, epsilon, y, X) 
    alpha = acc_MMALA(beta_old = beta, beta_prop, y, X, epsilon)
    Acc = min(1,alpha)
    if (runif(1)<=Acc){
      beta= beta_prop;
    }
    epsilon = epsilon + (1/jj^0.7)*(Acc - 0.5);
    Res[jj,] = beta
  }
  return(Res)
}

MMALA_logRegr = function(Niter, y, X, epsilon){
  p=length(X[1,]); n=length(y);
  Res = matrix(NA,ncol=p,nrow = Niter);
  beta = double(p); #initialize with zeros
  lpi = logpi(beta,y,X);
  for ( jj in 1:Niter){
    beta_prop =  beta_prop_M2MALA(beta, epsilon, y, X) 
    alpha = acc_MMALA(beta_old = beta, beta_prop, y, X, epsilon)
    Acc = min(1,alpha)
    if (runif(1)<=Acc){
      beta= beta_prop;
    }
    epsilon = epsilon + (1/jj^0.7)*(Acc - 0.5);
    Res[jj,] = beta
  }
  return(Res)
}


# random walk metropolis
Gam_logRegr = function(Niter, y, X){
  p=length(X[1,]); n=length(y);
  Res = matrix(NA,ncol=p,nrow = Niter);
  beta = double(p); #initialize with zeros
  lpi = logpi(beta,y,X);
  for ( jj in 1:Niter){
    mC = gen_moments(beta, y, X)
    
    beta_prop = as.numeric(rmvnorm(1, mean = mC$m, sigma = mC$C))
    
    mC_star = gen_moments(beta_prop, y, X)
    
    #beta_prop =  beta + epsilon*rnorm(p); #new beta
    #lpi_prop = logpi(beta_prop,y,X);
    alpha = acc_fun(beta, beta_prop, moments_t = mC, moments_star = mC_star, y, X)
    
    Acc = min(1,alpha)
    if (runif(1)<=Acc){
      beta= beta_prop;
      #lpi=lpi_prop;
    }
    #epsilon = epsilon + (1/jj^0.7)*(Acc - 0.5);
    Res[jj,] = beta
  }
  return(Res)
}
