set.seed(3)
library(mvtnorm)

# compute the log-posterior dist. -- assumes N(0,I) prior on beta vector
logpi = function(beta,y,X, dstar){
  eta = X%*%beta;
  return(sum(y*eta -log(1+exp(eta))) - .5 * sum( dstar * beta^2 ) )
}

# compute gradient based on a N(0,I) prior on beta
grad = function(beta,y,X, dstar) {
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

# mMALA acceptance probability
acc_MMALA = function(beta_old, beta_prop, y, X, epsilon, dstar) {
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

mMALA_logRegr = function(Niter, y, X, epsilon){
  
  p=length(X[1,]); n=length(y);
  Res = matrix(NA,ncol=p,nrow = Niter);
  beta = double(p); #initialize with zeros
  lpi = logpi(beta,y,X, dstar);
  
  for (jj in 1:Niter) {
    beta_prop =  beta_prop_MMALA(beta, epsilon, y, X, dstar) 
    alpha = acc_MMALA(beta_old = beta, beta_prop, y, X, epsilon, dstar)
    Acc = min(1,alpha)
    if (runif(1)<=Acc){
      beta= beta_prop;
    }
    epsilon = epsilon + (1/jj^0.7)*(Acc - 0.5);
    Res[jj,] = beta
  }
  return(Res)
}

