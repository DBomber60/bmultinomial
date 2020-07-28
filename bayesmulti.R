library(nnet)
library(ks) # vec
set.seed(1)

# log likelihood
# input: Y (vector response); eta = X %*% beta.current: p x (k-1) matrix of coefficients

logpi_m = function(Y, eta) { # y_int, 
  n = nrow(eta)
  eta_ = cbind(rep(0,n), eta)
  sum ( rowSums(eta_ * Y) - log(rowSums(exp(eta_))) )
  #sum( eta_[cbind(1:n, y_int)] - log(rowSums(exp(eta_))) )
}

# gradient: vector of length (k-1) times p (p covariates for each of K-1 classes)
# dstar is the inverse of the prior covariance

grad_multi = function(B, phat_m, X, Y, K, w, dstar) {
  n = nrow(Y)
  -B * dstar + (kronecker(diag(K-1), t(X))) %*% array(vec((Y[,2:K] - w * phat_m)), dim = c(n*(K-1),1)) 
}

# hessian: matrix - kp times kp
H = function(phat_m, X, K, w, dstar) {
  p = dim(X)[2]
  n = dim(X)[1]
  t1 = array(0, dim = c( p * (K-1), p* (K-1) ) )
  for(i in 1:n) {
    t1 = t1 + kronecker( diag(w[i] * phat_m[i,]) - outer(w[i] * phat_m[i,], phat_m[i,]), outer(X[i,], X[i,])) #diag
  }
  return(solve(dstar + t1) ) # diag((K-1)*p) +
}

# input: two lists - old derivatives(d) and new (q)
alpha_fun = function(d, q, beta_old, beta_star, epsilon, p) {
  mu_old = beta_old + .5 * epsilon^2 * d$Ginv %*% d$grad  #* d$Ginv %*% d$grad
  mu_new = beta_star + .5 * epsilon^2 * q$Ginv %*% q$grad # q$Ginv %*% q$grad
  
  new_given_old = dmvnorm(as.numeric(beta_star), mean= mu_old, sigma = epsilon^2 * diag((K-1)*p), log = T)  # d$Ginv
  old_given_new = dmvnorm(as.numeric(beta_old), mean= mu_new, sigma = epsilon^2 * diag((K-1)*p), log = T) # q$Ginv
  
  alpha = exp(q$lpi + old_given_new - (d$lpi + new_given_old))
  
  return(alpha)
}

# input: beta (p x (k-1) matrix)
# output: lpi, gradient, Ginv, sqrt(Ginv)
get_derivs = function(Y, beta.current, eta, X) {
  #beta_mat = matrix(beta, nrow = p)
  #eta = X %*% beta_mat
  lpi = logpi_m(Y, eta)
  
  # phat matrix
  phat_m = exp(eta)/(1+rowSums(exp(eta)))
  grad = grad_multi(beta.current, phat_m, X, Y)
  
  Ginv = H(phat_m, X)
  rootG = chol(Ginv)
  
  return(list(lpi = lpi,
              grad = grad,
              Ginv = Ginv,
              rootG = rootG))  
}


sample_betas = function(Zcurr, beta.current, X, epsilon) {
  p = dim(X)[2]
  
  eta_current = (X %*% matrix(beta.current, ncol = (K-1)))
  
  # get loppi, gradient, hessian based on new beta
  d = get_derivs(Zcurr, beta.current, eta = eta_current, X)
  
  # now propose a new beta 
  beta_prop = beta.current + as.numeric(.5 * epsilon^2 * d$Ginv %*% d$grad + epsilon * d$rootG %*% rnorm((K-1)*p))
  
  # accept it?
  
  # first need derivatives at new location
  q = get_derivs(Zcurr, beta_prop, eta = X %*% matrix(beta_prop, ncol = (K-1)), X)
  
  alpha = alpha_fun(d = d, q = q, beta.current, beta_star = beta_prop, epsilon, p)
  
  Acc = min(1,alpha)
  if (runif(1)<=Acc){
    beta.current = beta_prop;
  }
  
  return(list(Acc = Acc, beta.current = beta.current, epsilon = epsilon))
}



# functions for higher order MALA

# creates dPdB_ij vector for i=1,...k-1 and j=1,...p ( (k-1)*p of them)
dP_fun = function(phat_matrix, subject_index, class_index, coef_index) {
  result = array(0, dim = c(k-1,1))
  p_i = phat_matrix[subject_index,]
  for (i in 1:(k-1)) {
    if(i==class_index) {result[i] = p_i[i]*(1-p_1[i])
    } else {result[i] = -p_i[class_index]*p_i[i] }
  }
  return(as.numeric(result * X[subject_index, coef_index] ))
}

# get the whole thing for class 1, coefficient one
dGdB = function(class_index, coef_index, phat_m) {
  res = array(0, dim = c((k-1)*p, (k-1)*p ))
  for (subj in 1:n) {
    dP = dP_fun(phat_m, subj, class_index, coef_index)
    t1 = outer(dP, phat_m[subj,]) + outer(phat_m[subj,], dP)
    res = res + kronecker(t1, outer(X[subj,], X[subj,]))
  }
  return(res)
}

# dP_fun(phat_matrix, 1, 1, 1)
# dP = c(p_1[1]*(1-p_1[1]), -p_1[1]*p_1[2]) * X[1,1]
# 
# dG3 = array(0, dim = rep( (k-1)*p, 3) )
# 
# 
# for (class in 1:(k-1)) {
#   for (coef in 1:p) {
#     dG3[,,(class-1)*p + coef] = dGdB(class, coef, phat_matrix)
#   }
# }

