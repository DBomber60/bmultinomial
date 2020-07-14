rm(list=ls())
set.seed(1)

n = 1000
p = 5
X = cbind(matrix(rnorm(n*p), nrow = n))
beta.true = c(1, 1, 2, 3, 0)
y = X%*%beta.true + rnorm(n)
m0 = lm(y ~ X - 1)
summary(m0)

solve ( t(X) %*% X ) %*% t(X) %*% y


lambda = 1

solve ( t(X) %*% X + lambda * diag(p)) %*% t(X) %*% y

# now collect results
lambda_seq = seq(from = 0, to = 50000, by = 100)
res = matrix(nrow = p, ncol = length(lambda_seq))
colnames(res) = paste("lambda",lambda_seq)
a = data.frame(res)
a$var = factor(1:5) #paste("beta",1:(p+1))

# now, fill in the data

for (m in 1:length(lambda_seq)) {
  a[,m] = solve ( t(X) %*% X + lambda_seq[m] * diag(p)) %*% t(X) %*% y
}

long = pivot_longer(a, cols = starts_with("lambda"), values_to = "beta")
long$nu = rep(lambda_seq, p)

ggplot(data = long, aes(x=nu, y=beta[,1])) +
  geom_line(aes(colour=var)) +
  #scale_color_manual(values=c( rep('#E69F00', 6), rep('#999999', 5)) ) +
  theme_bw() + ylab("lambda")

# now, for logistic regression
source('emvs.R')
set.seed(1)

n = 1000
p = 5
X = cbind(matrix(rnorm(n*p), nrow = n))
beta.true = c(1, 1, 2, 3, 0)
Y = rbinom(n,1, prob = plogis(X%*%beta.true))

#lambda = 0

#beta.hat = nr.beta(beta.init = rep(0,p), phat.init = rep(.5, n), dstar = lambda * rep(1,p))

# now collect results
lambda_seq = seq(from = 0, to = 30, by = 1)
res = matrix(nrow = p, ncol = length(lambda_seq))
colnames(res) = paste("lambda",lambda_seq)
a = data.frame(res)
a$var = factor(1:5) #paste("beta",1:(p+1))

# now, fill in the data

for (m in 1:length(lambda_seq)) {
  a[,m] = nr.beta(beta.init = rep(0,p), phat.init = rep(.5, n), dstar = lambda_seq[m] * rep(1,p))
}

long = pivot_longer(a, cols = starts_with("lambda"), values_to = "beta")
long$nu = rep(lambda_seq, p)

ggplot(data = long, aes(x=nu, y=beta[,1])) +
  geom_line(aes(colour=var)) +
  #scale_color_manual(values=c( rep('#E69F00', 6), rep('#999999', 5)) ) +
  theme_bw() + ylab("lambda")





