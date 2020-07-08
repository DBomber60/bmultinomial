# spike and slab regularization plot
theta_0 = 0.5
nu_0_seq = seq(from = 0.1, to = 2, by = .1)
res = matrix(nrow = p+1, ncol = length(nu_0_seq))
colnames(res) = paste("nu0",nu_0_seq)
a = data.frame(res)
a$var = factor(1:11) #paste("beta",1:(p+1))

# now, fill in the data

for (m in 1:length(nu_0_seq)) {
  em.iter = EMVS(beta_0, theta_0, nu_0_seq[m], nu_1=1000)
  print(em.iter[[3]] )
  a[,m] = em.iter[[1]]
}

long = pivot_longer(a, cols = starts_with("nu0"), values_to = "beta")
long$nu = rep(nu_0_seq, p+1)

ggplot(data = long, aes(x=nu, y=beta[,1])) +
  geom_line(aes(colour=var)) +
  scale_color_manual(values=c( rep('#E69F00', 6), rep('#999999', 5)) ) +
  theme_bw()

