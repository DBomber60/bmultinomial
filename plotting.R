library(tidyverse) # pivot_longer
library(latex2exp)

# spike and slab regularization plot

# input: nu_0 sequence, initial beta and theta; p, number of predictors (including intercept)
# output: regularization plot

plot_emvs = function(p, nsig, beta_0, theta.init, nu0_low, nu0_high, increment, nu_1) {
  
  # sequence of nu_0 values
  nu_0_seq = seq(from = nu0_low, to = nu0_high, by = increment)
  
  # dataframe to hold results
  res = matrix(nrow = p, ncol = length(nu_0_seq))
  colnames(res) = paste("nu0",nu_0_seq)
  res_df = data.frame(res)
  res_df$var = factor(1:p) #paste("beta",1:(p+1))
  
  # fill in 
  for (m in 1:length(nu_0_seq)) {
    em.iter = EMVS(beta_0, theta.init, nu_0_seq[m], nu_1 = nu_1)
    res_df[,m] = em.iter[[1]]
  }
  
  long = pivot_longer(res_df, cols = starts_with("nu0"), values_to = "beta")
  long$nu = rep(nu_0_seq, p)
  
  ggplot(data = long, aes(x=nu, y=beta[,1])) +
    geom_line(aes(colour=var)) +
    scale_color_manual(values=c( rep('#E69F00', 6), rep('#999999', (p-5))) ) +
    theme(legend.position = "none") + ylab(TeX("$\\hat{\\beta}$")) +
    xlab(TeX("$\\nu_0$")) + theme_bw()
}

# now, fill in the data



