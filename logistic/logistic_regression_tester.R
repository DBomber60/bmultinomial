# for 6-5-2020 LC meeting
source('MALA.R')
#source('BayesLRsampler.R')

#Generate data
n=200;
#X=cbind(rep(1,n), matrix(rnorm((p-1)*n),ncol=p-1));

sig = matrix(c(1,.7,.7,1), nrow = 2, ncol = 2)
X=cbind(rep(1,n), rmvnorm(n, sigma = sig));

#X = rmvnorm(n, sigma = )

beta_true=c(-2,1,2);
#Generate the Y’s
mp = X%*%beta_true;
prob = exp(mp)/(1+exp(mp))
y = rbinom(n,1,prob)


par(mfrow=c(3,3))

#Call the MCMC Sampler
epsilon = .05;
Niter = 10000;
burnin = floor(0.2 * Niter)
whichBeta = 2

#ResRW = RWM_logRegr(Niter,y,X,epsilon = 0.07)
#ResMALA = MALA_logRegr(Niter,y,X,epsilon)
#ResmMALA = mMALA_logRegr(Niter,y,X,epsilon)
#ResGam = Gam_logRegr(Niter,y,X)
ResMM = MMALA_logRegr(Niter,y,X,epsilon)

plotter = function(output, whichbeta, burnin, Niter, type) {
  vec = output[burnin:Niter,whichBeta]
  plot(vec,type='l',col = rgb(red = 1, green = 0, blue = 0, alpha = 0.5), main = type)
  acf(vec,lag.max=100,col='blue', main = "")
  hist(vec,nclass=50,prob=T,col='blue', main = "")
}

# plot
par(mfrow=c(4,3))

#plotter(ResRW, whichbeta = 2, burnin, Niter, type = "Random Walk")
plotter(ResGam, whichbeta = 2, burnin, Niter, type = "Gamerman")
plotter(ResMALA, whichbeta = 2, burnin, Niter, type = "MALA")
plotter(ResmMALA, whichbeta = 2, burnin, Niter, type = "mMALA")
plotter(ResMM, whichbeta = 2, burnin, Niter, type = "MMALA")

summary(glm(y~X-1, family = binomial))
