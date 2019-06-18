
# ams 207 - take home 1

rm(list=ls(all=TRUE))
set.seed(7)
require(lme4)

data(Dyestuff)
Dyestuff$Batch = as.numeric(Dyestuff$Batch)

N = 6
n = 5
Nn = N*n

dyes = (as.data.frame(Dyestuff))
aggregate(dyes, by=list(dyes$Batch), FUN=mean)
aggregate(dyes, by=list(dyes$Batch), FUN=var)

# gibbs sampling

## functions to sample each parameter from its full conditionals

### update b_i
fn_update_b = function(y, N, n, mu, sig2y, sig2b)
{
  b<-NULL
  for(j in 1:N)
  {
    b[j] = rnorm(1, ((sum(dyes$Yield[dyes$Batch==j] - mu)/sig2y)*(1/((n/sig2y)+(1/sig2b)))), sqrt(1/((n/sig2y)+(1/sig2b))))
  }
  return(b)
}

## udpate mu
fn_update_mu = function(y_b, sig2y,Nn)
{
  va = sig2y/Nn
  mme = va*(y_b/sig2y)
  
  mu = rnorm(1, mme, sqrt(va))
  return(mu)
}

## update sig2y
fn_update_sig2y = function(Nn, y_mu_b_sq)
{
  sigy = 1/rgamma(1, ((Nn-2)/2), (y_mu_b_sq/2))  ## mean a/b
  return(sigy)
}

## update sig2b
fn_update_sig2b = function(N, y_b_sq)
{
  sig2b = 1.0/rgamma(1, ((N-2)/2), (y_b_sq/2))
  return(sig2b)
}

## initial values
par_sam = NULL
par_sam$b = rep(mean(dyes$Yield), N)
par_sam$mu = 1000
par_sam$sig2y = 10^3
par_sam$sig2b = 10^3

## variables for the MCMC
ns = 51000

## save simulated para
sample = NULL
sample$b = array(NA, dim=c(N, ns))
sample$mu = rep(NA, ns)
sample$sig2y = rep(NA, ns)
sample$sig2b = rep(NA, ns)

## MCMC modeling

for(i_iter in 1:ns)
{
  if((i_iter%%1000)==0)
  {
    print(paste("i.iter=", i_iter))
    print(date())
  }
  
  ## udpate theta
  par_sam$b = fn_update_b(dyes, N, n, par_sam$mu, par_sam$sig2y, par_sam$sig2b)
  
  for(i in 1:N){
    if(i==1){
      par_sam$b_list = rep(par_sam$b[i],n)
    } else{
      p_the = rep(par_sam$b[i],n)
      par_sam$b_list = c(par_sam$b_list,p_the)}
  }
  
  dyes$b_list = par_sam$b_list
  
  par_sam$mu = fn_update_mu((sum(dyes$Yield - dyes$b_list)), par_sam$sig2y, Nn)
  
  ## update sigma y
  par_sam$sig2y = fn_update_sig2y(Nn, sum((dyes$Yield - par_sam$mu - dyes$b_list)^2))
  
  ## udpate sigma b
  par_sam$sig2b = fn_update_sig2b(N, sum(par_sam$b^2))
  
  ## save cur_sam
  sample$b[,i_iter] = par_sam$b
  sample$mu[i_iter] = par_sam$mu
  sample$sig2y[i_iter] = par_sam$sig2y
  sample$sig2b[i_iter] = par_sam$sig2b
}

n_bur = 1000
thin = 10

#Distribution plots
plot(density(sample$b[1, seq(n_burn+1, ns, by=thin)]), main = expression('Distribution of' ~ b[i]))

plot(density((sample$sig2y[seq(n_bur+1, ns,by=thin)])), main = expression('Distribution of' ~ sigma[y]^2))

plot(density((sqrt(sample$sig2b[seq(n_bur+1, ns,by=thin)]))), main = expression('Distribution of' ~ sigma[b]^2))

plot(density((sample$mu[seq(n_bur+1, ns,by=thin)])), main = expression('Distribution of' ~ mu))


#posterior mean for the parameters
mean(sample$mu[seq(n_bur+1, ns,by=thin)])
mean(sample$sig2y[seq(n_bur+1, ns,by=thin)])
mean(sample$sig2b[seq(n_bur+1, ns,by=thin)])
mean(sqrt(sample$sig2y[seq(n_bur+1, ns,by=thin)]))
mean(sqrt(sample$sig2b[seq(n_bur+1, ns,by=thin)]))

apply(sample$b[,seq(n_bur+1, ns,by=thin)], 1, mean)

#credibility intervals
quantile(sample$mu[seq(n_bur+1, ns,by=thin)], probs=c(0.025, .5, 0.975))
quantile(sample$sig2y[seq(n_bur+1, ns,by=thin)], probs=c(0.025, .5, 0.975))
quantile(sample$sig2b[seq(n_bur+1, ns,by=thin)], probs=c(0.025, .5, 0.975))
quantile(sqrt(sample$sig2y[seq(n_bur+1, ns,by=thin)]), probs=c(0.025, .5, 0.975))
quantile(sqrt(sample$sig2b[seq(n_bur+1, ns,by=thin)]), probs=c(0.025, .5, 0.975))

apply(sample$b[,seq(n_bur+1, ns,by=thin)], 1, quantile, probs=c(.025,.5, .975), na.rm=TRUE)

# the crude densities for mu and sigma2
par(mfrow=c(1,3))
plot.ts(sample$b[3,seq(n_bur+1, ns,by=thin)])

plot.ts(sample$mu[seq(n_bur+1, ns,by=thin)])
plot.ts(sample$sig2y[seq(n_bur+1, ns,by=thin)])
plot.ts(sample$sig2b[seq(n_bur+1, ns,by=thin)])

par(mfrow=c(1,3))
plot(density(sample$mu[seq(n_bur+1, ns,by=thin)]), main = expression(paste(mu)), xlab="")
plot(density(sqrt(sample$sig2y[seq(n_bur+1, ns,by=thin)])), xlab="", main = expression(paste(sigma[y])))
plot(density(sqrt(sample$sig2b[seq(n_bur+1, ns,by=thin)])), xlab="", main = expression(paste(sigma[b])))

par(mfrow=c(3,1))
plot(density(sample$mu[seq(n_bur+1, ns,by=thin)]), cex.lab=1.5,cex.main = 3,cex.axis=2,main = expression(paste(mu)), xlab="")
plot(density(log(sqrt(sample$sig2y[seq(n_bur+1, ns,by=thin)]))), cex.lab=1.5,cex.main =3,cex.axis=2, xlab="", main = expression("log(" * sigma[y] * ")"))
plot(density(log(sqrt(sample$sig2b[seq(n_bur+1, ns,by=thin)]))), cex.lab=1.5,cex.main = 3,cex.axis=2, xlab="", main = expression("log(" * sigma[b] * ")"))


# graphs for densities for the mu and transformation log(sigma)
# just to compare with the results given by the normal approximation
# for the transformed parameters
par(mfrow=c(1,3))
plot.ts(sample$mu[seq(n_bur+1, ns,by=thin)])
plot.ts(log(sqrt(sample$sig2y[seq(n_bur+1, ns,by=thin)])))
plot.ts(log(sqrt(sample$sig2b[seq(n_bur+1, ns,by=thin)])))

par(mfrow=c(1,3))
plot(density(sample$mu[seq(n_bur+1, ns,by=thin)]), main = expression(paste(mu)))
plot(density(log(sqrt(sample$sig2y[seq(n_bur+1, ns,by=thin)]))), main = expression(paste(sigma[y])))
plot(density(log(sqrt(sample$sig2b[seq(n_bur+1, ns,by=thin)]))), main = expression(paste(sigma[b])))

