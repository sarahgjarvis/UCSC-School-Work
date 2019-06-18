
# The real one this time 

library(LearnBayes)

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

n_burn = 1000
thin = 10

#Distribution plots
par(mfrow=c(2,2))

plot(density((sample$mu[seq(n_burn+1, ns,by=thin)])), main = expression('Distribution of' ~ mu), xlab=expression(mu))

plot(density(sqrt(sample$sig2y[seq(n_burn+1, ns,by=thin)])), main = expression('Distribution of' ~ sigma[y]), xlab=expression(sigma[y]))

plot(density((sqrt(sample$sig2b[seq(n_burn+1, ns,by=thin)]))), main = expression('Distribution of' ~ sigma[b]), xlab=expression(sigma[b]))

plot(density(sample$b[1, seq(n_burn+1, ns, by=thin)]), main = expression('Distribution of' ~ b[i]), xlab=expression(b[i]))

dev.off()

#posterior mean for the parameters
mu.mean=mean(sample$mu[seq(n_burn+1, ns,by=thin)])
sig2y.mean=mean(sample$sig2y[seq(n_burn+1, ns,by=thin)])
sig2b.mean=mean(sample$sig2b[seq(n_burn+1, ns,by=thin)])
sigy.mean=mean(sqrt(sample$sig2y[seq(n_burn+1, ns,by=thin)]))
sigb.mean=mean(sqrt(sample$sig2b[seq(n_burn+1, ns,by=thin)]))
mean(sample$b[,seq(n_burn+1, ns,by=thin)])
#SD
mu.sd=sd(sample$mu[seq(n_bur+1, ns,by=thin)])
sig2y.sd=sd(sample$sig2y[seq(n_bur+1, ns,by=thin)])
sig2b.sd=sd(sample$sig2b[seq(n_bur+1, ns,by=thin)])
sigy.sd=sd(sqrt(sample$sig2y[seq(n_bur+1, ns,by=thin)]))
sigb.sd=sd(sqrt(sample$sig2b[seq(n_bur+1, ns,by=thin)]))

#credibility intervals
quantile(sample$mu[seq(n_burn+1, ns,by=thin)], probs=c(0.025, 0.975))
quantile(sample$sig2y[seq(n_burn+1, ns,by=thin)], probs=c(0.025, 0.975))
quantile(sample$sig2b[seq(n_burn+1, ns,by=thin)], probs=c(0.025, 0.975))
quantile(sqrt(sample$sig2y[seq(n_burn+1, ns,by=thin)]), probs=c(0.025, 0.975))
quantile(sqrt(sample$sig2b[seq(n_burn+1, ns,by=thin)]), probs=c(0.025, 0.975))
quantile(sample$b[,seq(n_burn+1, ns,by=thin)], probs=c(0.025, 0.975))

apply(sample$b[,seq(n_burn+1, ns,by=thin)], 1, quantile, probs=c(.025,.5, .975), na.rm=TRUE)

par(mfrow=c(3,1))
plot.ts(sample$mu[seq(n_burn+1, ns,by=thin)], ylab=expression(mu))
plot.ts(sample$sig2y[seq(n_burn+1, ns,by=thin)], ylab=expression(sigma[y]^2))
plot.ts(sample$sig2b[seq(n_burn+1, ns,by=thin)], ylab=expression(sigma[b]^2))

#------------(4)------------

y=matrix(c(1545,1440,1440,1520,1580,
              1540,1555,1490,1560,1495,
              1595,1550,1605,1510,1560,
              1445,1440,1595,1465,1545,
              1595,1630,1515,1635,1625,
              1520,1455,1450,1480,1445),c(6,5),byrow=TRUE)

#marginal posterior density
log.post.var.comp=function(theta,y)
{
  mu = theta[1]
  sigma.y = exp(theta[2])
  sigma.b = exp(theta[3])
  Y=apply(y,1,mean)
  n=dim(y)[2]
  S=apply(y,1,var)*(n-1)
  loglike=sum(dnorm(Y,mu,sqrt(sigma.y^2/n+sigma.b^2),log=TRUE)+dgamma(S,shape=(n-1)/2,rate=1/(2*sigma.y^2),log=TRUE))
  return(loglike+theta[2]+theta[3])
}

##################Finding the posterior mode###################

fit=laplace(log.post.var.comp, c(1500, 10, 10), y)
theta.hat.mle<-fit$mode
se.hat.theta.hat <- sqrt( diag( fit$var ) )

############Plotting the log posterior density############

n.grid <- 25

c( theta.hat.mle[ 1 ] - 3 * se.hat.theta.hat[ 1 ], 
   theta.hat.mle[ 1 ] + 3 * se.hat.theta.hat[ 1 ] )

mu.grid <- seq( theta.hat.mle[ 1 ] - 3 * se.hat.theta.hat[ 1 ], 
                theta.hat.mle[ 1 ] + 3 * se.hat.theta.hat[ 1 ], length = n.grid )

c( theta.hat.mle[ 2 ] - 3 * se.hat.theta.hat[ 2 ], 
   theta.hat.mle[ 2 ] + 3 * se.hat.theta.hat[ 2 ] )


sig.y.grid <- seq( theta.hat.mle[ 2 ] - 3 * se.hat.theta.hat[ 2 ], 
                   theta.hat.mle[ 2 ] + 3 * se.hat.theta.hat[ 2 ], length = n.grid )

c( theta.hat.mle[ 3 ] - 3 * se.hat.theta.hat[ 3 ], 
   theta.hat.mle[ 3 ] + 3 * se.hat.theta.hat[ 3 ] )


sig.b.grid <- seq( theta.hat.mle[ 3 ] - 3 * se.hat.theta.hat[ 3 ], 
                theta.hat.mle[ 3 ] + 3 * se.hat.theta.hat[ 3 ], length = n.grid )

###############Fixed sigma y###################
log.likelihood.grid <- matrix( NA, n.grid, n.grid )

for ( i in 1:n.grid ) {
  
  for ( j in 1:n.grid ) {
    
    log.likelihood.grid[ i, j ] <- 
      log.post.var.comp( c( mu.grid[ i ],
                                      theta.hat.mle[2], sig.b.grid[ j ]), y )
    
  }
  
}
par( mfrow = c( 2, 2 ) )

contour( mu.grid, sig.b.grid, log.likelihood.grid,
         xlab = 'mu', ylab = 'log(sigma[b])', nlevels = 10, col = 'blue',
         lwd = 2 )

persp( mu.grid, sig.b.grid, log.likelihood.grid,
       xlab = 'mu', ylab = 'log(sigma[b])', theta = 30, phi = 30,
       zlab = 'log likelihood', col = 'red' )

contour( mu.grid, sig.b.grid, exp( log.likelihood.grid ),
         xlab = 'mu', ylab = 'log(sigma[b])', nlevels = 10, col = 'blue',
         lwd = 2 )

persp( mu.grid, sig.b.grid, exp( log.likelihood.grid ),
       xlab = 'mu', ylab = 'log(sigma[b])', theta = 30, phi = 30,
       zlab = 'likelihood', col = 'red' )

par( mfrow = c( 1, 1 ) )

dev.off( ) 

############Fixed sig.b#################
log.likelihood.grid <- matrix( NA, n.grid, n.grid )

for ( i in 1:n.grid ) {
  
  for ( j in 1:n.grid ) {
    
    log.likelihood.grid[ i, j ] <- 
      log.post.var.comp( c( mu.grid[ i ],
                            sig.y.grid[ j ], theta.hat.mle[3]), y )
    
  }
  
}
par( mfrow = c( 2, 2 ) )

contour( mu.grid, sig.y.grid, log.likelihood.grid,
         xlab = 'mu', ylab = 'log(sigma[y])', nlevels = 10, col = 'blue',
         lwd = 2 )

persp( mu.grid, sig.y.grid, log.likelihood.grid,
       xlab = 'mu', ylab = 'log(sigma[y])', theta = 30, phi = 30,
       zlab = 'log likelihood', col = 'red' )

contour( mu.grid, sig.y.grid, exp( log.likelihood.grid ),
         xlab = 'mu', ylab = 'log(sigma[y])', nlevels = 10, col = 'blue',
         lwd = 2 )

persp( mu.grid, sig.y.grid, exp( log.likelihood.grid ),
       xlab = 'mu', ylab = 'log(sigma[y])', theta = 30, phi = 30,
       zlab = 'likelihood', col = 'red' )

par( mfrow = c( 1, 1 ) )

dev.off( ) 

#####################Fixed Mu######################
log.likelihood.grid <- matrix( NA, n.grid, n.grid )

for ( i in 1:n.grid ) {
  
  for ( j in 1:n.grid ) {
    
    log.likelihood.grid[ i, j ] <- 
      log.post.var.comp( c(theta.hat.mle[1], sig.y.grid[ i ],
                            sig.b.grid[ j ]), y )
    
  }
  
}
par( mfrow = c( 2, 2 ) )

contour( sig.y.grid, sig.b.grid, log.likelihood.grid,
         xlab = 'log(sigma[y])', ylab = 'log(sigma[b])', nlevels = 10, col = 'blue',
         lwd = 2 )

persp( sig.y.grid, sig.b.grid, log.likelihood.grid,
       xlab = 'log(sigma[y])', ylab = 'log(sigma[b])', theta = 30, phi = 30,
       zlab = 'log likelihood', col = 'red' )

contour( sig.y.grid, sig.b.grid, exp( log.likelihood.grid ),
         xlab = 'log(sigma[y])', ylab = 'log(sigma[b])', nlevels = 10, col = 'blue',
         lwd = 2 )

persp( sig.y.grid, sig.b.grid, exp( log.likelihood.grid ),
       xlab = 'log(sigma[y])', ylab = 'log(sigma[b])', theta = 30, phi = 30,
       zlab = 'likelihood', col = 'red' )

par( mfrow = c( 1, 1 ) )

dev.off( ) 


#Gibbs values to compare
c(mean(sample$mu[seq(n_bur+1, ns,by=thin)]), log(mean(sqrt(sample$sig2y[seq(n_bur+1, ns,by=thin)]))), log(mean(sqrt(sample$sig2b[seq(n_bur+1, ns,by=thin)]))))

##################Normal approximation####################

###mu and sigy###
contour( mu.grid, sig.y.grid, log.likelihood.grid,
         xlab = 'mu', ylab = 'log(sigma[y])', nlevels = 10,
         lwd = 2 )
npar=list(m=fit$mode[1:2],v=fit$var[c(1,2), c(1,2)])
mycontour(lbinorm,c(1450,1630,0,6),npar,xlab=expression(mu), ylab=expression(log(sigma[y])), add=T, main = "Normal Approximation to the Posterior Distribution", col='hotpink')

###Mu and sigb###
contour( mu.grid, sig.b.grid, log.likelihood.grid,
         xlab = 'mu', ylab = 'log(sigma[b])', nlevels = 10,
         lwd = 2 )
qpar=list(m=cbind(fit$mode[1], fit$mode[3]),v=fit$var[c(1,3), c(1,3)])
mycontour(lbinorm,c(1450,1630,0,6),qpar,xlab=expression(mu), ylab=expression(log(sigma[b])),main = "Normal Approximation to the Posterior Distribution", add=T, col='hotpink')

###sig and sig###
contour( sig.y.grid, sig.b.grid, log.likelihood.grid,
         xlab = 'log(sigma[y])', ylab = 'log(sigma[b])', nlevels = 10,
         lwd = 2 )
zpar=list(m=fit$mode[2:3],v=fit$var[c(2,3), c(2,3)])
mycontour(lbinorm,c(2,6,2,6),zpar,xlab=expression(log(sigma[y])), ylab=expression(log(sigma[b])),main = "Normal Approximation to the Posterior Distribution", col='hotpink', add=T)

#construct 95% interval estimates for mu, log sigma.y, and log sigma.b
sds=sqrt(diag(fit$var))[1:3]
modes=fit$mode[1:3]
int.ests=cbind(modes-1.96*sds, modes+1.96*sds)

sds=sqrt(diag(fit$var))[2:3]
modes=fit$mode[2:3]
int.ests=cbind(modes-1.96*sds, modes+1.96*sds)
exp(int.ests)

#-------(5)----------
#Rejection Sampling

#Changing to multivariate t
logfT=function(theta,y)
{
  data=datapar$data
  tpar=datapar$par
  d=log.post.var.comp(theta,y)-dmt(theta,mean=c(tpar$m),
                              S=tpar$var,df=tpar$df,log=TRUE)
  return(d)
}
tpar=list(m=fit$mode,var=2*fit$var,df=4)
datapar=list(data=y,par=tpar)
lapt=laplace(logfT,c(1500, 10, 10),y)
tmode=lapt$mode
tmax=logfT(tmode, y)


library(hexbin)
library(RColorBrewer)
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)

###########Sampling Importance Re-Sampling###############
sir=sir(log.post.var.comp, tpar, 10000, y)

#####mu and sigy#####
contour( mu.grid, sig.y.grid, log.likelihood.grid,
         xlab = 'mu', ylab = 'log(sigma[y])', nlevels = 10,
         lwd = 2 )
points(sir[,1],sir[,2])
plot(hexbin(sir[,1],sir[,2]), colramp=rf, xlab=expression(mu),
     ylab=expression(log(sigma[y])), main="Sampling Importance Resampling")

###mu and sigb####
contour( mu.grid, sig.b.grid, log.likelihood.grid,
         xlab = 'mu', ylab = 'log(sigma[b])', nlevels = 10,
         lwd = 2, main='Sampling Importance Resampling')
points(sir[,1], sir[,3])
plot(hexbin(sir[,1], sir[,3]), colramp=rf, xlab=expression(mu),
     ylab=expression(log(sigma[b])), main="Sampling Importance Resampling")

####sig and sig####
contour( sig.y.grid, sig.b.grid, log.likelihood.grid,
         xlab = 'log(sigma[y])', ylab = 'log(sigma[b])', nlevels = 10,
         lwd = 2 )
points(sir[,2],sir[,3])
plot(hexbin(sir[,2],sir[,3]), colramp=rf, xlab=expression(log(sigma[y])),
     ylab=expression(log(sigma[b])), main="Sampling Importance Resampling")

#To compare with results from Gibbs
plot(hexbin(sample$mu,log(sqrt(sample$sig2y))), colramp=rf, xlab=expression(mu),
     ylab=expression(log(sigma[y])), main="Gibbs Sampling")


plot(hexbin(sample$mu,log(sqrt(sample$sig2b))), colramp=rf, xlab=expression(mu),
     ylab=expression(log(sigma[y])), main="Gibbs Sampling")


plot(hexbin(log(sqrt(sample$sig2y)),log(sqrt(sample$sig2b))), colramp=rf, xlab=expression(log(sigma[y])), ylab=expression(log(sigma[b])), main="Gibbs Sampling")

#Results
#yield (mu)
mean(sir[,1])
plot(density(sir[,1]), main='Distribution of Yield', xlab='Yield')
legend(1340, .015, legend = 'mean=15.27.19',lty = 1, col = 'hotpink', lwd=2, cex = 0.8)
abline(v=mean(sir[,1]), col='hotpink', lwd=2)

#Yield variability
change<-exp(sir[,2])
plot(density(change), main='Distribution of Yield Variability', xlab=expression(sigma[y]))
abline(v=median(change), col='hotpink', lwd=2)     
legend(80, .04, legend = 'median=51.51',lty = 1, col = 'hotpink', lwd=2, cex = 0.8)
median(change)

#Batches variability (sigb)
change<-exp(sir[,3])
plot(density(change), main='Distribution of Batches Variability', xlab=expression(sigma[b]))
abline(v=median(change), col='hotpink', lwd=2)     
legend(200, .015, legend = 'median=50.12',lty = 1, col = 'hotpink', lwd=2, cex = 0.8)
median(change)

dev.off()
plot(density(y), main='Data Distribution', xlab='Yield')
