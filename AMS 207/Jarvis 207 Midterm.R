library(invgamma)
time <- log(volcano$Intereventtime[1:41])
hist(time)
plot(ts(volcano$Intereventtime))
plot(density(time))


#Gibbs Sampler
nu <- 5
n <- length(time)
a <- 10
b <- 10
c <- 10
d <- 10
s2 <- 10
m <- 10

sigmasq <- tausq <- mu <- rep(NA, 11000)
lambdai=matrix(0,11000,41)
lambdai[1,]=rgamma(n,0,1000)

mui=matrix(0,11000,41)
mui[1,]=rnorm(n,0,1000)

T <- 1000    # burnin
 # initialisation
sigmasq[1] <- 1
tausq[1] <- 1
mu[1] <- 1
for(i in 2:11000)
{   
  
  lambdai[i,] <- rgamma(n = 41, (nu+1)/2, (nu/2)+((time - mui[i-1,])^2/(2*sigmasq[i-1]))) 
  
  sigmasq[i] <- rinvgamma(n = 1, a + (n/2), b + ((lambdai[i,]*(time-mui[i,])^2)/2))
  
  tausq[i] <- rinvgamma(n = 1, c + (n/2), d + 0.5*sum((mui[i-1,]-mu[i-1])^2))
  
  mu[i] <- rnorm(n = 1, ((s2*n*mean(mui[i-1,])+tausq[i]*m)/(n*s2+tausq[i]))
                 , sqrt((s2*tausq[i])/(n*s2+tausq[i])))
  
  mui[i,] <- rnorm(n = 41, (tausq[i]*time*lambdai[i,]+
                             sigmasq[i]*mu[i])/(tausq[i]*lambdai[i,]+sigmasq[i])
                  , sqrt((tausq[i]*sigmasq[i])/(tausq[i]*lambdai[i,]+sigmasq[i])))
}

# remove burnin
mu <- mu[-(1:T)]
tau2 <- tausq[-(1:T)]
sigma2 <- sigmasq[-(1:T)]
mui <- mui[-(1:T),]
lambdai <- lambdai[-(1:T),]

#Gibbs Plots
par(mar=c(1,1,1,1))
par(mfrow=c(3,1))
plot(ts(mu))
plot(ts(tau2))
plot(ts(sigma2))
plot(ts(mui))
plot(ts(lambdai))
dev.off

#Distribution plots
plot(density(mu), main = expression('Distribution of'~mu), lwd=2, xlab = expression(mu))
polygon(density(mu), border = 'darkturquoise', col = 'darkslategray2', lwd=2)

plot(density(tau2), main = expression('Distribution of'~tau^2), lwd=2, xlab = expression(tau^2))
polygon(density(tau2), border = 'darkviolet', col = 'orchid1', lwd=2)

plot(density(sigma2), main = expression('Distribution of'~sigma^2), lwd=2, xlab = expression(sigma^2))
polygon(density(sigma2), border = 'deeppink', col = 'palevioletred1', lwd=2)

plot(density(mui), main = expression('Distribution of'~mu[i]), lwd=2, xlab = expression(mu[i]))
polygon(density(mui), border = 'springgreen3', col = 'palegreen', lwd=2)

plot(density(lambdai), main = expression('Distribution of'~lambda[i]), lwd=2, xlab = expression(lambda[i]))
polygon(density(lambdai), border = 'sienna2', col = 'tan1', lwd=2)

## Model 2 ##

sig2 <- mu2 <- rep(NA, 11000)

T <- 1000    # burnin
# initialisation
sig2[1] <- 1
mu2[1] <- 1
for(i in 2:11000)
{
  mu2[i] <- rnorm(n = 1, (s2*sum(time)+sig2[i-1]*m)/(n*s2+sig2[i-1]), (s2*sig2[i-1])/(n*s2+sig2))
  sig2[i] <- rinvgamma(n = 1, a+(n/2), b+0.5*sum((time-mu2[i])^2))
}
#remove burnin
mu2 <- mu2[-(1:T)]
sig2 <- sig2[-(1:T)]

plot(density(mu2), main = expression('Distribution of'~mu), lwd=2, xlab = expression(mu))
polygon(density(mu2), border = 'darkturquoise', col = 'darkslategray2', lwd=2)

plot(density(sig2), main = expression('Distribution of'~sigma^2), lwd=2, xlab = expression(sigma^2))
polygon(density(sig2), border = 'deeppink', col = 'palevioletred1', lwd=2)


### Compare the models ###

### DIC M1 ###
deviance = matrix(0, 10000, n)
for (i in 1:n)
  deviance[,i] = dnorm(time[i], mui[,i], sqrt(sigma2 / lambdai[,i]), log = TRUE)
deviance = -2*apply(deviance, 1, sum)
print(m1.DIC <- mean(deviance) + 0.5*var(deviance))

### DIC M2 ###
deviance2 = matrix(0, 10000, n)
for (i in 1:n)
  deviance2[,i] = dnorm(time[i], mu2, sqrt(sig2), log = TRUE)
deviance2 = -2*apply(deviance2, 1, sum)
print(m2.DIC <- mean(deviance2) + 0.5*var(deviance2))


### Bayes Factor with MCMC ###
## Model 1 ##
N <- 1000000 #Need a large n b/c law of large numbers
prior.mean = rnorm(N, rnorm(N, m, sqrt(s2)), rinvgamma(N, c, d))
prior.var = rinvgamma(N, a, b)/rgamma(N, nu/2, nu/2)
likeatprior <- matrix(0, N, n)
for (i in n) 
  {
  likeatprior[,i] = dnorm(time[i], prior.mean, sqrt(prior.var), log = T)
  }
bf1 <- apply(likeatprior, 1, sum)


## Model 2 ##
prior.mean.m2 <- rnorm(N, m, sqrt(s2))
prior.var.m2 <- rinvgamma(N, a, b)
likeatprior.2 <- matrix(0, N, n)
for(i in n)
{
  likeatprior.2[,i] = dnorm(time[i], prior.mean.m2, sqrt(prior.var.m2), log = T)
}
bf2 <- apply(likeatprior.2, 1, sum)

mean(exp(bf1))/mean(exp(bf2))
mean(exp(bf2))/mean(exp(bf1))
  
### G & G ###

# Model 1 #
predicted.y.m1 = matrix(0, 10000, n)
for(i in 1:n)
{
  predicted.y.m1[,i]= rnorm(10000, mui[,i], sqrt(sigma2/lambdai[,i]))
}
print(gg.m1 <- sum((time - apply(predicted.y.m1,2,mean))^2) + sum(apply(predicted.y.m1, 2, var)))


# Model 2 #
predicted.y.m2 = matrix(0, 10000, n)
for(i in 1:n)
{
  predicted.y.m2[,i]= rnorm(10000, mu2[i], sqrt(sig2[i]))
}
print(gg.m2 <- sum((time - apply(predicted.y.m2,2,mean))^2) + sum(apply(predicted.y.m2, 2, var)))


### Posterior predictive checking ###

hist(time, main = "Interevent Times", xlab = 'log(time)', col = 'aliceblue')
plot(density(time),  main = "Interevent Times", xlab = 'log(time)')
polygon(density(time), col = 'aliceblue')

#Each row of predicted should look like data
plot(density(predicted.y.m1[1,]), main = 'Posterior Predictive Sample Distribution', xlab = 'Sample')
polygon(density(predicted.y.m1[1,]), col = 'aliceblue')
hist(predicted.y.m1[1,], col = 'aliceblue', main = 'Posterior Predictive Sample Distribution', xlab = 'Sample')

plot(density(predicted.y.m2[5,]), main = 'Posterior Predictive Sample Distribution', xlab = 'Sample')
polygon(density(predicted.y.m2[5,]), col = 'aliceblue')

hist(predicted.y.m2[5,], col = 'aliceblue', main = 'Posterior Predictive Sample Distribution', xlab = 'Sample')


## Model 1 ##
par(mfrow=c(2, 2))
predict.mean = apply(predicted.y.m1, 1, mean)
hist(predict.mean, main = 'Posterior Predicted Means', xlab = "Predicted Mean")
abline(v=mean(time), lwd=2, col='darkturquoise')
#p-value
length(which(predict.mean>mean(time))==TRUE)/10000
text(4.75,2500, labels = 'p=0.5344')

predict.sd = apply(predicted.y.m1, 1, sd)
hist(predict.sd, main = 'Posterior Predicted SDs', xlab = "Predicted SD")
abline(v=sd(time), lwd=2, col='deeppink')
#p-value
length(which(predict.sd>sd(time))==TRUE)/10000
text(3,2300, labels = 'p=0.7837', cex = 0.8)

predict.min = apply(predicted.y.m1, 1, min)
hist(predict.min, main = "Posterior Predicted Minimums", xlab='Predicted Min')
abline(v=min(time), lwd=2, col='darkviolet')
#p-value
length(which(predict.min>min(time))==TRUE)/10000
text(-8,3000, labels = 'p=0.6443')

predict.max = apply(predicted.y.m1, 1, max)
hist(predict.max, main = "Posterior Predicted Maximums", xlab='Predicted Max')
abline(v=max(time), lwd=2, col='springgreen3')
#p-value
length(which(predict.max>max(time))==TRUE)/10000
text(19,3000, labels = 'p=0.4604')

dev.off()

### Model 2 ###
par(mfrow=c(2, 2))
predict.mean = apply(predicted.y.m2, 1, mean)
hist(predict.mean, main = 'Posterior Predicted Means', xlab = "Predicted Mean")
abline(v=mean(time), lwd=2, col='darkturquoise')
#p-value
length(which(predict.mean>mean(time))==TRUE)/10000
text(5,1500, labels = 'p=0.519')


predict.sd = apply(predicted.y.m2, 1, sd)
hist(predict.sd, main = 'Posterior Predicted SDs', xlab = "Predicted SD")
abline(v=sd(time), lwd=2, col='deeppink')
#p-value
length(which(predict.sd>sd(time))==TRUE)/10000
text(0.85,1200, labels = 'p=0.3104')

predict.min = apply(predicted.y.m2, 1, min)
hist(predict.min, main = "Posterior Predicted Minimums", xlab='Predicted Min')
abline(v=min(time), lwd=2, col='darkviolet')
#p-value
length(which(predict.min>min(time))==TRUE)/10000
text(0.85,2500, labels = 'p=0.9321')

predict.max = apply(predicted.y.m2, 1, max)
hist(predict.max, main = "Posterior Predicted Maximums", xlab='Predicted Max')
abline(v=max(time), lwd=2, col='springgreen3')
#p-value
length(which(predict.max>max(time))==TRUE)/10000
text(10.5,2500, labels = 'p=0.1398')
dev.off()
