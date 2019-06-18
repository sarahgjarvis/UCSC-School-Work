data <- data.frame(matrix(c(551, 6,
                            651, 4,
                            832, 17,
                            375, 9,
                            715, 14,
                            868, 8,
                            271, 5, 
                            630, 7, 
                            491, 7,
                            372, 7,
                            645, 6,
                            441, 8,
                            895, 28,
                            458, 4,
                            642, 10,
                            492, 4,
                            543, 8, 
                            842, 9, 
                            905, 23, 
                            542, 9, 
                            522, 6, 
                            122, 1, 
                            657, 9, 
                            170, 4, 
                            738, 9, 
                            371, 14, 
                            735, 17, 
                            749, 10, 
                            495, 7, 
                            716, 3, 
                            952, 9, 
                            417, 2), ncol= 2, byrow=TRUE))
names(data) <- c("length", "faults")

#fit a Bayesian Poisson GLM with the logarithmic link, log(µi) = β1 + β2xi

#to check results
model <- glm(faults ~ length, data = data, family = poisson(link = 'log'))

library(MASS)

prior <- function(x){
  return(log(1))
} 

k <- function(beta){
  return(exp(prior(beta) + sum(data$faults*(beta[1] + beta[2]*data$length) - exp(beta[1] + beta[2]*data$length))))
}

jinv <- vcov(model)

#Metropolis-Hastings Algorithm for beta1 and beta2
beta <- matrix(NA, ncol = 2, nrow = 10000)
#initialize at MLE
beta[1,] <- model$coefficients
for (i in 2:10000){
  beta.squig <- mvrnorm(n = 1, beta[i-1,], jinv)
  q <- min(1, k(beta.squig)/k(beta[i-1,]))
  if(runif(1) < q){
    beta[i,] <- beta.squig
  }
  else{
    beta[i,] <- beta[i-1,]
  }
  print(i)
}


mean(beta[,1])
mean(beta[,2])
#0.9704817, 0.00193221 YAY
par(mar=c(2,4,1,1))
par(mfrow = c(2,1))
plot(ts(beta[,1]), ylab = expression(beta[1]))
plot(ts(beta[,2]), ylab = expression(beta[2]))

#posterior distributions
par(mar=c(4,4,1,1))
par(mfrow = c(2,1))
plot(density(beta[,1]), main='Posterior Distributions', xlab=expression(beta[1]))
plot(density(beta[,2]), main='', xlab=expression(beta[2]))

#response mean as a function of the covariate
q <- matrix(NA, 500, 3)
x <- seq(100, 1000, length = 500)
for (i in 1:500){
  dist <- exp(beta[,1] + beta[,2]*x[i])
  q[i,] <- quantile(dist, probs = c(0.05, 0.5, 0.95))
}

plot(x, q[,2], type = 'l', ylim = c(0,20), ylab = 'faults', xlab = 'length', lwd=2, lty=3, main="Response Mean")
points(data$length, data$faults, pch=20, lwd=2)
polygon(c(rev(x), x), 
        c(rev(q[ ,3]), q[ ,1]), 
        col =rgb(0,1,1,alpha=0.3), border = NA)

#posterior predictive distribution
post.pred <- apply(beta, 1, function(theta){rpois(32, exp(theta[1] + theta[2]*data$length))})

#plot some posterior predictive samples
plot(density(data$faults), main = 'Distribution of Faults', xlab='faults')
for(i in 1:5){lines(density(post.pred[,i]), col='mediumturquoise', main=NULL)}
lines(density(data$faults), col='black', main="Data Distribution")
legend('topright', legend = c('predictive sample', 'data'), cex=0.8, col = c('mediumturquoise', 'black'), lty=1)

#plot of prediction intervals with data superimposed
pred.int <- apply(post.pred, 1, quantile, c(.05,.95))
ind <- 1:32
matplot(rbind(ind,ind), pred.int, type = 'l', lty = 1, xlab="Observation",ylab="Faults", col='darkturquoise', main='Posterior Prediction Intervals')
points(ind, data$faults, pch=19, cex=0.5)

#90% interval for residuals
resids <- matrix(NA, 2, 32)
for (i in 1:32){
  resids[1,i] <- quantile(data$faults[i] - post.pred[i,], 0.05)
  resids[2,i] <- quantile(data$faults[i] - post.pred[i,], 0.95)
}

matplot(rbind(data$length, data$length), resids, type = "l", lty = 1, xlab = 'length', ylab = 'residuals', col = 'darkturquoise', main = 'Posterior Predictive Residual Distributions')
abline(h=0)


#PART B 

#Gibbs sampler
y <- data$faults
x <- data$length
n <- nrow(data)
N <- 51000
T <- 1000 #burn in
t <- 1.01 #for prior sensitivity test

k1 <- function(mui, lambda, beta){
  return(sum(-lambda*mui/exp(beta[1] + beta[2]*x) - lambda*(beta[1] + beta[2]*x)))
}

k2 <- function(mui, lambda, beta){
  return(sum(lambda*log(mui) - lambda*mui/exp(beta[1] + beta[2]*x) + lambda*log(lambda) - lambda*(beta[1]+beta[2]*x)-lgamma(lambda))-t*log(1+lambda) + log(lambda))
}

beta <- matrix(0, N, 2)
mui <- gammai <- matrix(0, N, n)
lambda <- rep(NA, N)

#initialize
mui[1,] <- rnorm(n, 0, 1000)
beta[1,] <- model$coefficients
lambda[1] <- 8

for(i in 2:N)
  {

  mui[i,] <- rgamma(n, y + lambda[i-1], (lambda[i-1]/exp(beta[i-1,1] + beta[i-1,2]*x)) + 1)
  
  #M-H step for beta
  beta.squig <- mvrnorm(n = 1, beta[i-1,], jinv)
  q <- min(1, exp(k1(mui[i,], lambda[i-1], beta.squig) - k1(mui[i,], lambda[i-1], beta[i-1,])))
  if(runif(1) < q){
    beta[i,] <- beta.squig
  }
  else{
    beta[i,] <- beta[i-1,]
  }
  
  #M-H step for lambda
  lambda.squig <- exp(rnorm(n = 1, log(lambda[i-1]), 1))
  q <- min(1, exp(k2(mui[i,], lambda.squig, beta[i-1,]) - k2(mui[i,], lambda[i-1], beta[i-1,])))
  if(runif(1) < q){
    lambda[i] <- lambda.squig
  }
  else{
    lambda[i] <- lambda[i-1]
  }
}

# remove burnin
mui <- mui[-(1:T),]
beta <- beta[-(1:T),]
lambda <- lambda[-(1:T)]

par(mar=c(2,4,1,1))
par(mfrow = c(2,1))
plot(ts(beta[,1]), xlab='iteration', ylab=expression(beta[1]))
plot(ts(beta[,2]), xlab='iteration', ylab=expression(beta[2]))

plot(ts(lambda))
plot(ts(lambda), xlab='iteration', ylab=expression(lambda))
plot(ts(mui[,1]), xlab='iteration', ylab=expression(mu[1]))

colMeans(beta)
colMeans(mui)
mean(lambda)
median(lambda)

plot(density(beta[,1]), main ='Distribution of' ~ beta[1]) 
abline(v=mean(beta[,1]), col='red')
plot(density(beta[,2]), main = 'Distribution of' ~ beta[2])
abline(v=mean(beta[,2]), col='red')

plot(density(lambda), main ='Distribution of' ~ lambda) 
abline(v=median(lambda), col='red')
plot(density(mui[,1]), main = 'Distribution of' ~ mu[1])
abline(v=mean(mui[,1]), col='red')

#response mean as a function of the covariate
q <- matrix(NA, 500, 3)
x <- seq(100, 1000, length = 500)
for (i in 1:500){
  dist <- exp(beta[,1] + beta[,2]*x[i])
  q[i,] <- quantile(dist, probs = c(0.05, 0.5, 0.95))
}

plot(x, q[,2], type = 'l', ylim = c(0,20), ylab = 'faults', xlab = 'length', lwd=2, lty=3, main="Response Mean")
points(data$length, data$faults, pch=20, lwd=2)
polygon(c(rev(x), x), 
        c(rev(q[ ,3]), q[ ,1]), 
        col =rgb(0.5,0.25,1,alpha=0.3), border = NA)

#posterior predictive distribution for observed data
post.pred.b <- apply(mui, 1, function(mui){rpois(32, mui)})

#model checking
#plot posterior predictive samples
plot(density(data$faults), main = 'Distribution of Faults', xlab='faults')
for(i in 1:5){lines(density(post.pred.b[,i]), col='blueviolet', main=NULL)}
lines(density(data$faults), col='black', main="Data Distribution")
legend('topright', legend = c('predictive sample', 'data'), cex=0.8, col = c('blueviolet', 'black'), lty=1)

#plot of prediction intervals with data superimposed
pred.int <- apply(post.pred.b, 1, quantile, c(.05,.95))
ind <- 1:32
matplot(rbind(ind,ind), pred.int, type = 'l', lty = 1, xlab="Observation",ylab="Faults", col='blueviolet', main='Posterior Prediction Intervals')
points(ind, data$faults, pch=19, cex=0.5)

#90% interval for residuals
resids <- matrix(NA, 2, 32)
for (i in 1:32){
  resids[1,i] <- quantile(data$faults[i] - post.pred.b[i,], 0.05)
  resids[2,i] <- quantile(data$faults[i] - post.pred.b[i,], 0.95)
}

matplot(rbind(data$length, data$length), resids, type = "l", lty = 1, xlab = 'length', ylab = 'residuals', col = 'blueviolet', main = 'Posterior Predictive Residual Distributions')
abline(h=0)

#Gelfand and Gosh
print(gg.m1 <- sum((y - apply(post.pred,1,mean))^2) + sum(apply(post.pred, 1, var)))
#953.7096

print(gg.m2 <- sum((y - apply(post.pred.b,1,mean))^2) + sum(apply(post.pred.b, 1, var)))
#577.2691

print(gg.m3 <- sum((y - apply(post.pred.b,1,mean))^2) + sum(apply(post.pred.b, 1, var)))
#560.6615

print(gg.m4 <- sum((y - apply(post.pred.b,1,mean))^2) + sum(apply(post.pred.b, 1, var)))
#540.5264

print(gg.m5 <- sum((y - apply(post.pred.b,1,mean))^2) + sum(apply(post.pred.b, 1, var)))
#600.9412