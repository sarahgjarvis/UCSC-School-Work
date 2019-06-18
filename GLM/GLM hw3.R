#GLM hw 3

#Problem 1

#xi
logdose <- c(1.6907, 1.7242, 1.7552, 1.7842, 1.8113, 1.8369, 1.8610, 1.8839)
#mi
n.beetles <- c(59, 60, 62, 56, 63, 59, 62, 60)
#yi 
n.killed <- c(6, 13, 18, 28, 52, 53, 61, 60)

# (a) Fit 3 binomial GLMs to these data corresponding to 3 link functions, logit, probit and complementary log-log. 

prob <- n.killed/n.beetles

#logit
model1 <- glm(prob ~ logdose, weights = n.beetles, family = binomial(link = 'logit'))
summary(model1)

#probit
model2 <- glm(prob ~ logdose, weights = n.beetles, family = binomial(link = 'probit'))
summary(model2)

#cloglog
model3 <- glm(prob ~ logdose, weights = n.beetles, family = binomial(link = 'cloglog'))
summary(model3)

# Perform residual analysis for each model, using the deviance residuals.

dev_resid <- lapply(models,function(x) resid(x,type="deviance"))
sapply(dev_resid,function(x) sqrt(mean(x^2)))

d.resid.m1 <- resid(model1, type = "deviance")
d.resid.m2 <- resid(model2, type = "deviance")
d.resid.m3 <- resid(model3, type = "deviance")

d.resid <- rbind(d.resid.m1, d.resid.m2, d.resid.m3)

#RMS deviance residuals
apply(d.resid, 1, function(x) sqrt(mean(x^2)))

#Obtain fitted values, ˆπi, under each model and compare with observed proportions, yi/mi

# predict pi
pi.m1 <- predict(model1, type = 'response')
pi.m2 <- predict(model2, type = 'response')
pi.m3 <- predict(model3, type = 'response')

#Obtain the estimated dose-response curve under each model by evaluating ˆπ(x) = g−(βˆ1 +βˆ2x) over a grid of values x for log dose in the interval (1.65, 1.9). 

plot(logdose, prob, lwd = 2, main = 'Observed vs. Fitted Values', ylab = expression(y[i]/m[i]))
points(logdose, pi.m1, col = 'deeppink', lwd = 2)
points(logdose, pi.m2, col = 'darkmagenta', lwd = 2)
points(logdose, pi.m3, col = 'darkturquoise', lwd = 2)
legend('topleft', legend = c('observed','logit', 'probit', 'c log-log'), col = c('black','deeppink','darkmagenta', 'darkturquoise'), lwd = 2, cex = 0.8)

#grid predictions
grid <- seq(1.65, 1.9, length = 500)

grid.m1 <- predict(model1, newdata = list(logdose = grid), type = 'response')
grid.m2 <- predict(model2, newdata = list(logdose = grid), type = 'response')
grid.m3 <- predict(model3, newdata = list(logdose = grid), type = 'response')

plot(grid, grid.m1, type = 'l', col = 'deeppink', lwd = 2, main = 'Estimated Dose-Response Curves', ylab = expression(ginv(beta[1] + beta[2]*x)), xlab = 'logdose')
points(grid, grid.m2, type = 'l', col = 'darkmagenta', lwd = 2)
points(grid, grid.m3, type = 'l', col = 'darkturquoise', lwd = 2)
points(logdose, prob, type = 'l', lwd = 2)
legend('topleft', legend = c('observed', 'logit', 'probit', 'c log-log'), col = c('black','deeppink','darkmagenta', 'darkturquoise'), lwd = 2, cex = 0.8)


# (c)

special.invg <- function(x, a, b){
  eta <- c(x %*% b)
  return(exp(a*eta)/(1+exp(eta))^a)
}

beta.hat <- c(-113.625, 62.5)
alpha.hat <- 0.279

#estimated dose-response curve under the link
pred.special.invg <- special.invg(cbind(1, grid), alpha.hat, beta.hat)


plot(grid, grid.m1, type = 'l', col = 'deeppink', lwd = 2, main = 'Estimated Dose-Response Curves', ylab = expression(ginv(beta[1] + beta[2]*x)), xlab = 'logdose')
points(grid, grid.m2, type = 'l', col = 'darkmagenta', lwd = 2)
points(grid, grid.m3, type = 'l', col = 'darkturquoise', lwd = 2)
points(logdose, prob, type = 'l', lwd = 2)
points(grid, pred.special.invg, type = 'l', lwd = 2, col = 'plum')

legend('topleft', legend = c('observed', 'logit', 'probit', 'c log-log', 'modified logit'), col = c('black','deeppink','darkmagenta', 'darkturquoise', 'plum'), lwd = 2, cex = 0.8)

#obtain the fitted values pi hat
pi.hat <- special.invg(cbind(1, logdose), alpha.hat, beta.hat)
m <- n.beetles
y <- n.killed
y.hat <- pi.hat * m
n <- length(y)
res <- sign(y-y.hat) * sqrt(2*abs(y*log(y/y.hat) + (m-y)*log(1E-10+(m-y)/(m-y.hat))))
plot(res, ylab = "Deviance Residuals", main = "Residual Plot")
abline(h=0, col = "red")

# part (d)
logchoose <- function(m,y) lgamma(m+1) - lgamma(y+1) - lgamma(m-y+1)

K <- 3 # number of parameters

model1$aic
model2$aic
BIC(model1)
BIC(model2)
model3$aic
BIC(model3)
#modified logit AIC
aic <- -2 * sum(logchoose(m,y) + y*log(pi.hat) + (m-y)*log(1-pi.hat)) + 2*K
#BIC
aic + 3 * (log(n)-2)

#Problem 2

#Bayesian regression

#Consider a Bayesian binomial GLM with a complementary log-log link
invcloglog <- function(y, beta, x){
  return(1-exp(-exp(beta[1]+beta[2]*x)))
}

y = n.killed
m = n.beetles
N = length(y)
jinv <- 2*vcov(model3)

library(MASS)

#flat prior 
prior <- function(x){
  return(1)
}


prior <- function(x){
  return(dnorm(x[1], 0, 1, log = T) + dnorm(x[2], 0, 1, log = T))
}

k <- function(beta, m, x, y){
  return(exp(prior(beta) + sum(y * log(1-exp(-exp(beta[1] + beta[2]*x))) + (m-y)*-exp(beta[1] + beta[2]*x))))
}

beta <- matrix(NA, ncol = 2, nrow = 10000)
beta[1,] <- c(0, 0)
for (i in 2:10000){
  beta.squig <- mvrnorm(n = 1, beta[i-1,], jinv)
  q <- min(1, k(beta.squig, m, logdose, y)/k(beta[i-1,], m, logdose, y))
  if(runif(1) < q){
    beta[i,] <- beta.squig
  }
  else{
    beta[i,] <- beta[i-1,]
  }
}


#beta.flat <- beta
#beta.stnorm <- beta
#beta.norm <- beta

m.cloglog <- apply(beta, 2, mean)

par(mfrow = c(1,2))

plot(ts(beta.flat[,1]), ylim = c(-50, 0), ylab = expression(beta[1]))
lines(ts(beta.stnorm[,1]), col = 'turquoise')
lines(ts(beta.norm[,1]), col = 'magenta')

plot(ts(beta.flat[,2]), ylab = expression(beta[2]))
lines(ts(beta.stnorm[,2]), col = 'turquoise')
lines(ts(beta.norm[,2]), col = 'magenta')

par(mfrow = c(2,1))
plot(density(beta.flat[5000:10000,1]), xlim = c(-55, 0), ylim = c(0, 0.6), lwd = 2, main = "Parameter Distributions", xlab = expression(beta[1]))
lines(density(beta.stnorm[5000:10000,1]), col = 'turquoise', lwd = 2)
lines(density(beta.norm[5000:10000,1]), col = 'magenta', lwd = 2)
legend('topleft', legend = c(expression(p(beta[1]) ~ "~1"), expression(p(beta[1]) ~ "~N(0, 1)"), expression(p(beta[1]) ~ "~N(0,100)")), cex = 0.8, col = c('black', 'turquoise', 'magenta'), lwd = 2)

plot(density(beta.flat[5000:10000,2]), xlim = c(0, 30), ylim = c(0, 0.8), lwd = 2, xlab = expression(beta[2]), main = "")
lines(density(beta.stnorm[5000:10000,2]), col = 'turquoise', lwd = 2)
lines(density(beta.norm[5000:10000,2]), col = 'magenta', lwd = 2)
legend('topright', legend = c(expression(p(beta[2]) ~ "~1"), expression(p(beta[2]) ~ "~N(0, 1)"), expression(p(beta[2]) ~ "~N(0,100)")), cex = 0.8, col = c('black', 'turquoise', 'magenta'), lwd = 2)

q <- matrix(NA, 500, 3)
x <- seq(1.65, 1.9, length = 500)
for (i in 1:500){
  dist <- 1-exp(-exp(beta.flat[,1] + beta.flat[,2]*x[i]))
  q[i,] <- quantile(dist, probs = c(0.05, 0.5, 0.95))
}


dist.flat <- q
plot(x, dist.flat[,1], type = "l", main = "Dose-Response Curve (Flat Prior)", xlab = "logdose", ylab = expression(pi(x)), lty = 2, col = 'darkturquoise', lwd = 2)
lines(x, dist.flat[,2], lwd = 2)
lines(x, dist.flat[,3], lty = 2, col = 'darkturquoise', lwd = 2)
#legend('topleft', legend = c('point estimate', '90% interval estimates'), col = c('black', 'darkturquoise'), lwd = 2, cex = 0.8)


#obtain the posterior distribution for the median lethal dose, LD50, that is, the dose level at which the probability of response is 0.5

x.post <- (-0.3665 - beta.flat[500:10000,1])/beta.flat[500:10000,2]
plot(density(x.post), main = "Posterior Distribution for the Median Lethal dose", xlab = "logdose")
mean(x.post)


# Part b - new pi function - logit link
jinv <- 2*vcov(model1)

prior <- function(beta){
  return(1)
}

k2 <- function(beta, m, x, y){
  return(exp(log(prior(beta)) + sum(y*(beta[1]+beta[2]*x) - y*log(1 + exp(beta[1] + beta[2]*x)) - (m - y)*log(1 + exp(beta[1] + beta[2]*x)))))
}


#obtain MCMC samples from the posterior distributions
beta <- matrix(NA, ncol = 2, nrow = 10000)
beta[1,] <- c(0, 0)
for (i in 2:10000){
  beta.squig <- mvrnorm(n = 1, beta[i-1,], jinv)
  q <- min(1, k2(beta.squig, m, logdose, y)/k2(beta[i-1,], m, logdose, y))
  if(runif(1) < q){
    beta[i,] <- beta.squig
  }
  else{
    beta[i,] <- beta[i-1,]
  }
}

#beta.logit <- beta

m.logit <- apply(beta.logit, 2, mean)


plot(ts(beta.logit), main = 'MCMC Samples logit Link')

par(mar=c(4,1,1,1))
par(mfrow = c(2,1))
plot(density(beta.logit[,1]), main = "Posterior Density Estimates", xlab = expression(beta[1]))
abline(v=mean(beta.logit[,1]), col = 'red')
legend('topright', legend = 'mean = -16.25', col = 'red', lwd = 1, bty = 'n')

plot(density(beta.logit[,2]), main = "", xlab = expression(beta[2]))
abline(v=mean(beta.logit[,2]), col = 'red')
legend('topleft', legend = 'mean = 34.56', col = 'red', lwd = 1, bty = 'n')

#lethal dose 0.5
x.post <- beta.logit[500:10000,1]/beta.logit[500:10000,2]
plot(density(x.post), main = "Posterior Distribution for the Median Lethal dose", xlab = "logdose")
mean(x.post)

#point and interval estimates for the dose-response curve π(x)
q <- matrix(NA, 500, 3)
x <- seq(1.65, 1.9, length = 500)
for (i in 1:500){
  dist <- exp(beta.logit[,1]+beta.logit[,2]*x[i])/((1 + exp(beta.logit[,1] + beta.logit[,1]*x[i])))
  q[i,] <- quantile(dist, probs = c(0.05, 0.5, 0.95))
}


dist.logit <- q
plot(x, dist.logit[,1], type = "l", main = "Dose-Response Curve (logit)", xlab = "logdose", ylab = expression(pi(x)), lty = 2, col = 'darkturquoise', lwd = 2)
lines(x, dist.logit[,2], lwd = 2)
lines(x, dist.logit[,3], lty = 2, col = 'darkturquoise', lwd = 2)
#legend('topleft', legend = c('point estimate', '90% interval estimates'), col = c('black', 'darkturquoise'), lwd = 2, cex = 0.8)

#Part (c)
#Consider the binomial GLM with the parametric link given in (1.1). Develop an MCMC method to sample from the posterior distribution of (β1, β2, α).

#expand the data
## Number of Steps for MCMC
N = 100000

## Two vectors storing the beta coefficients
b1 <- numeric()
b2 <- numeric()
alpha <- numeric()

## Set starting values for b1 and b2
b1[1] = 0
b2[1] = 0
alpha[1] = 1


###################################################################
#### This is the part we modify (i.e. logit to modified logit) ####
loglikelihood <- function(b1, b2, alpha){
  ## Calculate p through link function
  p = exp(alpha*(b1+b2*dosage))/((1 + exp(b1 + b2*dosage))^alpha)
  # print(p)
  ## sum dbinom at each level of dosage
  sum(dbinom(dead,1, p, log=TRUE))
}
###################################################################

## Generating N draws from random uniform for acceptance/rejection
## This is more computationally efficient than drawing from runif in every realization


## Running the algorithm for N realizations
for (i in 2:N) {
  ## Draw from the candidate, which is a normal dist. centered around previously accepted value
  ## SDs were approximated using trial and error (they were tuned)
  b1.prime <- rnorm(1, b1[i-1], 1)  
  b2.prime <- rnorm(1, b2[i-1], 1)
  alpha.prime <- rgamma(1, alpha[i-1], 2)
  
  ## Priors
  ## I didn't have prior knowledge of beetle deaths, so I used a flat prior
  ## i.e. a normal distribution with a very large standard deviation
  current.prior.dens = dnorm(b1.prime,0,100000,log = TRUE) + dnorm(b2.prime,0,100000,log=TRUE) + dgamma(alpha.prime, 1,1)
  prime.prior.dens = dnorm(b1[i-1],0,100000, log = TRUE) + dnorm(b2[i-1],0,100000, log = TRUE) + dgamma(alpha[i-1], 1,1)
  
  ## Calculating log likelihood for both the proposed and current values of b1 and b2
  current.loglik = loglikelihood(b1[i-1], b2[i-1], alpha[i-1])
  prime.loglik = loglikelihood(b1.prime, b2.prime, alpha.prime)
  
  ## Posterior
  current.target.dens = current.loglik + current.prior.dens
  prime.target.dens = prime.loglik + prime.prior.dens
  
  ## Alpha denotes the acceptance probability 
  acceptance =  exp(prime.target.dens - current.target.dens)
  
  ## Comparing alpha against a draw from U[0,1]
  if (runif(1) < acceptance) {
    ## Updating the chain if proposed value is accepted
    b1[i] = b1.prime
    b2[i] = b2.prime
    alpha[i] = alpha.prime
  } else {
    ## Updating the chain if chain proposed value is not accepted (chain doesn't move)
    b1[i] = b1[i-1]
    b2[i] = b2[i-1]
    alpha[i] = alpha[i-1]
  }
  print(i)
}

#remove burnin
b1.new <- b1[-(1:60000)]
b2.new <- b2[-(1:60000)]
alpha.new <- alpha[-(1:60000)]


par(mar=c(2,4,1,1))
par(mfrow = c(3,1))
plot(ts(b1.new), ylab = expression(beta[1]))
plot(ts(b2.new), ylab = expression(beta[2]))
plot(ts(alpha.new), ylab = expression(alpha))

m.b1 <- mean(b1.new)
m.b2 <- mean(b2.new)
m.a <- mean(alpha.new)

par(mar=c(4,1,1,1))
plot(density(b1.new), main = "Parameter Distributions", xlab = expression(beta[1]))
plot(density(b2.new), main = "", xlab = expression(beta[2]))
plot(density(alpha), main = "", xlab = expression(alpha))


#LD50 again
blah <- function(b1, b2, alpha){
  return((-log(exp(b1.new) * 2^(1/alpha.new) - exp(b1.new)))/b2.new)
}
plot(density(blah(b1.new, b2.new, alpha.new)), main = expression(LD[50]), xlab = 'logdose')

#point and interval estimates for pi

q <- matrix(NA, 500, 3)
x <- seq(1.65, 1.9, length = 500)
for (i in 1:500){
  dist <- (exp(b1.new + b2.new*x[i])/(1 + exp(b1.new + b2.new*x[i])))^alpha.new
  q[i,] <- quantile(dist, probs = c(0.05, 0.5, 0.95))
}

plot(x, q[,1], type = "l", main = "Dose-Response Curve (modified logit)", xlab = "logdose", ylab = expression(pi(x)), lty = 2, col = 'darkturquoise', lwd = 2)
lines(x, q[,2], lwd = 2)
lines(x, q[,3], lty = 2, col = 'darkturquoise', lwd = 2)

#perform a residual analysis for each model using the Bayesian residuals:
#(yi/mi) − π(xi; β1, β2) for the first two models, and (yi/mi) − π(xi; β1, β2, α) for the third

post.m1 <- (1-exp(-exp(m.cloglog[1] + m.cloglog[2]*logdose)))
post.m2 <- (exp(m.logit[1]+m.logit[2]*logdose)/((1 + exp(m.logit[1] + m.logit[1]*logdose))))
post.m3 <- (exp(m.b1 + m.b2*logdose)/(1 + exp(m.b1 + m.b2*logdose)))^m.a

resid.m1 <- (y/m) - post.m1
resid.m2 <- (y/m) - post.m2
resid.m3 <- (y/m) - post.m3

par(mar=c(2,4,1,1))
par(mfrow = c(3,1))
plot(resid.m1, main = "Bayesian Residuals", lwd = 2)
plot(resid.m2, lwd = 2)
plot(resid.m3, lwd = 2)

#use the quadratic loss L measure for formal comparison of the three models

loss <- function(data, postpred, K){
  return(sum(var(postpred) + K/(K+1) * (data - mean(postpred,2))^2))
}

loss(logdose, post.m1, 2)
loss(logdose, post.m2, 2)
loss(logdose, post.m3, 3)

