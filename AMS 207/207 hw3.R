# Consider the data available as "birthweight" form the package "LearnBayes". Fit a linear regression that considers age and gender as explanatory variables for birth weight. Describe the posterior distribution of the regression parameters using a sample-based approach. Explore the predictive posterior distribution for the birth weight of children in the following four cases: (a) 36 week female/male; (b) 40 week female/male. Compare.

library(LearnBayes)
attach(birthweight)


#The objective is to fit a model that describes the variation in the weight of babies in terms of the covariates age and gender.

#EDA
hist(weight, main = "Distribution of Weight", col = 'aquamarine')
hist(age, main = 'Distribution of Age', col = 'darkturquoise')
table(sex)
#No skewness

plot(age, weight)
#Since the categorical variables gender takes only two values, we use the jitter function in R to jitter the horizontal location of the points so we can see any overlapping points.
plot(jitter(gender), weight)

#We write the regression model as:
#E(weight|x,theta) = beta_0 + beta_1*age + beta_2*gender

#As the covariate gender is categorical with two levels, it can be represented by binary indicators; in the datafile birthweight, gender is coded 0 (1) for male (female).

#sex <- as.factor(gender)
#birthweight <- data.frame(age, sex, weight)
#str(birthweight)
#model=lm(weight ~ age + sex)

#We first perform the traditional least-squares fit using the lm command.
fit=lm(weight ~ age + gender, x=T, y=T)
summary(fit)

par(mfrow=c(1,2))
plot(fit, which = 1:2)

#We see from the output that AGE is strong effect; babies with a larger gestational age tend to have larger weights. The GENDER appears to not be significant.

#Describe the posterior distribution of the regression parameters using a sample-based approach.
#Bayesian fit - posterior means and sds of beta

#The function blinreg is used to sample from the joint posterior distribution of β and σ. The inputs to this function are the vector of values of the response variable y, the design matrix of the linear regression fit X, and the number of simulations m. We used the optional arguments x = TRUE, and y = TRUE in the function lm so that the design matrix and response vector are available as components of the structure fit.

theta.sample = blinreg(fit$y, fit$x, 10000)

#Posterior distributions of regression parameters
apply(theta.sample$beta, 2, mean)
apply(theta.sample$beta, 2, sd)
mean(theta.sample$sigma)
sd(theta.sample$sigma)
#The function blinreg returns two components: beta is a matrix of simulated draws from the marginal posterior of β, where each row is a simulated draw, and sigma is a vector of simulated draws from the marginal posterior of σ.

#The following R commands construct density plots of the simulated posterior draws of the individual regression coefficients β1 and β2 and the error standard deviation σ:
par(mfrow=c(2,2))
hist(theta.sample$beta[,1], main = 'Intercept', xlab = expression(beta[0]))
hist(theta.sample$beta[,2], main = 'Age', xlab = expression(beta[1]))
hist(theta.sample$beta[,3], main = 'Gender', xlab = expression(beta[2]))
hist(theta.sample$sigma, main = 'Error SD', xlab = expression(sigma))     

par(mar=c(1,1,1,1))
par(mfrow=c(4,1))
plot(ts(theta.sample$beta[,1]))
plot(ts(theta.sample$beta[,2]))
plot(ts(theta.sample$beta[,3]))
plot(ts(theta.sample$sigma))
dev.off()


#We can summarize each individual parameter by computing the 5th, 50th, and 95th percentiles of each collection of simulated draws. 
apply(theta.sample$beta,2,quantile,c(.05,.5,.95))

quantile(theta.sample$sigma,c(.05,.5,.95))

#As expected, the posterior medians of the regression parameters are similar in value to the ordinary regression estimates. Any small differences between the posterior medians and the least-squares estimates are due to small errors inherent in the simulation.

#predictive posterior distribution for the birth weight of children in the following four cases: (a) 36 week female/male; (b) 40 week female/male. Compare.


#We define the four sets of covariates and stack these sets in the matrix X1
cov1 = c(1, 36, 0)
cov2 = c(1, 36, 1)
cov3 = c(1, 40, 0)
cov4 = c(1, 40, 1)
X1=rbind(cov1,cov2,cov3,cov4)

#The function blinregpred will produce a simulated sample of future response values for a regression model. The inputs to the function blinregpred are a matrix X1 where each row corre- sponds to a covariate set and the structure of simulated values of the parameters β and σ. 
pred.response = blinregpred(X1,theta.sample)
par(mar=c(1,1,1,1))
par(mfrow=c(4,1))
plot(ts(pred.response[,1]), main="A")
plot(ts(pred.response[,2]), main="B")
plot(ts(pred.response[,3]), main="C")
plot(ts(pred.response[,4]), main="D")
dev.off()
apply(pred.response,2,quantile,c(.05, .5, .95))
mean(pred.response[,4])
sd(pred.response[,4])

#Displays histograms of the simulated draws from the predictive distribution for the same four sets of covariates.
c.labels=c("A","B","C","D")
par(mfrow=c(2,2))
for (j in 1:4) hist(pred.response[,j], main=paste("Covariate set",c.labels[j]),xlab="Weight")

#Posterior predictive check: checking if the observations are consistent with the fitted model

#Let yi∗ denote the density of a future weight for a baby with covariate vector xi. Using the function binregpred, we can simulate draws of the posterior predictive distributions for all y_1∗,...,y_24* by using fit$x as an argument.

#In the R code, we summarize each predictive distribution by the 5th and 95th quantiles and graph these distributions as line plots using the matplot command. We place the actual weights y1,...,y24 as solid dots in the figure. We are looking to see if the observed response values are consistent with the corresponding predictive distributions; any points that fall outside of the corresponding 90% interval band are possible outliers. There is one point (labeled in the figure) that exceed the 95th percentile, corresponding to a female baby age 36.

pred.draws=blinregpred(fit$x, theta.sample)
pred.sum=apply(pred.draws, 2, quantile, c(.05,.95))
par(mfrow=c(1,1))
ind=1:length(weight)
matplot(rbind(ind,ind), pred.sum, type="l", lty=1,col=1, xlab="Observation",ylab="Weight")
points(ind, weight, pch=19)
out = (weight > pred.sum[2,])
text(ind[out], weight[out], label=ind[out], pos = 4)


dim(X1)
array(0, c(length(theta.sample$sigma), 4))
