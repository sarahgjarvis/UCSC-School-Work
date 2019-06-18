#Problem 2: alligator food choice
data <- data.frame(matrix(c(1.30, 'I', 'M',
                    1.80, 'F', 'M', 
                    1.24, 'I', 'F', 
                    2.56, 'O', 'F',
                    1.32, 'F', 'M',
                    1.85, 'F', 'M',
                    1.30, 'I', 'F',
                    2.67, 'F', 'F',
                    1.32, 'F', 'M',
                    1.93, 'I', 'M',
                    1.45, 'I', 'F',
                    2.72, 'I', 'F',
                    1.40, 'F', 'M',
                    1.93, 'F', 'M',
                    1.45, 'O', 'F',
                    2.79, 'F', 'F',
                    1.42, 'I', 'M',
                    1.98, 'I', 'M',
                    1.55, 'I', 'F',
                    2.84, 'F', 'F',
                    1.42, 'F', 'M',
                    2.03, 'F', 'M',
                    1.60, 'I', 'F',
                    1.47, 'I', 'M',
                    2.03, 'F', 'M',
                    1.60, 'I', 'F',
                    1.47, 'F', 'M',
                    2.31, 'F', 'M',
                    1.65, 'F', 'F',
                    1.50, 'I', 'M',
                    2.36, 'F', 'M',
                    1.78, 'I', 'F',
                    1.52, 'I', 'M',
                    2.46, 'F', 'M',
                    1.78, 'O', 'F',
                    1.63, 'I', 'M',
                    3.25, 'O', 'M',
                    1.80, 'I', 'F',
                    1.65, 'O', 'M',
                    3.28, 'O', 'M',
                    1.88, 'I', 'F',
                    1.65, 'O', 'M',
                    3.33, 'F', 'M',
                    2.16, 'F', 'F',
                    1.65, 'I', 'M',
                    3.56, 'F', 'M',
                    2.26, 'F', 'F',
                    1.65, 'F', 'M',
                    3.58, 'F', 'M',
                    2.31, 'F', 'F',
                    1.68, 'F', 'M',
                    3.66, 'F', 'M',
                    2.36, 'F', 'F',
                    1.70, 'I', 'M',
                    3.68, 'O', 'M',
                    2.39, 'F', 'F',
                    1.73, 'O', 'M',
                    3.71, 'F', 'M',
                    2.41, 'F', 'F',
                    1.78, 'F', 'M',
                    3.89, 'F', 'M',
                    2.44, 'F', 'F',
                    1.78, 'O', 'M'), ncol = 3,byrow = T))
names(data) <- c('length', 'choice','gender')
data$length <- as.numeric(as.character(data$length))

#length is the single covariate
#Develop a Bayesian multinomial regression model, using the baseline-category logits formulation with “fish” as the baseline category, to estimate (with point and interval estimates) the response probabilities as a function of length.

female <- data[data$gender == 'F',]
female$length <- as.numeric(as.character(female$length))
male <- data[data$gender == "M",]
male$length <- as.numeric(as.character(male$length))

library(MCMCpack)

#bayesian model
m0 <- MCMCmnl(choice ~ length, baseline = 'F', mcmc.method="IndMH", data = data, burnin = 1000, mcmc = 20000)

# #model for female alligators
# m1 <- MCMCmnl(choice ~ length, baseline = 'F', mcmc.method="IndMH", data = female, burnin = 1000, mcmc = 20000)
# 
# #model for male alligators
# m2 <- MCMCmnl(choice ~ length, baseline = 'F', mcmc.method="IndMH", data = male, burnin = 1000, mcmc = 20000)

#model with gender as a covariate
m3 <- MCMCmnl(choice ~ length + gender, baseline = 'F', mcmc.method="IndMH", data = data, burnin = 1000, mcmc = 20000)

summary(m0)
# summary(m1)
# summary(m2)
summary(m3)

apply(m0, 2, quantile, probs=c(0.025, 0.5, 0.975))
# apply(m1, 2, quantile, probs=c(0.025, 0.5, 0.975))
# apply(m2, 2, quantile, probs=c(0.025, 0.5, 0.975))
apply(m3, 2, quantile, probs=c(0.025, 0.5, 0.975))

par(mar=c(4,1,1,1))
par(mfrow = c(3,2))
traceplot(m3[,"length.I"], main="Invertebrate Coefficient", ylab = expression(beta[1]))
traceplot(m3[,"(Intercept).I"], main="Invertebrate Intercept", ylab = expression(alpha[1]))
traceplot(m3[,"length.O"], main="Other Coefficient", ylab = expression(beta[2]))
traceplot(m3[,"(Intercept).O"], main="Other Intercept", ylab = expression(alpha[1]))
traceplot(m3[,"genderM.I"], main="Gender Invertebrate Coefficient")
traceplot(m3[,"genderM.O"], main="Gender Other Coefficient")

par(mar=c(4,4,1,1))
par(mfrow = c(4,1))
traceplot(m3[,"length.I"], main="Invertebrate Coefficient", ylab = expression(beta[1]))
traceplot(m3[,"(Intercept).I"], main="Invertebrate Intercept", ylab = expression(alpha[1]))
traceplot(m3[,"length.O"], main="Other Coefficient", ylab = expression(beta[2]))
traceplot(m3[,"(Intercept).O"], main="Other Intercept", ylab = expression(alpha[2]))


par(mar=c(4,4,1,1))
par(mfrow = c(4,2))

# traceplot(m1[,"length.I"], main="Invertebrate Coefficient", ylab = expression(beta[1]), col = "darkolivegreen")
# traceplot(m2[,"length.I"], main="Invertebrate Coefficient", ylab = expression(beta[1]), col = "darkmagenta")
# 
# traceplot(m1[,"(Intercept).I"], main="Invertebrate Intercept", ylab = expression(alpha[1]), col = "darkolivegreen")
# traceplot(m2[,"(Intercept).I"], main="Invertebrate Intercept", ylab = expression(alpha[1]), col = "darkmagenta")
# 
# traceplot(m1[,"length.O"], main="Other Coefficient", ylab = expression(beta[2]), col = "darkolivegreen")
# traceplot(m2[,"length.O"], main="Other Coefficient", ylab = expression(beta[2]), col = "darkmagenta")
# 
# traceplot(m1[,"(Intercept).O"], main="Other Intercept", ylab = expression(alpha[2]), col = "darkolivegreen")
# traceplot(m2[,"(Intercept).O"], main="Other Intercept", ylab = expression(alpha[2]), col = "darkmagenta")

beta1.mean <- mean(m0[,"length.I"])
beta2.mean <- mean(m0[,"length.O"])
alpha1.mean <- mean(m0[,"(Intercept).I"])
alpha2.mean <-mean(m0[,"(Intercept).O"])

fbeta1.mean <- mean(m1[,"length.I"])
fbeta2.mean <- mean(m1[,"length.O"])
falpha1.mean <- mean(m1[,"(Intercept).I"])
falpha2.mean <-mean(m1[,"(Intercept).O"])

mbeta1.mean <- mean(m2[,"length.I"])
mbeta2.mean <- mean(m2[,"length.O"])
malpha1.mean <- mean(m2[,"(Intercept).I"])
malpha2.mean <-mean(m2[,"(Intercept).O"])

#save values
beta1 <- m0[,"length.I"]
beta2 <- m0[,"length.O"]
alpha1 <- m0[,"(Intercept).I"]
alpha2 <- m0[,"(Intercept).O"]

fbeta1 <- m1[,"length.I"]
fbeta2 <- m1[,"length.O"]
falpha1 <- m1[,"(Intercept).I"]
falpha2 <- m1[,"(Intercept).O"]
mbeta1 <- m2[,"length.I"]
mbeta2 <- m2[,"length.O"]
malpha1 <- m2[,"(Intercept).I"]
malpha2 <- m2[,"(Intercept).O"]

gbeta1 <- m3[,"length.I"]
gbeta2 <- m3[,"length.O"]
galpha1 <- m3[,"(Intercept).I"]
galpha2 <- m3[,"(Intercept).O"]
gbeta31 <- m3[,"genderM.I"]
gbeta32 <- m3[,"genderM.O"]

gbeta1.mean <- mean(gbeta1)
gbeta2.mean <- mean(gbeta2)
galpha1.mean <- mean(galpha1)
galpha2.mean <- mean(galpha2)
gbeta31.mean <- mean(gbeta31)
gbeta32.mean <- mean(gbeta32)

#check values with frequentist approach
require(nnet)
multinom(choice ~ length + gender, data = data)
#values are similar/close

#Plot posterior densities
par(mar=c(4,4,1,1))
par(mfrow = c(4,1))
plot(density(beta1), xlab = expression(beta[1]), main = "Invertebrate Intercept")
abline(v=beta1.mean, col = "red")
plot(density(alpha1), xlab = expression(alpha[1]), main = "Invertebrate Coefficient")
abline(v=alpha1.mean, col = "red")
plot(density(beta2), xlab = expression(beta[2]), main = "Other Intercept")
abline(v=beta2.mean, col = "red")
plot(density(alpha2), xlab = expression(alpha[2]), main = "Other Coefficient")
abline(v=alpha2.mean, col = "red")


#for model with gender covariate
par(mar=c(4,4,1,1))
par(mfrow = c(3,2))
plot(density(gbeta1), xlab = expression(beta[1]), main = "Invertebrate Intercept")
abline(v=gbeta1.mean, col = "red")
plot(density(galpha1), xlab = expression(alpha[1]), main = "Invertebrate Coefficient")
abline(v=galpha1.mean, col = "red")
plot(density(gbeta2), xlab = expression(beta[2]), main = "Other Intercept")
abline(v=gbeta2.mean, col = "red")
plot(density(galpha2), xlab = expression(alpha[2]), main = "Other Coefficient")
abline(v=galpha2.mean, col = "red")
plot(density(gbeta31), xlab = expression(beta[31]), main = "Male, Invertebrate Coefficient")
abline(v=gbeta31.mean, col = "red")
plot(density(gbeta32), xlab = expression(beta[32]), main = "Male, Other Coefficient")
abline(v=gbeta32.mean, col = "red")

plot(density(m1[,"length.I"]), main = "Posterior Density of" ~ beta[1])
abline(v=beta1, col = 'red')

plot(density(m1[,"length.O"]))
abline(v=beta2, col = 'red')

#Response curves
response <- function(beta1, beta2, alpha1, alpha2){
  q <- matrix(NA, 500, 3)
  x <- seq(0, 5, length = 500)
  for (i in 1:500){
    #uncomment for whichever pi_i you want to evalutate
    #dist <- exp(alpha1 + beta1*x[i])/(1 + exp(alpha1 + beta1*x[i]) + exp(alpha2 + beta2*x[i])) #pi1
    #dist <- exp(alpha2 + beta2*x[i])/(1 + exp(alpha1 + beta1*x[i]) + exp(alpha2 + beta2*x[i])) #pi2
    dist <- 1/(1 + exp(alpha1 + beta1*x[i]) + exp(alpha2 + beta2*x[i])) #pi3
    q[i,] <- quantile(dist, probs = c(0.05, 0.5, 0.95))
  }
  return(q)
}

#two covariates
#0 for female 1 for male
x2 <- 1
q <- matrix(NA, 500, 3)
x <- seq(0, 5, length = 500)
for (i in 1:500){
  #uncomment for whichever pi_i you want to evalutate
  #dist <- exp(galpha1 + gbeta1*x[i] + gbeta31*x2)/(1 + exp(galpha1 + gbeta1*x[i] + gbeta31*x2) + exp(galpha2 + gbeta2*x[i] + gbeta32*x2)) #pi1
  #dist <- exp(galpha2 + gbeta2*x[i] + gbeta32*x2)/(1 + exp(galpha1 + gbeta1*x[i] + gbeta31*x2) + exp(galpha2 + gbeta2*x[i]) + gbeta32*x2) #pi2
  dist <- 1/(1 + exp(galpha1 + gbeta1*x[i] + gbeta31*x2) + exp(galpha2 + gbeta2*x[i] + gbeta32*x2)) #pi3
  q[i,] <- quantile(dist, probs = c(0.05, 0.5, 0.95))
}

 female.response1 <- q
 female.response2 <- q
 female.response3 <- q
 
 male.response1 <- q
 male.response2 <- q
 male.response3 <- q

# female.response1 <- response(fbeta1, fbeta2, falpha1, falpha2)
# female.response2 <- response(fbeta1, fbeta2, falpha1, falpha2)
# female.response3 <- response(fbeta1, fbeta2, falpha1, falpha2)
# male.response1 <- response(mbeta1, mbeta2, malpha1, malpha2)
# male.response3 <- response(mbeta1, mbeta2, malpha1, malpha2)


par(mar=c(4,4,1,1))
par(mfrow = c(3,1))
#distp1 <- q
plot(x, distp1[,1], type = "l", main = "Response Probability - Invertebrate", xlab = "length", ylab = expression(pi[1](x)), lty = 2, col = 'red', lwd = 2, ylim = c(0, 1))
lines(x, distp1[,2], lwd = 2)
lines(x, distp1[,3], lty = 2, col = 'red', lwd = 2)

#distp2 <- q
plot(x, distp2[,1], type = "l", main = "Response Probability - Other", xlab = "length", ylab = expression(pi[2](x)), lty = 2, col = 'red', lwd = 2, ylim = c(0, 1))
lines(x, distp2[,2], lwd = 2)
lines(x, distp2[,3], lty = 2, col = 'red', lwd = 2)

#distp3 <- q
plot(x, distp3[,1], type = "l", main = "Response Probability - Fish", xlab = "length", ylab = expression(pi[3](x)), lty = 2, col = 'red', lwd = 2, ylim = c(0, 1))
lines(x, distp3[,2], lwd = 2)
lines(x, distp3[,3], lty = 2, col = 'red', lwd = 2)

#seperating by gender
par(mar=c(4,4,1,1))
par(mfrow = c(3,2))

plot(x, female.response1[,1], type = "l", main = "Response Probability - Invertebrate", xlab = "length", ylab = expression(pi[1](x)), lty = 2, col = 'forestgreen', lwd = 2, ylim = c(0, 1))
lines(x, female.response1[,2], lwd = 2)
lines(x, female.response1[,3], lty = 2, col = 'forestgreen', lwd = 2)

plot(x, male.response1[,1], type = "l", main = "Response Probability - Invertebrate", xlab = "length", ylab = expression(pi[1](x)), lty = 2, col = 'darkmagenta', lwd = 2, ylim = c(0, 1))
lines(x, male.response1[,2], lwd = 2)
lines(x, male.response1[,3], lty = 2, col = 'darkmagenta', lwd = 2)


plot(x, female.response2[,1], type = "l", main = "Response Probability - Other", xlab = "length", ylab = expression(pi[2](x)), lty = 2, col = 'forestgreen', lwd = 2, ylim = c(0, 1))
lines(x, female.response2[,2], lwd = 2)
lines(x, female.response2[,3], lty = 2, col = 'forestgreen', lwd = 2)

plot(x, male.response2[,1], type = "l", main = "Response Probability - Other", xlab = "length", ylab = expression(pi[2](x)), lty = 2, col = 'darkmagenta', lwd = 2, ylim = c(0, 1))
lines(x, male.response2[,2], lwd = 2)
lines(x, male.response2[,3], lty = 2, col = 'darkmagenta', lwd = 2)


plot(x, female.response3[,1], type = "l", main = "Response Probability - Fish", xlab = "length", ylab = expression(pi[3](x)), lty = 2, col = 'forestgreen', lwd = 2, ylim = c(0, 1))
lines(x, female.response3[,2], lwd = 2)
lines(x, female.response3[,3], lty = 2, col = 'forestgreen', lwd = 2)

plot(x, male.response3[,1], type = "l", main = "Response Probability - Fish", xlab = "length", ylab = expression(pi[3](x)), lty = 2, col = 'darkmagenta', lwd = 2, ylim = c(0, 1))
lines(x, male.response3[,2], lwd = 2)
lines(x, male.response3[,3], lty = 2, col = 'darkmagenta', lwd = 2)

#------------------------------------END PROBLEM-------------------------------------

#Problem 3: developmental toxicity study involving ordinal categorical outcomes

mice <- data.frame(matrix(c(0, 15, 1, 281, 297,
                            62.5, 17, 0, 225, 242,
                            125, 22, 7, 283, 312,
                            250, 38, 59, 202, 299,
                            500, 144, 132, 9, 285), nrow = 5, byrow = T))
names(mice) <- c("concentration", "dead", 'malformation', 'normal','m')


#proportion of dead mice
mice$p1 <- mice$dead/mice$m
#number of alive
mice$n2 <- mice$m - mice$dead
#proportion of malformation
mice$p2 <- mice$malformation/mice$n2

model1 <- glm(p1 ~ concentration, weights=m, family=binomial(link=logit), data=mice)
model2 <- glm(p2 ~ concentration, weights=n2, family=binomial(link=logit), data=mice)

#estimates and standard errors
#alpha1 = -3.247934, 0.1576602
#beta1 = 0.006389, 0.0004348
#alpha2 = -5.70190, 0.332248
#beta2 = 0.01737, 0.001227


x <- seq(0, 500, len=100)

p1 <- predict(model1,newdata=list(concentration=x),type="response")
p2 <- predict(model2,newdata=list(concentration=x),type="response")
p3 <- 1 - p2

pi1 <- p1
pi2 <- p2 * (1-pi1)
pi3 <- 1 - pi1 - pi2

plot(mice$concentration, mice$p1, ylim = c(0, 1), col = 'darkmagenta', lwd = 2, xlab = "concentration (mg/kg per day)", ylab = "probability", main = "Response Curves")
lines(x, pi1, col = 'darkmagenta', lwd = 2)
points(mice$concentration, mice$p2*(1-mice$p1), col = 'forestgreen', lwd = 2)
lines(x, pi2, col = 'forestgreen', lwd = 2)
points(mice$concentration, 1-mice$p1-mice$p2*(1-mice$p1), col = 'goldenrod', lwd = 2)
lines(x, pi3, col = 'goldenrod', lwd = 2)
legend('topright', legend = c("dead", 'malformation', 'normal'), col = c('darkmagenta', 'forestgreen', 'goldenrod'), lwd = 2, cex = 0.7)

#-------------------------BAYESIAN GLM------------------------
x = mice$concentration
y = mice$dead
m = mice$m
N = length(y)
jinv <- 2*vcov(model1)

library(MASS)


prior <- function(x){
  return(dnorm(x[1], 0, 100, log = T) + dnorm(x[2], 0, 100, log = T))
}

k <- function(beta){
  return(exp(prior(beta) + sum(y*(beta[1] + beta[2]*x) - m*log(1 + exp(beta[1] + beta[2]*x)))))
}

#Metropolis-Hastings Algorithm for alpha 1 and beta 1
beta <- matrix(NA, ncol = 2, nrow = 10000)
#initialize at MLE
beta[1,] <- model1$coefficients
for (i in 2:10000){
  beta.squig <- mvrnorm(n = 1, beta[i-1,], jinv)
  q <- min(1, k(beta.squig)/k(beta[i-1,]))
  if(runif(1) < q){
    beta[i,] <- beta.squig
  }
  else{
    beta[i,] <- beta[i-1,]
  }
}



#Metropolis-Hastings Algorithm for alpha 2 and beta 2
x = mice$concentration
y = mice$malformation
m = mice$n2
N = length(y)
jinv <- 2*vcov(model2)

k <- function(beta){
  return(exp(prior(beta) + sum(y*(beta[1] + beta[2]*x) - m*log(1 + exp(beta[1] + beta[2]*x)))))
}

beta2 <- matrix(NA, ncol = 2, nrow = 10000)
#initialize at MLE
beta2[1,] <- model2$coefficients
for (i in 2:10000){
  beta.squig <- mvrnorm(n = 1, beta2[i-1,], jinv)
  q <- min(1, k(beta.squig)/k(beta2[i-1,]))
  if(runif(1) < q){
    beta2[i,] <- beta.squig
  }
  else{
    beta2[i,] <- beta2[i-1,]
  }
}

apply(beta[500:1000,], 2, mean)
apply(beta2[500:1000,], 2, mean)
#YAY

alpha1 <- beta[,1]
beta1 <- beta[,2]
alpha2 <- beta2[,1]
beta2 <- beta2[,2]

par(mar=c(4,4,1,1))
par(mfrow = c(4,2))
plot(ts(alpha1), ylab = expression(alpha[1]))
plot(density(alpha1), main = "Dead Intercept", xlab = expression(alpha[1]))
plot(ts(beta1), ylab = expression(beta[1]))
plot(density(beta1), main = "Dead Coefficient", xlab = expression(beta[1]))
plot(ts(alpha2), ylab = expression(alpha[2]))
plot(density(alpha2), main = 'Malformation Intercept', xlab = expression(alpha[2]))
plot(ts(beta2), ylab = expression(beta[2]))
plot(density(beta2), main = 'Malformation Coefficient', xlab = expression(beta[2]))


#point and interval estimates for the dose-response curve π(x)
q <- matrix(NA, 1000, 3)
x <- seq(0, 500, length = 1000)
for (i in 1:1000){
  #dist <- exp(alpha1 + beta1*x[i])/(1 + exp(alpha1 + beta1*x[i]))
  #dist <- (exp(alpha2 + beta2*x[i])/(1+exp(alpha2 + beta2*x[i])))*(1/(1+exp(alpha1+beta1*x[i])))
  dist <- 1 - (exp(alpha1 + beta1*x[i])/(1 + exp(alpha1 + beta1*x[i]))) - ((exp(alpha2 + beta2*x[i])/(1+exp(alpha2 + beta2*x[i])))*(1/(1+exp(alpha1+beta1*x[i]))))
  q[i,] <- quantile(dist, probs = c(0.05, 0.5, 0.95))
}

pi1 <- q
pi2 <- q
pi3 <- q

par(mar=c(4,4,1,1))
par(mfrow = c(3,1))
plot(x, pi1[,1], type = "l", lty = 2, col = 'magenta', ylab = expression(pi[1](x)), main = "Response Curve", ylim = c(0,1))
lines(x, pi1[,2])
lines(x, pi1[,3], lty = 2, col = 'magenta')

lines(x, pi2[,1], type = "l", lty = 2, col = 'magenta', ylab = expression(pi[2](x)), ylim = c(0,1))
lines(x, pi2[,2])
lines(x, pi2[,3], lty = 2, col = 'magenta')

lines(x, pi3[,1], type = "l", lty = 2, col = 'magenta', ylab = expression(pi[3](x)), ylim = c(0,1))
lines(x, pi3[,2])
lines(x, pi3[,3], lty = 2, col = 'magenta')