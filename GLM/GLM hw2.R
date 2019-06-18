#GLM Homework 2

#Problem 3

length <- c(551, 651, 832, 375, 715, 868, 271, 630, 491, 372, 645, 441, 895, 458, 642, 492, 543, 842, 905, 542, 522, 122, 657, 170, 738, 371, 735, 749, 495, 716, 952, 417)
faults <- c(6, 4, 17, 9, 14, 8, 5, 7, 7, 7, 6, 8, 28, 4, 10, 4, 8, 9, 23, 9, 6, 1, 9, 4, 9, 14, 17, 10, 7, 3, 9, 2)

#Part A
#Use R to fit a Poisson GLM, with logarithmic link
model <- glm(faults ~ length, family = poisson(link = 'log'))
summary(model)

#Coefficients:
#               Estimate    Std.  Error z value Pr(>|z|)    
# (Intercept)   0.9717506   0.2124693   4.574   4.79e-06 ***
# length        0.0019297   0.0003063   6.300   2.97e-10 ***

#For every unit increase in fabric length, the number of faults increases 0.00193 on average. The disperson parameter for poisson family taken to be 1.


#Part B
#Use the quasipoisson “family” in R to fit the regression model for the response means in (1) using the quasi-likelihood estimation method, which allows for a dispersion parameter in the response variance function. Discuss the results.

model2 <- glm(faults ~ length, family = quasipoisson(link = 'log'))
summary(model2)

# Coefficients:
#             Estimate    Std. Error  t value  Pr(>|t|)    
# (Intercept) 0.9717506  0.3095033    3.140    0.003781 ** 
# length      0.0019297  0.0004462    4.325    0.000155 ***

#Using the quasi-likelihood estimation method, we obtain the same MLE estimates. The standard error's are slightly higher for this method and the dispersion parameter is 2.122.

#Part C
#Derive point estimates and asymptotic interval estimates for the linear predictor, η0 = β1+ β2x0, at a new value x0 for length of roll, under the standard (likelihood) estimation method from part (a), and the quasi-likelihood estimation method from part (b). Evaluate the point and interval estimates at x0 = 500 and x0 = 995.

#x0 = 500 Poisson
#point estimate for eta0
eta <- model$coefficients[1] + model$coefficients[2]*500
#1.936624 

#interval estimate for eta0
#vcov = phi^tilde j^-1(beta)
print(c(eta - 1.96* sqrt(c(1,500) %*% vcov(model) %*% c(1, 500)), eta + 1.96* sqrt(c(1,500) %*% vcov(model) %*% c(1, 500))))
#(1.783435, 2.089812)

#x0 = 500 Quasi-Poisson
#point estimate
eta <- model2$coefficients[1] + model2$coefficients[2]*500
#1.936624

#interval estimate
print(c(eta - 1.96* sqrt(c(1,500) %*% vcov(model2) %*% c(1, 500)), eta + 1.96* sqrt(c(1,500) %*% vcov(model2) %*% c(1, 500))))
#(1.713475, 2.159773)

#x0 = 995 Poisson
#point estimate
eta <- model$coefficients[1] + model$coefficients[2]*995
#2.891849

#interval estimate
print(c(eta - 1.96* sqrt(c(1, 995) %*% vcov(model) %*% c(1, 995)), eta + 1.96* sqrt(c(1,995) %*% vcov(model) %*% c(1, 995))))
#(2.662676, 3.121021)

#x0 = 995 Quasi-Poisson
#point estimate
model2$coefficients[1] + model2$coefficients[2]*995
#2.891849

#interval estimate
print(c(eta - 1.96* sqrt(c(1, 995) %*% vcov(model2) %*% c(1, 995)), eta + 1.96* sqrt(c(1,995) %*% vcov(model2) %*% c(1, 995))))
#(2.558014, 3.225683)


#Problem 4
#Build a Poisson GLM to study the effect of the covariates (jet fuel concentration and organism strain) on the number of Ceriodaphnia organisms.
names(ceriodaphnia) <- c('n.organisms', 'jetfuel.concen', 'strain')
plot(ceriodaphnia)

ceriodaphnia$strain = as.factor(ceriodaphnia$strain)
color_easy = c("mediumturquoise", "mediumvioletred")[ceriodaphnia$strain]

#origional plot
plot(ceriodaphnia$jetfuel.concen, ceriodaphnia$n.organisms, main = 'Jet Fuel Concentration vs. Number of Organisms', col = color_easy, xlab = "Jet Fuel Concentration", ylab = "Number of Organisms", lwd = 2)
legend('topright', legend = c("Strain 1", "Strain 0"), col = color_easy, lwd = 2)

#----------------------------------------------------------------------------

#log response
plot(ceriodaphnia$jetfuel.concen, log(ceriodaphnia$n.organisms), col = color_easy, xlab = "Jet Fuel Concentration", ylab = "log(number of organisms)", lwd = 2)
legend('topright', legend = c("Strain 1", "Strain 0"), col = color_easy, lwd = 2)
#log(response) produces a linear association between the two variables so we use a log link function.

#Model1: mu = exp(b0 + b1*xi1 + b2*xi2)
model <- glm(n.organisms ~ jetfuel.concen + strain, data = ceriodaphnia, family = poisson(link = 'log'))
summary(model)
BIC(model)
sum(residuals(model, "pearson")^2)
#AIC: 415.95
#BIC: 422.7
#Residual Deviance: 86.38
#average pearson residual: 1.14
#sum of squared pearson residuals: 79.83

#strain = 0
abline(4.18, - 1.543, col = "mediumturquoise")
#strain = 1
abline(4.18 + 0.28, -1.543, col = "mediumvioletred")

ceriodaphnia$strain <- as.numeric(as.character(ceriodaphnia$strain))

#----------------------------------------------------------------------------

plot(ceriodaphnia$jetfuel.concen, sqrt(ceriodaphnia$n.organisms), col = color_easy, xlab = "Jet Fuel Concentration", ylab = "sqrt(number of organisms)", lwd = 2)
#legend('topright', legend = c("Strain 1", "Strain 0"), col = color_easy, lwd = 2)

model2 <-  glm(n.organisms ~ jetfuel.concen + strain, data = ceriodaphnia, family = poisson(link = 'sqrt'))

#strain = 0
abline(7.7267, -3.59, col = "mediumturquoise")
#strain = 1
abline(7.7276 + 0.5868, -3.59, col = "mediumvioletred")

summary(model2)
BIC(model2)
sum(residuals(model2, "pearson")^2)
mean(residuals(model2, "pearson")^2)

#AIC: 473.75
#BIC: 480.5
#Residual deviance: 144.18
#X^2: 149.94
#average pearson residuals: 2.142 

#----------------------------------------------------------------------------

plot(sqrt(ceriodaphnia$jetfuel.concen), log(ceriodaphnia$n.organisms), col = color_easy, xlab = "sqrt(Jet Fuel Concentration)", ylab = "log(number of organisms)", lwd = 2)
#legend('topright', legend = c("Strain 1", "Strain 0"), col = color_easy, lwd = 2)

model3 <-  glm(n.organisms ~ sqrt(jetfuel.concen) + strain, data = ceriodaphnia, family = poisson(link = 'log'))

#strain = 0
abline(4.253, -1.658, col = "mediumturquoise")
#strain = 1
abline(4.253 + 0.275, -3.59, col = "mediumvioletred")

BIC(model3)
sum(residuals(model3, "pearson")^2)
mean(residuals(model3, "pearson")^2)

#AIC: 493.9
#BIC: 500.6
#Residual deviance: 164.3
#X^2: 151.72
#Average pearson residual: 2.167

#----------------------------------------------------------------------------

plot(log(ceriodaphnia$jetfuel.concen+1), log(ceriodaphnia$n.organisms), col = color_easy, xlab = "log(Jet Fuel Concentration + 1)", ylab = "log(number of organisms)", lwd = 2)


model4 <-  glm(n.organisms ~ log(jetfuel.concen+1) + strain, data = ceriodaphnia, family = poisson(link = 'log'))

BIC(model4)
sum(residuals(model4, "pearson")^2)
mean(residuals(model4, "pearson")^2)

#strain = 0
abline(4.232, -2.38, col = "mediumturquoise")
#strain = 1
abline(4.232 + 0.275, -2.38, col = "mediumvioletred")

#AIC: 437.2
#BIC: 443.95
#Residual deviance: 164.3
#X^2: 96.96
#Average pearson residual: 1.385


#----------------------------------------------------------------------------

#only quantitative variables
model5 <- glm(n.organisms ~ jetfuel.concen, data = ceriodaphnia, family = poisson(link = 'log'))

plot(ceriodaphnia$jetfuel.concen, log(ceriodaphnia$n.organisms), col = color_easy, xlab = "Jet Fuel Concentration", ylab = "log(number of organisms)", lwd = 2)
abline(4.327, -1.543)

BIC(model5)
sum(residuals(model5, "pearson")^2)
mean(residuals(model5, "pearson")^2)

#AIC: 446.6
#BIC: 451.07
#Residual Deviance: 119
#X^2: 115.135
#Average pearson residual: 1.645

#----------------------------------------------------------------------------

model6 <- glm(n.organisms ~ jetfuel.concen^2 + strain, data = ceriodaphnia, family = poisson(link = 'log'))

plot(ceriodaphnia$jetfuel.concen^2, log(ceriodaphnia$n.organisms), col = color_easy, xlab = "Jet Fuel Concentration^2", ylab = "log(number of organisms)", lwd = 2)

#strain = 0
abline(4.18, -1.543, col = "mediumturquoise")
#strain = 1
abline(4.18 + 0.275, -1.543, col = "mediumvioletred")

BIC(model6)
sum(residuals(model6, "pearson")^2)
mean(residuals(model6, "pearson")^2)

#AIC: 416
#BIC: 422.7
#Residual Deviance: 86.38
#X^2: 79.83
#Average pearson residual: 1.14



m <- matrix(c(1,2,3,4,5,6,7,7,7),nrow = 3,ncol = 3,byrow = TRUE)

layout(mat = m,heights = c(0.4,0.4,0.2))

plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend('top', inset = 0, legend = c("Strain 1", "Strain 0"), col = color_easy,horiz = T, lwd = 2)
par(mar = c(2,2,1,1))