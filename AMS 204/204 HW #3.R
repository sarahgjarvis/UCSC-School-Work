#6.2
attach(nyc.marathon)
women.marathon=subset(nyc.marathon, Gender=="female", select = c(Minutes, Age))
t.test(women.marathon$Age, mu=36.1)
hist(women.marathon$Age)

wilcox.test(women.marathon$Age, mu=36.1, conf.level = .9)

t.test(women.marathon$Age, mu=36.1, conf.level = .90, conf.int=T)


#6.3
men.marathon=subset(nyc.marathon, Gender=="male", select = c(Minutes, Age))
hist(men.marathon$Age)
t.test(men.marathon$Age, women.marathon$Age, alternative = "greater", conf.int=T)

length(women.marathon$Age)

boxplot(women.marathon$Age, men.marathon$Age)
t.test(men.marathon$Age, women.marathon$Age, conf.int=T, var.equal = T)
#dont know which pvalue to use
var.test(women.marathon$Age, men.marathon$Age)

#6.5
#part a
diff=buffalo.cleveland.snowfall$Buffalo-buffalo.cleveland.snowfall$Cleveland
t.test(diff)
hist(diff, col = "purple", main = "Histogram of Differences", xlab = "Differences")
qqnorm(diff, col = "purple")
qqline(diff)
max(diff)

#REMOVES OUTLIERS
remove_outliers <- function(PATIENT, na.rm = TRUE, ...) {
  qnt <- quantile(PATIENT, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(PATIENT, na.rm = na.rm)
  y <- diff
  y[diff < (qnt[1] - H)] <- NA
  y[diff > (qnt[2] + H)] <- NA
  y
}
y <- remove_outliers(PATIENT)

hist(y, col = "purple", main = "Histogram of Differences", xlab = "Differences")
qqnorm(y, col="purple")
qqline(y)
t.test(y)


#6.7
qqnorm(prezheight$Height.cm., main = "Normal Q-Q Plot of Winner Height")
qqline(prezheight$Height.cm.)

qqnorm(prezheight$Oheight.cm., main = "Normal Q-Q Plot of Loser Height")
qqline(prezheight$Oheight.cm.)

hist(prezheight$Height.cm., main = "Winning President Height", xlab = "Height")
hist(prezheight$Oheight.cm., main = "Losing President Height",xlab = "Height")

t.test(prezheight$Height.cm., prezheight$Oheight.cm., paired = T)


library(MASS)
attach(mammals)
plot(mammals$body, mammals$brain)

#7.1
#fit a regression line
line=lm(body ~ brain)
abline(line$coef)
lm(body ~brain)
#plot residuals
plot(line, which=1, add.smooth=FALSE)

#7.2
plot(log(brain), log(body))
line=lm(log(brain) ~ log(body)) 

#estimate of the error variance
(summary(line)$sigma)^2

abline(line$coefficients)
#plot residuals
plot(line, which=1, add.smooth=FALSE)
qqnorm(line$residuals)
summary(line)

(cor(log(brain), log(body)))^2

#7.7
L = lm(cars$dist ~ cars$speed)
summary(L)
#plot residuals
plot(L, which=1, add.smooth=FALSE)

#7.8
#quadratic model
speed2=(cars$speed)^2
quadratic.model= lm(cars$dist ~ cars$speed + speed2)


quadratic.model
plot(cars$speed,cars$dist, col="blue",
     xlab="speed",ylab="distance",
     main="Distance vs Speed")
curve(2.47014 + 0.91329*x +0.09996*x^2, add=TRUE)