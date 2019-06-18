1.6
x=rpois(1000, .61)
hist(x, col="pink", main ="Histogram of Poisson λ=.61 n=1000")
hist(y, col="pink", main ="Histogram of Poisson λ=.61 n=10000")
y=rpois(10000, .61)

mean(x)
mean(y)
var(x)
var(y)

#compare with theoretical poisson
theo=dpois(0:5, lambda=.61)
Sample=table(x)/1000
cbind(theo, Sample)
mean(x)
var(x)

Sample=table(y)/10000
cbind(theo, Sample)

#1.8
#comparing the CDF with the empirical CDF
cdf=ppois(0:4, .61)
empcdf=cumsum(Sample)

#2.3
mtcars
?mtcars
boxplot(log(mtcars), main="Boxplot of mtcars")
pairs(mtcars, main = "mtcars data pairwise plot")\

#2.4
?mammals
library(MASS)
data()
mammals
r=mammals$brain/mammals$body
head(order(r))
tail(order(r))

#2.5
plot(r, mammals$body)

#2.6
LakeHuron
plot(LakeHuron, main="Timeplot of level of Lake Huron in feet")

plot(LakeHuron, ylab="Mean annual Level in feet")
abline(h=mean(LakeHuron))
lines(lowess(LakeHuron))
#transforming so the mean is stable over time
#do this by obtaining a time series of first differences x2-x1, x3-x2,...
d=diff(LakeHuron)
plot(d, ylab="First differences of mean annual level")
abline(h=0, lty=3)
lines(lowess(d))  
#the plot suggests that the mean of the differenced series is stable over time; notice that the lowess curve (solid line) is nearly horizontal and very close to the dotted horizontal line through 0.

#2.7
#create a matrix out of random numbers from uniform distribution
m=matrix(runif(9), nrow=3, ncol=3)
m
#compute the mean and variance of each row of the matrix
rowmeans=c(mean(m[1,]), mean(m[2,]), mean(m[3,]))
rowvar=c(var(m[1,]), var(m[2,]), var(m[3,]))

### can just use rowmeans(m)

hist(rowMeans(m), main = "Histogram of Sample Means", col = "plum1", xlab="Sample Means")

#The off-diagonal elements in the variance-covariance matrix are the sample covariance, and theoretically the covariances should be zero: the numbers in columns should be uncorrelated if in fact the table represents independent and identically distributed (iid) numbers. The sample correlations are each close to zero in absolute value:

#sample covariance
diag(var(m))
cor(m)

#the runif data have moderately strong linear correlations (close to -1 or 1)


#We are interested in the distribution of the row means
#larger matrix to get 500 row means

hist(rowMeans(m), main="Histogram of sample means from unif", xlab = "Sample Means")

hist(rowMeans(m), prob=TRUE, main="Probability Histogram of Sample Means", xlab = "Sample Means", col="orchid1")
plot(density(rowMeans(m)), main="Density Plot of the Sample Means", col = "orchid1")

#another way to compare CURVE DOESNT WORK
truehist(rowMeans(m), main = "True Histogram of Sample Means", xlab = "Sample Means")
curve(dnorm(x, 1/2, sd=sqrt(1/36)), add=TRUE)

#another way to check normality
qqnorm(rowMeans(m))
qqline(rowMeans(m), col = "orchid1")


#2.8
#now make a 400 x 10 matrix for a sample size of n=10
largern=rowMeans(matrix(runif(4000), nrow=10, ncol=400))
hist(largern, main="Histogram of Sample Means with n=10", xlab="Sample Means", col = "seagreen1")

hist(largern, prob=TRUE, main="Probability Histogram of Sample Means from Uniform Dist", xlab = "Sample Means", col = "seagreen1")
plot(density(largern), main = "Density Estimate for the Sample Means", col="seagreen1")

truehist(largern, main = "True Histogram of Sample Means", xlab = "Sample Means")
curve(dnorm(x, 1/2, sd=sqrt(1/420)), add=TRUE)

qqnorm(largern)
qqline(largern, col="seagreen1")

#2.12
mam=order(mammals$brain)
head(mam)
tail(mam)
mammals$brain[15]

#how to find which animal it is
cbind(1:62, mammals)
mammals[30,]

#2.13

plot(mammals$body, mammals$brain, main = "Mammals Body vs. Brain (Orininal Scale)", xlab="Body", ylab="Brain")
want=mammals[c("Cat", "Cow", "Human"),]
polygon(want)
text(want, rownames(want))

#3.2
#simulate 1000 rolls of a fair die
die1=sample(6, 1000, replace = TRUE)
die2=sample(6, 1000, replace = TRUE)
die.sum = die1 + die2
table(die.sum)/1000

actualprob=c(1/36, 2/36, 3/36, 4/36, 5/36, 6/36, 5/36, 4/36, 3/36, 2/36, 1/36)
plot(table(die.sum)/1000, main = "Sum of Two Die", xlab = "Sum", ylab="Sample Proportion", col = "orchid1")
points(actualprob, col = "seagreen1")
legend("topright",c("Simulated","Actual prob"), lty = c(1,1), col = c("orchid1", "seagreen1"))


#3.3
#part A
expected=c(dbinom(0, 4, .312)*70, dbinom(1, 4, .312)*70, dbinom(2, 4, .312)*70, dbinom(3, 4, .312)*70+dbinom(4, 4, .312)*70)
expected

observed=c(17, 31, 17, 5)
chisq.test(expected, observed)

#partB
zero=dbinom(0, 5, .312)*25
one=dbinom(1, 5, .312)*25
two=dbinom(2, 5, .312)*25
threeormore=dbinom(3, 5, .312)*25+dbinom(4, 5, .312)*25+dbinom(5, 5, .312)
expected=c(one, two, zero+threeormore)#combining first and last bin to make expected counts>5
observed=c(5, 4, 16)
chisq.test(expected, observed)


#3.4
twins.dat_$AGE
#cut function categorizes
a=table(cut(twins.dat_$AGE, c(0, 30, 40, 50, 100)))/183
barplot(a, main = "Age of Twin 1", xlab = "Age group", ylab = "Proportion", col = 101)
age=subset(twins.dat_, select = AGE, subset = !is.na(AGE), drop = T)

#3.5
#PART A
twins$AGE
Age=cut(twins.dat_$AGE, c(0, 30, 40, 50, 100))
HourlyWage=cut(twins$HRWAGEL, c(0, 7, 13, 20, 150))
length(HourlyWage)

#PARTB
tab=table(Age, HourlyWage)

#PARTC
p=prop.table(tab)

#PARTD
mosaicplot(p, main = "Hourly Wage in Relation to Age", color = "seagreen1" )
plot(p)
barplot(p)

#3.6
#PARTA
test=chisq.test(tab)

#PARTB
plot(as.vector(test$residuals), main = "Plot of the Residuals", ylab = "Residuals")
abline(h=0)
identify(as.vector(test$residuals), n=2, labels = "+", col="hotpink")
legend("topleft", legend="+ |residual|>2", text.col = "hotpink")


#PARTC
mosaicplot(tab, shade=T, main = "Plot of the Residuals")

#3.7
#PARTA
die1=sample(6, 1000, replace = T)
die2=sample(6, 1000, replace = T)
#PARTB
max.rolls=pmax(die1, die2)
sum.rolls=die1+die2
#PARTC
j=table(max.rolls, sum.rolls)
#PARTD
t=chisq.test(j)
t$residuals
plot(as.vector(t$residuals), main = "Plot of the Residuals", ylab = "Residuals")
abline(h=0)
abline(h=2, col="hotpink")
abline(h=-2, col="hotpink")


#3.8
pidigits =
  read.table("http://www.itl.nist.gov/div898/strd/univ/data/PiDigits.dat",
             skip=60)
#PARTA
p=table(pidigits, exclude = 0)
#PARTB
barplot(p, main = "First 5000 Didgits of Pi", col = "hotpink")
#PARTC
t=chisq.test(p)

t$residuals
plot(as.vector(t$residuals), main = "Plot of the Residuals", ylab = "Residuals")
abline(h=0)
abline(h=2, col="hotpink")
abline(h=-2, col="hotpink")
