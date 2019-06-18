#8.1
#here's the data
y1 = c(22, 26)
y2 = c(28, 24, 29)
y3 = c(29, 32, 28)
y4 = c(23, 24)
y = c(y1, y2, y3, y4)
Model = c(rep("A", 2), rep("B", 3), rep("C", 3), rep("D", 2))
mileages=data.frame(y, Model)

#visualize
boxplot(y ~ Model, xlab = "Model")
stripchart(y ~ Model, vertical=T, xlab = "Model")

#comparing means of the 4 groups
meansd = function(x) c(mean=mean(x), sd=sd(x))
by(y, Model, FUN=meansd)
#SD are close for the groups

#are the variances equal?
var.test(mileages$y, mileages$Model)

#test, don't assume equal var
oneway.test(y ~ Model)
#Fail to reject. There is no evidence that the means are different.

##########################################################################

#8.3
attach(iris)
boxplot(Sepal.Length ~ Species, ylab = "Sepal Length", xlab = "Species")
stripchart(Sepal.Length ~ Species, vertical = T)

#similar variances, different centers
# Let us compare the means and standard deviations of the three groups.
by(Sepal.Length, Species, FUN=meansd)

#check variances
var.test(Sepal.Length, Species)
#use test with equal variances

#test
oneway.test(Sepal.Length ~ Species, var.equal=TRUE)
#small p-value=at least one of the means differs

#Effects model
#general
#y_ij = mu + tau_j + epsilon_ij
#y_ij denotes the ith flower in the jth species
#Set tau_1=0 restriction j= 1, 2, 3
#mu=5.006 the sample mean for the setosa group
#tau_2=0.93 the difference between the mean of versicolor and the mean of setosa
#tau_3=1.582 the difference between the mean of virginica and the mean of setosa
#epsilon_ij is our random error variable
#with this restriction model for the three species, we have four parameters, mu, tau_2, tau_3, and sigma^2
L = lm(Sepal.Length ~ Species)
L

#The intercept is equal to the sample mean for setosa, but the other two coefficients correspond to the differences between the sample means
5.006-6.588
predict(L)

#ANOVA table
anova(L)$sigma^2
summary(L)$sigma^2
#we estimate the common variance (sigma^2) with sigma^2 hat, the residual MSE, which is 0.265 with 147 degrees of freedom


M = aov(Sepal.Length ~ Species)
model.tables(M, type="means")
#the differences between the group means and the grand mean
model.tables(M)

##########################################################################
#8.4
#Checking ANOVA model assumptions:
#Independence of cases – independent because different types of flowers
#Normality – the distributions of the residuals are normal. yes, below
#errors (epsilon_ij) need to be iid with mean 0 and constant variance sigma^2


#plot residuals vs fits
plot(L$fit, L$res)
abline(h=0) #add horizontal line through 0

#Normal-QQ plot of residuals with reference line
qqnorm(L$res)
qqline(L$res)
 
#The residuals are approximately symmetric about zero and have approximately equal variance. In the Normal-QQ plot the residuals lie approximately along the reference line. The plots do not reveal any severe departure from the assumptions for this model.

#########################################################################

#8.5
attach(PATIENT)
PATIENT = na.omit(PATIENT)

boxplot(PATIENT)
stripchart(PATIENT, vertical=T)

L=lm(PATIENT)
L

#plot residuals vs fits
plot(L$fit, L$res)
abline(h=0) #add horizontal line through 0

qqnorm(L$res)
qqline(L$res)


#The residuals are not symmetric about zero and do not have approximately equal variance. In the Normal-QQ plot not all of the residuals lie along the reference line. The plots reveal severe departure from the assumptions for this model.

#we will use a log transformation
boxplot(log(PATIENT))
stripchart(log(PATIENT), vertical=T)
L=lm(log(PATIENT))
L

#plot residuals vs fits
plot(L$fit, L$res)
abline(h=0) #add horizontal line through 0

qqnorm(L$res)
qqline(L$res)

plot(L, which = 1, add.smooth = F)

#########################################################################
#9.1)
attach(personal.bgsu.edu)
str(personal.bgsu.edu)

player = factor(personal.bgsu.edu$block)
time = data.frame(times, method, player)
str(time)

#Analyze the data assuming a randomized block design with time to round first base as the response, method of rounding first as the treatment, and player as the block variable. 
L = aov(times ~ method + player)
summary(L)
#Here the F statistics for both treatment (exam) and block (student) are very large (the p-values are very small) and both tests are significant.
#We can conclude that there are significant differences among the mean time to round first base by method and significant differences between mean player times.

CIs = TukeyHSD(L, which=1)
CIs

#From the table displayed by the TukeyHSD function we can identify significant differences between wide angle and narrow angle as well as wide angle and round out, but note that there is no significant difference between round out and narrow angle.

#check with a plot
plot(CIs, las=1)
#We can see from the plot that the interval for round out and narrow anlge inclues zero and therefore are not significantly different

#plot residuals
plot(L, which=1:2)



#The Normal QQ plot does not reveal a serious departure from the normality assumptions, although there are some outliers labeled on the QQ plot. 
#The plot of residuals vs fits in has a ‘M' shape, which may indicate nonconstant error variance.

boxplot(times ~ player, xlab="Student Number", ylab="Score", data = time)
