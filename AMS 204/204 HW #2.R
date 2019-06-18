#4.1
plot(cars$speed, cars$dist, xlab = "Speed (mph)", ylab = "Stopping Distance (ft)", main="Stopping Distance vs. Speed", col="red", pch=17)

#4.4
par(mfrow=c(2, 2))
attach(mtcars)
plot(disp, mpg, xlab = "Displacement")
plot(wt, mpg, xlab = "Weight")
plot(drat, mpg, xlab = "Rear Axle Ratio" )
plot(hp, mpg, xlab = "Horsepower")

#4.5
house=function(x, y, ...){
  lines(c(x - 1, x + 1, x + 1, x - 1, x - 1),
        c(y - 1, y - 1, y + 1, y + 1, y - 1), ...)
  lines(c(x - 1, x, x + 1), c(y + 1, y + 2, y + 1), ...)
  lines(c(x - 0.3, x + 0.3, x + 0.3, x - 0.3, x - 0.3),
        c(y - 1, y - 1, y + 0.4, y + 0.4, y - 1), ...)
}
plot.new()
plot.window(xlim = c(0, 10), ylim = c(0, 10))
lines(house(1, 1), house(4, 2), house(7,6))
house(0, 7, col="purple", lty = 17)
box()

#4.6
curve(dbeta(x, 2, 6), from=0, to=1)
curve(dbeta(x, 4, 4), from=0, to=1, add = T)
curve(dbeta(x, 6, 2), from=0, to=1, add = T)
title(expression(f(y)==frac(1,B(a,b))*y^{a-1}*(1-y)^{b-1}))
text(.05, 2.7, labels = "Beta(2, 6)")
text(.4, 1.5, labels = "This is how to make words")

#4.7
faithful$leng = ifelse(faithful$eruptions < 3.2, "short", "long")
faithful
library(lattice)
bwplot(waiting ~ leng, data = faithful, main = "Waiting Time by Length")
densityplot(~waiting, group = leng, data = faithful, auto.key=list(space="top"), xlab = "Waiting Time", main = "Faithful Waiting Time Density Plot")

#5.1
#PART A
col=subset(college, complete.cases(college))
stripchart(Pct.20 ~ Tier, method = "stack", main="Percentage of Small Classes in National Universities", xlab="Small Class Percentage", data = col)
#PART B
identify(col$Pct.20, col$Tier, n=1, labels=col$School)
#cant get identify to work
#PARTC
median(col$Pct.20)
abline(v=45, col="hotpink")

#5.2
#PART A
big=college$Pct.50[!is.na(college$Pct.50)]
plot(college$Pct.20, college$Pct.50, main = "Percentage of Class Sizes", xlab = "Small Class Percentage", ylab = "Large Class Percentage")
#PART B
fit=line(col$Pct.20, col$Pct.50)
fit
abline(coef(fit))
#PART C
23.1-.2667*60
#PART D
plot(col$Pct.20, fit$residuals, xlab = "Small Class Percentage", ylab = "Residuals", main = "Residuals vs. Pct.20")
abline(h=0)
abline(h=10, col="hotpink")
abline(h=-10, col="hotpink")
identify(col$Pct.20, fit$residuals, n=7, labels = col$School)

plot(col$Pct.50, fit$residuals, xlab = "Large Class Percentage", ylab = "Residuals", main = "Residuals vs. Pct.50")
abline(h=0)
abline(h=10, col="hotpink")
abline(h=-10, col="hotpink")
identify(col$Pct.50, fit$residuals, n=7, labels = col$School)

#5.5
#PARTA
hist(col$Full.time, xlab = "Full Time Percentages")
froot=sqrt(col$Full.time) - sqrt(100 - col$Full.time)
flog = log(col$Full.time + 0.5) - log(100 - col$Full.time + 0.5)
hist(froot)
hist(flog)
truehist(flog)
curve(dnorm(x, m, stan), add = T)
m=mean(flog)
stan=sd(flog)
m+stan
m-stan


#5.7
stripchart(Alumni.giving ~ Tier, method = "stack", main="Alumni Giving Rates by Tier", xlab="Alumni Giving Rate", data = col, ylab="Tier")
identify(col$Alumni.giving, col$Tier, n=3, labels = col$School)

#PART D
stripchart(sqrt(Alumni.giving) ~ Tier, method = "stack", main="Alumni Giving Rates by Tier (square root)", xlab="Alumni Giving Rate", data = col, ylab="Tier")

boxplot(log(Alumni.giving) ~ Tier, method = "stack", main="Alumni Giving Rates by Tier (log)", xlab="log(Alumni Giving Rate)", data = col, ylab="Tier", horizontal =T)



#Extra questions
GSS=subset(`GSS2015_DST_09.(1)`, complete.cases(`GSS2015_DST_09.(1)`))
attach(GSS)

#Gets rid of old data
GSS$X2014olda=NULL

sci=subset(GSS, Field=="Science", select = c(X2010, X2011, X2012, X2013, X2014newa, X2015))
eng=subset(GSS, Field=="Engineering", select = c(X2010, X2011, X2012, X2013, X2014newa, X2015))
heal=subset(GSS, Field=="Health", select = c(X2010, X2011, X2012, X2013, X2014newa, X2015))
stat=subset(GSS, Field=="Statistics", select = c(X2010, X2011, X2012, X2013, X2014newa, X2015))

s=t(sci)
e=t(eng)
h=t(heal)
st=t(stat)
new=cbind(s, e, h)


barplot(new, beside = T, xlab = "Field", ylab = "Number of Graduate Students", names.arg = c("Science", "Engineering", "Health"), ylim = c(0,450000), col = cm.colors(6), main = "Number of Graduate Students Over the Years")
options(scipen = 7)
names=c("2010", "2011", "2012", "2013", "2014", "2015")
legend("topright", legend = names, pch = 15, col = cm.colors(6), bty = "n")

barplot(st, beside = T, col = cm.colors((6)), ylim = c(0,8000), main = "Graduate Students in Statistics", names.arg = "Statistics", ylab = "Number of Graduate Students")
legend("topleft", legend = names, pch = 15, col = cm.colors(6), bty = "n", horiz = T, cex=.8)

svst=cbind(s, st)
chisq.test(svst)

evst=cbind(e, st)
chisq.test(evst)
