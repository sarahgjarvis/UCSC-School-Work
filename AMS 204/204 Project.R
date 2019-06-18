str(adult.test)
attach(adult.test)

str(project)
attach(project)

#Model 1 
M1=glm(Bincome ~ Age + Work.Class + Education + edunum + Marrital.Status + Occupation + Relationship + Race + Gender + capital_gain + capital_loss + Hours.per.week + Native.Country, family	=	binomial("logit"))
s=summary(M1)
c=anova(M1,test="Chisq")
#AIC=10381

#EDA
boxplot(log(Age) ~ Income, ylab = "Age", main = "Boxplot of Age by Income")
plot(table(Gender, Income), shade = T, main = "Mosaic Plot Income by Gender")
plot(table(Race, Income), las=2, shade = T, main = "Mosaic Plot Income by Race")
plot(table(Education, Income), las = 2, shade = T, main = "Mosaic plot Income by Education")
plot(table(Work.Class, Income), las = 2, shade = T, main = "Mosaic plot Income by Work Class")
plot(table(Marrital.Status, Income), las = 2, shade = T, main = "Mosaic plot Income by Marrital Status")
plot(table(Occupation, Income), las = 2, shade = T, main = "Mosaic plot Income by Occupation")
plot(table(Relationship, Income), las = 2, shade = T, main = "Mosaic plot Income by Relationship")
boxplot(Hours.per.week ~ Income)
plot(Income, edunum, col="violet", main = "Boxplot of Education Number", ylab="edunum")


#printing tables
install.packages("xtable")
library(xtable)
newobject<-xtable(s)
print.xtable(newobject, type="html", file="filename.html")



#clustering by GDP
attach(ProjectGDP)
distGDP=dist(ProjectGDP)
h=hclust(distGDP, method="complete") 
plot(h, labels = Country, main = "Clustering Native Country by GDP")
legend("topleft", leg.txt, pch = 15, col = c("orange", "yellow", "purple", "dodgerblue"), horiz = T, cex=.8, bty = "n")
leg.txt=c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4")

#Changing cluster to a factor
attach(adult.test)
clust= factor(Cluster)
project = data.frame(Age, Work.Class, Education, edunum, Marrital.Status, Occupation, Relationship, Race, Gender, capital_gain, capital_loss, Hours.per.week, clust)
str(project)

#Model 2 clustering Native Country
M2=glm(Bincome ~ Age + Work.Class + Education + edunum + Marrital.Status + Occupation + Relationship + Race + Gender + capital_gain + capital_loss + Hours.per.week + clust, family	=	binomial("logit"))
s=summary(M2)
c=anova(M2,test="Chisq")
#AIC=10148

#Model 3 because cluster variable was not significant
M3=glm(Bincome ~ Age + Work.Class + Education + edunum + Marrital.Status + Occupation + Relationship + Race + Gender + capital_gain + capital_loss + Hours.per.week, family	=	binomial("logit"))
s=summary(M3)
c=anova(M3,test="Chisq")

#M3 has the lowest AIC of 10145


