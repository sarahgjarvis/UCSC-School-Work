install.packages("ISwR")
library(ISwR)
tb.dilute


#1
#Perform a two-way ANOVA on tb.dilute
#The two-way ANOVA compares the mean differences between groups that have been split on two independent variables (called factors).
#Independent variables: animal, logdose
#Dependent variable: reaction

#Interaction plots
with(data = tb.dilute, expr = {
  interaction.plot(logdose, animal, response = reaction, lwd=2)
  interaction.plot(animal, logdose, response = reaction)
  })
#Interaction occurs between animals 5 and 6, 3 and 4, 3 and 5, 4 and 5. 
#Interaction also occurs between logdose 0.5 and 0.
#We will use an ANOVA test to evaluate the statistical significance of the effects.

#Two-way ANOVA
L=aov(reaction ~ animal + logdose, data = tb.dilute)
anova(L)

#The animal and logdose effects are significant.

model.tables(L, type = "means")
#The means looks slightly different for each group, we will do a tukey multiple comparison to see if the means are significantly different.

model_L=lm(reaction~animal+logdose, data = tb.dilute)
summary(model_L)
#All variables of our model are significant except animal 2.

TukeyHSD(L, which = c("animal", "logdose"))

CIs = TukeyHSD(L, which=1)
plot(CIs, las=1)
plot(TukeyHSD(L, which = c("animal", "logdose"))
)

#plot residuals vs fits
plot(L, which=1:2)
#Something about the residual plot
#The qq-plot reveals 3 outliers.

#Does this deviate from the assumptions?

#H01:The mean reaction is the same for all 6 animals
#H02:The mean reaction is the same for all dosages
#H03:There is no interaction between factors

#With a f-stat of 8.2641 and a p-value of .002527, we conclude that at least one mean of animal reaction differs.
#With a f-stat of 36.199 and a p-value of 1.295e-05, we conclude that at least one mean dosage reaction differs.

#we don't have enough degrees of freedom to do hypothesis testing on interaction effects

#Based on the Tukey results, there are significant differences between the logdose groups; -0.5 and 0.5, and -0.5 and 0.
#For the difference in means between logdose 0 and 0.5, the p-value=0.0504. Although the p-value is very close to 0.05, we conclude that there are no significant differences between these group means. This is more easily seen in the interaction plot.

#The residuals are approximately symmetric about zero and have approximately equal variance. In the Normal-QQ plot the residuals lie approximately along the reference line. The plots do not reveal any severe departure from the assumptions for this model.
#________________________________________________________________________

#2
#changing group into factors.
rm(age)
attach(vitcap2)
grou = factor(vitcap2$group)
vcap2 = data.frame(grou, age, vital.capacity)
str(vcap2)

L2=aov(vital.capacity ~ age * grou, data = vcap2)
anova(L2)

with(data = vcap2, expr = {
  interaction.plot(grou, ag, response = vital.capacity, lwd=2)
  interaction.plot(ag, grou, response = vital.capacity)
})

#We conclude there are significant differences in vital capacity between ages.
#There is also interaction between age and group.


model_L2=lm(vital.capacity ~ age * grou, data = vcap2)
summary(model_L2)
#In our model, the intercept is significant. The intercept corresponds to group 1 and age interacting with group one.
#Age, group 3, as well as age interacting with group 3 are all significant in predicting vital capacity.
TukeyHSD(L2, which = c("grou", "age"))

CIs = TukeyHSD(L2, which=1)
plot(CIs, las=1)
plot(TukeyHSD(L2, which = c("grou"))
)

#plot residuals vs fits
plot(L2, which=1:2)
#The residuals are approximately symmetric about zero and have approximately equal variance. In the Normal-QQ plot the residuals lie approximately along the reference line. The plots do not reveal any severe departure from the assumptions for this model.

#________________________________________________________________________
#3
attach(malaria)
#explanatory variables: age, log(ab)
#response: mal
M1=glm(mal ~  age + log(ab), family = binomial("logit"), data = malaria)
summary(M1)
#The intercept and log(ab) is significant. Our results make sense sense because people of any age can get maleria. Therefore, age is not a predictior of malaria.

#Null deviance corresponds to the deviance of a model that contains only the interval (and so, a fixed probability of success)

anova(M1,test="Chisq")	
#chi squared tests the difference between the null deviance and the residual deviance. Shows how our model is doing agains the null model (the model with only the intercept)

M2=glm(mal ~  age * log(ab), family = binomial("logit"), data = malaria)
summary(M2)
#We will use M1 with the + because AIC is lower.

#Age is not significant so we remove it from the model and use M3
M3=glm(mal ~ log(ab), family = binomial("logit"), data = malaria)
summary(M3)
#The estimate of alpha is 2.1552 and of beta1 is -0.7122. The odds of having malaria are less if the respondent has higher antibody count.
#The equation for logit is estimated by 2.1552-0.7122(log(ab))
anova(M3,test="Chisq")	

boxplot(log(ab) ~ mal, col = c(4,5), ylab = "log(ab)", xlab = c("Malaria"), xaxt = "n")
axis(side = 1, at = c(1, 2), labels = c("No", "Yes"))
#We notice 1 outlier with a low antibody count who doesnt have malaria and 2 outliers with high antibody counts who do have malaria. From the boxplot, it appears that those with malaria have a lower antibody count.


exp(cbind(OR=coef(M3),confint(M3)))
#The odds of malaria are estimated as 0.49 of what they would be if a person had one unit less of antibody level. Odds of having malaria increases with less antibodies.
#________________________________________________________________________
#4
#response variable: gvhd 0=no, 1=yes
attach(graft.vs.host)
rm(time)
M1=glm(gvhd ~ rcpage + donage + log(index), family	=	binomial("logit"))
summary(M1)	
anova(M1,test="Chisq")

M2=glm(gvhd ~ rcpage + log(index), family	=	binomial("logit"))
summary(M2)	
anova(M2,test="Chisq")

M3=glm(gvhd ~ donage + index, family	=	binomial("logit"))
summary(M3)	
anova(M3,test="Chisq")
#M3 has all variables significant and the lowest AIC

ntype = factor(graft.vs.host$type)
prego = factor(graft.vs.host$preg)
theydead = factor(graft.vs.host$dead)
ngraft = data.frame(pnr, rcpage, donage, ntype, prego, index, gvhd, time, theydead)
str(ngraft)

#The estimate of alpha is -6.086 and of beta1 is -.14075 an of beta 2 is 0.94439. The odds of having gvhd are more if the respondent has higher index and donor age.
#The equation for logit is estimated by -6.086+0.14075(donage)+0.94439(index)

#odds ratio
exp(cbind(OR=coef(M3),confint(M3)))
#The odds of graft vs. host are estimated as 1.15 times what they would be if a person was one year younger. Odds of having gvhd increases with older donor age.
boxplot(donage~gvhd)
#Increasing age increases the effect of graft vs. host

#The odds of graft vs. host are estimated as 2.57 times what they would be if a person had an index of 1 unit less. Odds of having gvhd increase with higher index. 
boxplot(index~gvhd)
