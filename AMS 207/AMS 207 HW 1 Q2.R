install.packages('LearnBayes')
library(LearnBayes)
cancermortality  


initial <- c(-6.82, 7.57)

H.star <- function(theta1, theta2, cancermortality)
{
  yj = cancermortality[,1]
  nj = cancermortality[,2]
  logf = function(yj, nj) (-1/20)*log(theta1+1000)-(1/20)*(lbeta((exp(theta1+theta2)/(1+exp(theta1)))+yj, exp(theta2)/(1+exp(theta1))+nj+yj))/lbeta(exp(theta1+theta2)/(1+exp(theta1)), exp(theta2)/(1+exp(theta1))-(1/20))*log((1+exp(theta1)^2)/exp(theta1)*(1+exp(theta2))^2)-(1/20)*log(exp(theta1+theta2)/(1+exp(theta1))^2)
}



optim(initial, H.star(-6.82, 7.57))


######Question 3########

t1=23721
t=15962989
suff_stat=c(t, t1)

logf = function(theta, suff_stat) {(theta[1]+theta[2]-8*theta[1])-((suff_stat[1]-8*suff_stat[2]+8*exp(theta[2]))/exp(theta[1]))}

logf1 = function(theta) {(theta[1]+theta[2]-8*theta[1])-((suff_stat[1]-8*suff_stat[2]+8*exp(theta[2]))/exp(theta[1]))}

#3
mycontour(logf(theta,suff_stat), c(13, 18, 4, 18), suff_stat, main="Contour Plot of LogPosterior", xlab='Theta1', ylab='Theta2')

itit_value=c(15, 12)

optim(itit_value,logf)

#4
#Normal approximation
npar=list(m=fit$mode,v=fit$var)
mycontour(lbinorm,c(13,18,4,18),npar,xlab="Theta1", ylab="Theta2",
          main = "Normal Approximation to the Posterior Distribution")

mycontour(logf(theta,suff_stat), c(13, 18, 4, 18), suff_stat, main="LogPosterior vs. Normal Approximation", xlab='Theta1', ylab='Theta2')
mycontour(lbinorm,c(13,18,4,18),npar, add=T, col='hotpink')
legend('bottomright', legend = c('Actual Posterior', 'Normal Approx'), lty=1, lwd=2, col=c('black', 'hotpink'))

#5
#Rejection sampling
#Changing to Multivariate-t
logfT=function(theta,suff_stat)
{
  d=logf(theta,suff_stat)-dmt(theta,mean=c(tpar$m),
                                S=tpar$var,df=tpar$df,log=TRUE)
  return(d)
}

lapt=laplace(logfT,c(15,13),suff_stat)
tmode=lapt$mode
tmax=logfT(tmode, suff_stat)

set.seed(3923)

tpar=list(m=fit$mode,var=2*fit$var,df=4)
RS<-rejectsampling(logf, tpar, tmax, 10000, suff_stat)
yay=dim(RS)
yay[1]/10000

mycontour(logf,c(13,18,4,18),suff_stat,xlab=expression(theta[1]),
          ylab=expression(theta[2]), main="Rejection Sampling")
points(RS[,1],RS[,2])

plot(hexbin(RS[,1],RS[,2]), colramp=rf, xlab=expression(theta[1]),
     ylab=expression(theta[2]), main="Rejection Sampling")

#Sampling Importance Resampling 
sir=sir(logf, tpar, 10000, suff_stat)

mycontour(logf,c(13,18,4,18),suff_stat,xlab=expression(theta[1]),ylab=expression(theta[2]), main="Sampling Importance Resampling")
points(sir[,1],sir[,2])

plot(hexbin(sir[,1],sir[,2]), colramp=rf, xlab=expression(theta[1]),
     ylab=expression(theta[2]), main="Sampling Importance Resampling")

#6
#Laplace
fit=laplace(logf, c(15, 13), suff_stat)
#posterior mode
mode=fit$mode
#variance-covariance matrix
fit$var

#Importance sampling
set.seed(273)
qpar=list(m=fit$mode,var=fit$var,df=4)
myfunc2=function(theta)return(theta[2])
s2=impsampling(logf,qpar,myfunc2,10000,suff_stat)
cbind(s2$est,10000*s2$se^2)

myfunc1=function(theta)return(theta[1])
s1=impsampling(logf,qpar,myfunc1,10000,suff_stat)
cbind(s1$est,s1$se^2)

#Mean and Variance
mu=s1$est
myvar1=function(theta)return((theta[1]-mu)^2)
svar1=impsampling(logf, qpar, myvar1, 10000, suff_stat)
print(svar1$est)

mu2=s2$est
myvar2=function(theta)return((theta[2]-mu2)^2)
svar2=impsampling(logf, qpar, myvar2, 10000, suff_stat)
print(svar2$est)

mycontour(logf,c(13,18,4,18),suff_stat,xlab=expression(theta[1]),ylab=expression(theta[2]), main="Importance Sampling")
points(s1$theta)


plot(hexbin(s1$theta), colramp=rf, xlab=expression(theta[1]),main="Importance Sampling", ylab=expression(theta[2]))

install.packages("RColorBrewer")
library(RColorBrewer)
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)

install.packages('hexbin')
library(hexbin)


#7
#Changing back to b & m
b=function(theta1)(exp(theta1))
m=function(theta2)(t1-exp(theta2))
nowitsb=b(sir[,1])
nowitsm=m(sir[,2])

realiability=exp(-(10^6 -nowitsm)/nowitsb)
mean(realiability)
var(realiability)
plot(density(realiability), main = expression('Reliability Distribution at'~ 10^6), xlab = 'Reliability')
abline(v=mean(realiability), col='hotpink')
text(0.3, 3.0, labels = 'Mean=0.622')
text(0.3, 2.5, labels='Variance=0.013')
