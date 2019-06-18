attach(sst)
library(invgamma)
library(MASS)
pairs(~ temp + lat + lon + Type)
str(sst)
library(mvtnorm)

install.packages('GGally')
library(GGally)

ggpairs(sst[, c("lat", "lon", "temp", "N", "Type")], aes(color=Type))

fit <- lm(temp ~ lat + lon + Type, x = T)
summary(fit)
X<-model.matrix(temp ~ lat + lon +Type, sst)
xtx <- t(X)%*%X
invxtx <- solve(xtx)

#Gibbs Sampler
T <- 1000 #Burnin
n <- 86
mui = matrix(0, 11000, n)
mui[1,] = temp #initialize

sigi = matrix(0, 11000, n)
sigi[1,] = 1 #initialize

beta = matrix(0, 11000, 6)
beta[1,] = c(fit$coefficients)

sig <- tau <- rep(NA, 11000)
sig[1] <- tau[1] <- 1 #initialize

alpha <- 3
a <- 2
b <- 2
temp <- sst$temp
for (i in 2:11000)
{  
  
  V = tau[i-1] * invxtx
  M = invxtx %*% t(X) %*% mui[i-1,]
  beta[i,] <- rmvnorm(1, mean=M, sigma = V)
  
  tempo = mui[i-1,] - (X %*% beta[i,])
  tau[i] = rinvgamma(1, (n/2)-1, (t(tempo) %*% tempo)/2)
  
  v = (N/sigi[i-1,] + 1/tau[i])^(-1)
  m = v*(temp*N/sigi[i-1,] + (X %*% beta[i,])/tau[i])
  mui[i,] <- rnorm(86, m, sqrt(v))

  sigi[i,] = rinvgamma(n, alpha+1.5, alpha*sig[i-1]+(N*(temp-mui[i,])^2)/2)
  
  sig[i] = rgamma(1, n*alpha+n+a, b + alpha * sum(1/sigi[i,]))
  
}
 #YAY#




################################ Remove burnin ONLY RUN ONCE ###################################
beta <- beta[-(1:T),]
mui <- mui[-(1:T),]
sigi <- sigi[-(1:T),]
tau <- tau[-(1:T)]
sig <- sig[-(1:T)]
################################################################################################

######Trace plots#####

par(mfrow=c(2,3))
for(i in 1:6)
{plot(ts(beta[,i]), ylab=expression(beta[i]), xlab=i)}
mtext(expression(beta[i]~'Samples'), outer = T, line = -3)
dev.off()

#One column of mui and sigi samples
par(mfrow=c(2,1))
plot(ts(mui[,86]), ylab=expression(mu[86]), xlab='iteration')
plot(ts(sigi[,86]), ylab=expression(sigma[86]^2), xlab='iteration')

par(mfrow=c(2,1))
plot(ts(tau), ylab=expression(tau^2), xlab='iteration')
plot(ts(sig), ylab=expression(sigma^2), xlab='iteration')

##Posterior Distributions##
par(mfrow=c(2,2))
hist(mui[,86],  xlab=expression(mu[86]), border = 'paleturquoise4', col = 'paleturquoise2', main = NULL)
hist(sigi[,86], xlab=expression(sigma[86]^2), border = 'paleturquoise4', col = 'paleturquoise2', main = NULL)
hist(tau, border = 'paleturquoise4', col = 'paleturquoise2', xlab=expression(tau^2), main = NULL)
hist(sig, border = 'paleturquoise4', col = 'paleturquoise2', xlab=expression(sigma^2), main = NULL)
mtext('Gibbs Samples', outer = T, line = -3)


par(mfrow=c(2,3))
for (i in 1:6) {
  hist(beta[,i], border = 'paleturquoise4', col = 'paleturquoise2', main = NULL, xlab = i)
}
mtext(expression(beta[i]~'Samples'), outer = T, line = -3)
dev.off()

######Null Model######
x.null <- as.matrix(rep(1,86))
xtx.null <- t(x.null)%*%x.null
invxtx.null <- solve(xtx.null)

mui.null = matrix(0, 11000, n)
mui.null[1,] = temp #initialize

sigi.null = matrix(0, 11000, n)
sigi.null[1,] = 1 #initialize

beta.null = matrix(0, 11000, 1)
beta.null[1] = 1#c(fit$coefficients)

sig.null <- tau.null <- rep(NA, 11000)
sig.null[1] <- tau.null[1] <- 1 #initialize

for (i in 2:11000)
{  
  
  V = tau.null[i-1] * invxtx.null
  M = invxtx.null %*% t(x.null) %*% mui.null[i-1,]
  beta.null[i,] <- rmvnorm(1, mean=M, sigma = V)
  
  tempo.null = mui.null[i-1,] - (x.null %*% beta.null[i,])
  tau.null[i] = rinvgamma(1, (n/2)-1, (t(tempo.null) %*% tempo.null)/2)
  
  v = (N/sigi.null[i-1,] + 1/tau.null[i])^(-1)
  m = v*(temp*N/sigi.null[i-1,] + (x.null %*% beta.null[i,])/tau.null[i])
  mui.null[i,] <- rnorm(86, m, sqrt(v))
  
  sigi.null[i,] = rinvgamma(n, alpha+1.5, alpha*sig.null[i-1]+(N*(temp-mui.null[i,])^2)/2)
  
  sig.null[i] = rgamma(1, n*alpha+n+a, b + alpha * sum(1/sigi.null[i,]))
  print(i)
}


################################ Remove burnin ONLY RUN ONCE ###################################
beta.null <- beta.null[-(1:T),]
mui.null <- mui.null[-(1:T),]
sigi.null <- sigi.null[-(1:T),]
tau.null <- tau.null[-(1:T)]
sig.null <- sig.null[-(1:T)]
################################################################################################

#Trace plots
#One column of mui and sigi samples
par(mfrow=c(2,1))
plot(ts(mui.null[,86]), ylab=expression(mu[86]~'.null'), xlab='iteration')
plot(ts(sigi.null[,86]), ylab=expression(sigma[86]~'.null'^2), xlab='iteration')

par(mfrow=c(3,1))
plot(ts(beta.null), ylab=expression(beta[null]), xlab='iteration')
plot(ts(tau.null), ylab=expression(tau^2), xlab='iteration')
plot(ts(sig.null), ylab=expression(sigma^2), xlab='iteration')

### null posterior distributions ###
par(mfrow=c(2,2))
hist(mui.null[,86],  xlab=expression(mu[86]), border = 'palegreen4', col = 'palegreen2', main = NULL)
hist(sigi.null[,86], xlab=expression(sigma[86]^2), border = 'palegreen4', col = 'palegreen2', main = NULL)
hist(tau.null, border = 'palegreen4', col = 'palegreen2', xlab=expression(tau^2), main = NULL)
hist(sig.null, border = 'palegreen4', col = 'palegreen2', xlab=expression(sigma^2), main = NULL)
mtext('Null Model Samples', outer = T, line = -3)

hist(beta.null, border = 'palegreen4', col = 'palegreen2', main = expression(beta~'Samples'), xlab = expression(beta))
dev.off()


### G & G to compare the models###
pred.draws.m2=matrix(0, n, 10000)
for (i in 1:10000)
{
  pred.draws.m2[,i]=rnorm(n, mui.null[i,], sqrt(sigi.null[i,]/N))
}

### G&G ###
print(gg.m1 <- sum((temp - apply(pred.draws,1,mean))^2) + sum(apply(pred.draws, 1, var)))
print(gg.m2 <- sum((temp - apply(pred.draws.m2,1,mean))^2) + sum(apply(pred.draws.m2, 1, var)))



### DIC M1 ###
mui.hat <- apply(mui, 2, mean)
sigi.hat <- apply(sigi, 2, mean)

lnormlike=function(y, mu, sigi){
  val=sum(dnorm(y, mu, sigi, log=T))
  return(val)
}

hlp<-NULL
for(t in 1:10000){
  hlp[t]<-lnormlike(temp, mui[t,], sqrt(sigi[t,]/N))
}
lph<-lnormlike(temp, mui.hat, sqrt(sigi.hat/N))
pdic=2*(lph-mean(hlp))
print(DIC.m1 <- -2*lph+2*pdic)

### DIC M2 ###
mui.hat.null <- apply(mui.null, 2, mean)
sigi.hat.null <- apply(sigi.null, 2, mean)

hlp<-NULL
for(t in 1:10000){
  hlp[t]<-lnormlike(temp, mui.null[t,], sqrt(sigi.null[t,]/N))
}
lph<-lnormlike(temp, mui.hat.null, sqrt(sigi.hat.null/N))
pdic=2*(lph-mean(hlp))
print(DIC.m1 <- -2*lph+2*pdic)

####Posterior Predictive Check####
pred.draws=matrix(0, n, 10000)
for (i in 1:10000)
{
  pred.draws[,i]=rnorm(n, mui[i,], sqrt(sigi[i,]/N))
}

pred.int=apply(pred.draws, 1, quantile, c(.05,.95))


ind=1:86
matplot(rbind(ind,ind), pred.int, type="l", lty=1, xlab="Observation",ylab="SST", col='darkmagenta')
points(ind, temp, pch=19, cex=.4)

par(mfrow=c(2, 2))
predict.mean <- apply(pred.draws, 2, mean)
hist(predict.mean, border = 'darkmagenta', main = NULL, xlab = 'mean')
abline(v=mean(temp), lty=2, lwd=2)
length(which(predict.mean>mean(temp))==TRUE)/11000
text(20, 2500, labels = 'p=0.513')

pred.sd <- apply(pred.draws, 2, sd)
hist(pred.sd, border  = 'darkmagenta', main = NULL, xlab = 'sd')
abline(v=sd(temp), lty=2, lwd=2)
length(which(pred.sd>sd(temp))==TRUE)/11000
text(2, 2500, labels = 'p=0.573')

pred.min <- apply(pred.draws, 2, min)
hist(pred.min, border = 'darkmagenta', main = NULL, xlab = 'min')
abline(v=min(temp), lty=2, lwd=2)
length(which(pred.min>min(temp))==TRUE)/11000
text(10, 2000, labels = 'p=0.803')

pred.max <- apply(pred.draws, 2, max)
hist(pred.max, border = 'darkmagenta', main = NULL, xlab = 'max')
abline(v=max(temp), lty=2, lwd=2)
length(which(pred.max>max(temp))==TRUE)/11000
text(28, 2000, labels = 'p=0.946')
dev.off()

#Posterior predictive samples
for(i in 1:30){lines(density(pred.draws[,i]), col='lightpink', main=NULL)}
lines(density(temp), col='darkmagenta', main="Data Distribution")
legend('topleft', legend = c('predictive sample', 'data'), cex=0.8, col = c('lightpink', 'darkmagenta'), lty=1)

##### MAP THING #####
install.packages('ggmap')
library(ggplot2)
library(ggmap)

### Part 5 ###
#Consider a rectangle area that contains all the observational sites. Explore the posterior predictive dist of the SST on the grid.
lat.seq <- seq(min(lat), max(lat), by=0.2)
lon.seq <- seq(min(lon), max(lon), by=0.2)
grid <- expand.grid(lon.seq, lat.seq)
dev.type <- matrix(c(0,0,0,1,0,0,0,1,0,0,0,1), nrow = 4, ncol = 3, byrow = T)
data = merge(grid, dev.type, by = NULL)
data = cbind(rep(1,nrow(grid)),data)
n.l.l = nrow(data)

n.pred = 5000
#simulate new stuff from prior with samples
sigi.pred = matrix(NA, n.l.l, n.pred)
for (i in 1:n.l.l)
{sigi.pred[i,] = rinvgamma(n.pred, alpha+1, alpha*sig[1:n.pred])}

beta.pred = cbind(beta[1:n.pred,1], beta[1:n.pred,3], beta[1:n.pred,2], beta[1:n.pred,4], beta[1:n.pred,5],beta[1:n.pred,6])

mean.mu.i = as.matrix(data)%*%t(as.matrix(beta.pred))

mui.pred <-  matrix(0, n.l.l, n.pred)
for (i in 1:n.l.l) 
  {
  mui.pred[i,] = rnorm(n.pred, mean.mu.i[i,], sqrt(tau[1:n.pred]))
}


new.y.pred <-  matrix(0, n.l.l, n.pred)
for (i in 1:n.l.l) 
  {new.y.pred[i,] = rnorm(n.pred, mui.pred[i,], sqrt(sigi.pred[i,]))}


new.y.pred.mean = apply(new.y.pred, 1, mean)
new.y.pred.sd = apply(new.y.pred, 1, sd)
#say why just 1 n we assume

l = list(as.matrix(new.y.pred.mean), as.matrix(data[,-1]))
data.new = do.call(cbind,l)
colnames(data.new) <- c("y", "lon", "lat", "dev1", "dev2", "dev3" )
data.new = data.frame(data.new)

data.new.dev1 = data.new[data.new$dev1 == "0" & data.new$dev2 == "0" & data.new$dev3 == "0",]
data.new.dev2 = data.new[data.new$dev1 == "1" & data.new$dev2 == "0" & data.new$dev3 == "0",]
data.new.dev3 = data.new[data.new$dev1 == "0" & data.new$dev2 == "1" & data.new$dev3 == "0",]
data.new.dev4 = data.new[data.new$dev1 == "0" & data.new$dev2 == "0" & data.new$dev3 == '1',]


map <- get_map(location = c(lon = mean(sst$lat), lat = mean(sst$lon)), zoom = 5,
               maptype = "satellite", source = "google")
library(rworldmap)

ggmap(map) +
  geom_point(aes(x = lon, y = lat, size = temp,
                 colour = temp), data = sst.data, alpha = .5)


grid.arrange(p1, p2, p3, p4, nrow = 2)
ggplot(data = data.new.dev1, aes(x = lon, y = lat)) +
  scale_fill_gradient(low = "white",high = "aquamarine1") +
  geom_tile(aes(fill = y))+ labs(title = "Bucket")

ggplot(data = data.new.dev2, aes(x = lon, y = lat)) +
scale_fill_gradient(low = "white",high = "aquamarine1") +
geom_tile(aes(fill = y))+ labs(title = "d.buoy")

ggplot(data = data.new.dev3, aes(x = lon, y = lat)) +
  scale_fill_gradient(low = "white",high = "aquamarine1") +
  geom_tile(aes(fill = y))+ labs(title = "eri")

ggplot(data = data.new.dev4, aes(x = lon, y = lat)) +
  scale_fill_gradient(low = "white",high = "aquamarine1") +
  geom_tile(aes(fill = y))+ labs(title = "f.buoy")

###PLOT PIC ONTOP OF MAP###
require(png)
par(mfrow=c(2,2))
newmap <- getMap(resolution = "low")
plot(newmap, xlim = c(20, 35), ylim = c(30, 40), asp = 1, main='f.buoy')
img<-readPNG("/Users/SarJar/Desktop/f.png")
rasterImage(img, 25.1, 31.8, 33.8, 34.5)

q = list(as.matrix(new.y.pred.sd), as.matrix(data[,-1]))
data.new = do.call(cbind,q)
colnames(data.new) <- c("sd", "lon", "lat", "dev1", "dev2", "dev3" )
data.new = data.frame(data.new)

data.new.dev1 = data.new[data.new$dev1 == "0" & data.new$dev2 == "0" & data.new$dev3 == "0",]
data.new.dev2 = data.new[data.new$dev1 == "1" & data.new$dev2 == "0" & data.new$dev3 == "0",]
data.new.dev3 = data.new[data.new$dev1 == "0" & data.new$dev2 == "1" & data.new$dev3 == "0",]
data.new.dev4 = data.new[data.new$dev1 == "0" & data.new$dev2 == "0" & data.new$dev3 == '1',]

grid.arrange(p1, p2, p3, p4, nrow = 2)
ggplot(data = data.new.dev1, aes(x = lon, y = lat)) +
  scale_fill_gradient(low = "white",high = "orange2") +
  geom_tile(aes(fill = sd^2))+ labs(title = "Bucket")

ggplot(data = data.new.dev2, aes(x = lon, y = lat)) +
  scale_fill_gradient(low = "white",high = "orange2") +
  geom_tile(aes(fill = sd^2))+ labs(title = "d.buoy")

ggplot(data = data.new.dev3, aes(x = lon, y = lat)) +
  scale_fill_gradient(low = "white",high = "orange2") +
  geom_tile(aes(fill = sd^2))+ labs(title = "Eri")

ggplot(data = data.new.dev4, aes(x = lon, y = lat)) +
  scale_fill_gradient(low = "white",high = "orange2") +
  geom_tile(aes(fill = sd^2))+ labs(title = "f.buoy")
mean(sig)
sd(sig)
