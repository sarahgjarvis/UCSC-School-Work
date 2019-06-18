library(MASS)
library("geoR")
library(rworldmap)
data <- read.csv("~/Downloads/sarah_data (1).csv")

# To perform this paert you need to retrieve the variable "Probability Threshold" from the NetCDF file. This corresponds to an idex that indicates the quality of the retrieval of the information. The idenx has ten possible values that corrrespond to probabilities: 0.9500 0.9000 0.8000 0.7000 0.6000 0.5000 0.4000 0.3000 0.2000 0.1000.

data$probthresh[data$probthresh == 0] <- 0.95
data$probthresh[data$probthresh == 1] <- 0.9
data$probthresh[data$probthresh == 2] <- 0.8
data$probthresh[data$probthresh == 3] <- 0.7
data$probthresh[data$probthresh == 4] <- 0.6
data$probthresh[data$probthresh == 5] <- 0.5
data$probthresh[data$probthresh == 6] <- 0.4
data$probthresh[data$probthresh == 7] <- 0.3
data$probthresh[data$probthresh == 8] <- 0.2
data$probthresh[data$probthresh == 9] <- 0.1

# Obtain new albedo values multiplying the observations by their corresponding probabilities. You now have two variables: the recorded albedo and the weighted albedo

data$w.albedo = data$bhriso * data$probthresh

# random sample of 500
# data = data[sample(nrow(data), size = 500),]

# explore data

hist(data$w.albedo, main = 'Weighted Albedo', xlab = 'weighted albedo', col = 'slateblue1')

# very skewed
# log normalizes

box.cox = function(data){
  bc <- boxcoxfit(data, lambda2 = TRUE)
  lambda <- bc$lambda[1]
  return((data^lambda - 1)/lambda)
}

y = log(data$w.albedo)
# y = log(data$bhriso)

# normal?
par(mfrow = c(1,2))
hist(y, xlab = 'log albedo', main = 'Histogram of Transformed Weighted Albedo', col = 'slateblue1')
qqnorm(y)
qqline(y, col = 'red')

# map
# plot of 500 randomly selected points
color.gradient <- function(x, colors=c("blue","purple","red"), colsteps=10) {
  return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
}
newmap <- getMap(resolution = "low")
{layout.matrix <- matrix(c(1, 2), nrow = 1, ncol = 2)
  
  layout(mat = layout.matrix,
         widths = c(5, 2))} 
plot(newmap, xlim = c(-90, -85), ylim = c(5, 45), asp = 1, main = 'Transformed Weighted Albedo Measurements')
points(data$longitude, data$latitude, col = color.gradient(y), pch = 20)
#legend
{plot(x = rep(max(data$w.albedo) + 0.67, 100),
      y = seq(min(data$w.albedo), max(data$w.albedo), length = 100),
      pch = 15, col = color.gradient(seq(0, 1, length = 100)),
      cex = 2.5, axes=F, xlab = '', ylab = '')
  
  lines(rep(max(data$w.albedo) + 0.97, 2), range(data$w.albedo))
  for (i in 1:6){
    lines(max(data$w.albedo) + c(0.97, 1.07), rep(min(data$w.albedo)+(i-1)*diff(range(data$w.albedo)/5), 2))
    #   text(max(loc[,1])+1.15, min(loc[,2])+(i-1)*diff(range(loc[,2])/5),
    #       round(quantile(y, (i-1)/5), 3), pos=4)
    text(max(data$w.albedo)+1, min(data$w.albedo)+(i-1)*diff(range(data$w.albedo)/5),
         round(seq(min(y), max(y), length = 6)[i], 2), pos=4, cex = 0.9)
  }
}


# explore trend
plot(data$longitude, y, main = 'Longitude', pch = 20, xlab = 'longitude', ylab = 'transformed weighted albedo')
plot(data$latitude, y, main = 'Latitude', pch = 20, xlab = 'latitude', ylab = 'transformed weighted albedo')

loc = cbind(data$longitude, data$latitude)

# Is there evidence of a first or second order trend function of location?
# explore trend
# fit regression with linear and quadratic terms in each direction. Select covariates with a stepwise AIC procedure.
mod = step(lm(y ~ . + .^2 + I(loc^2), data = data.frame(loc)), scope = list("lower" = lm(y ~ 1), "upper" = lm(y ~ . + .^2 + I(loc^2), data = data.frame(loc))), direction = "both")
summary(mod)

# Step:  AIC=-245.68
# y ~ longitude + latitude + I(longitude^2) + I(latitude^2)
yhat = predict(mod)

# 2nd order model we will use 
mod2 = lm(log(data$w.albedo) ~ data$longitude + data$latitude + I(data$longitude^2) + I(data$latitude^2) + data$longitude*data$latitude)
summary(mod2)

par(mfrow = c(2, 2))
plot(data$longitude, y, main = 'Longitude', pch = 20, xlab = 'longitude', ylab = 'transformed albedo')
points(data$longitude, yhat, col = 'slateblue1', pch = 20)
plot(data$latitude, y, main = 'Latitude', pch = 20, xlab = 'latitude', ylab = 'transformed albedo')
points(data$latitude, yhat, col = 'slateblue1', pch = 20)
plot(data$longitude, mod$residuals, main = 'Longitude', pch = 20, xlab = 'longitude', ylab = 'residuals')
abline(h = 0, col = 'red', lwd = 2)
plot(data$latitude, mod$residuals, main = 'Latitude', pch = 20, xlab = 'latitude', ylab = 'residuals')
abline(h = 0, col = 'red', lwd = 2)

# detrend
y = y - yhat


# binned variogram
plot(variog(data = y, coords = loc, uvec = 16), col = 'slateblue1', main='Binned Semi-Variogram for Weighted Albedo', pch = 20)

# Explore possible anisotropies using a directional variogram.
plot(variog4(data = y, coords = loc), main = 'Directional Variogram for the Residuals', legend = F)
legend('bottomleft', legend = c(expression(paste("0",degree), paste("45", degree), paste("90",degree), paste("135",degree))), lty = c(1,2,3,4), col = c(1,2,3,4), cex = 0.7)

# no evidence of anisotropies
# Use least squares to fit the covariograms in the Matèrn family with smoothness equal to .5; 1; 1.5; 2.5. Plot the results. Use the plots and the values of the LSE to select the best fit.
var = variog(data = y, coords = loc, messages = F)

smoothness = list(0.5, 1, 1.5, 2.5)
vf = lapply(smoothness, function(x) {variofit(var, ini.cov.pars = c(0.2,1), kappa = x)})

for (i in 1:length(smoothness)) {
  plot(var, main = paste(expression(kappa), '=', smoothness[i], ', SS = ', round(vf[[i]]$value, 3)))
  lines(vf[[i]], lwd = 2)
  segments(0, 0, x1=0, y1=vf[[i]]$nugget, col = 'red', lwd = 2)
  abline(v = vf[[i]]$practicalRange, col = 3, lty = 2, lwd = 2)
  abline(h = vf[[i]]$cov.pars[1] + vf[[i]]$nugget, col = 'blue', lwd = 2, lty = 4)
  segments(0, 0, x1=vf[[i]]$cov.pars[2], y1=0, col = 'purple', lwd = 2, cex = 0.8)
  legend('bottomright', legend = c('fitted', 'sill', 'nugget', 'range', 'practical range'), lty = c(2, 4, 1, 1, 2), bty = 'n', cex = 0.8, col = c(1, 4, 2, 'purple', 3), lwd = 2)
}

# the best fit has the lowest SS (sum of squared loss)
vf.bestfit = vf[[which.min(sapply(vf, function(x){(x$value)}))]]

# MLE fit
mle = likfit(geodata, trend = '2nd', cov.model = 'matern', kappa = 1.5, ini.cov.pars = c(0.2,1))
summary(mle)

# MCMC
# not sure if this works use function later
{
D = as.matrix(mod$model)
k = ncol(D)
n = nrow(D)

dists = as.matrix(dist(loc, upper = TRUE, diag= TRUE))
K.fun = function(psi, gamma2, kappa){
  diag(length(bc)) + (1/gamma2) * geoR::matern(dists, psi, kappa)
}

# burnin
T = 10
# n samples
N = 1100


a = b = c = d = e = f = 1

# create parameters
beta = matrix(0, nrow = N, ncol = k)
tau2 = psi = gamma2 = kappa = rep(NA, N)

# initialize parameters
tau2[1] = psi[1] = gamma2[1] = 1

smoothness = c(0.5, 1, 1.5, 2.5)

kappa[1] = smoothness[2]


new.vals = function(psi, gamma2, kappa){
  K = K.fun(psi = psi, gamma2 = gamma2, kappa = kappa)
  K.inv = solve(qr.R(qr(K))) %*% t(qr.Q(qr(K)))
  DK.invD = t(D) %*% K.inv %*% as.matrix(D)
  DK.invD.inv = solve(qr.R(qr(DK.invD))) %*% t(qr.Q(qr(DK.invD)))
  b.hat = DK.invD.inv %*% (t(D) %*% K.inv %*% y)
  resid = y - (D %*% as.matrix(b.hat))
  S2 = t(resid) %*% K.inv %*% resid
  
  return(list(K = K, K.inv = K.inv, DK.invD = DK.invD, DK.invD.inv = DK.invD.inv, b.hat = b.hat, S2 = S2))
}
accept.count = 0
for (i in 2:N) {
  vals = new.vals(psi[i-1], gamma2[i-1], kappa[1])
  beta[i,] = as.matrix(mvrnorm(1, vals$b.hat, tau2[i-1] * vals$DK.invD))
  tau2[i] = 1/rgamma(1, (n-k)/2 + a, vals$S2/2 + b)
  
  # block sampling step/metropolis-hastings step?
  proposed = c(mvrnorm(1, c(psi[i-1], gamma2[i-1]), Sigma = 1e-3 * diag(2)))
  
  proposed.vals = new.vals(proposed[1], proposed[2], kappa[1])
  
  post.proposed = -0.5*determinant(proposed.vals$K, logarithm = T)$modulus[1] - 0.5*determinant(proposed.vals$DK.invD, logarithm = T)$modulus[1] - ((n+k)/2 + a)*log(proposed.vals$S2 + 2*b) + dgamma(proposed[1], c, d, log = T) + dgamma(proposed[2], e, f, log = T)
  
  post.current = -0.5*determinant(vals$K, logarithm = T)$modulus[1] - 0.5*determinant(vals$DK.invD, logarithm = T)$modulus[1] - ((n+k)/2 + a)*log(vals$S2 + 2*b) + dgamma(psi[i-1], c, d, log = T) + dgamma(gamma2[i-1], e, f, log = T)
  
  # acceptance ratio
  if(post.proposed - post.current > log(runif(1))){
    accept.count = accept.count + 1
    psi[i] = proposed[1]
    gamma2[i] = proposed[2]
  }
  else
    psi[i] = psi[i-1]
    gamma2[i] = gamma2[i-1]
  print(i)
}

plot.ts(beta[100:1100,])

plot.ts(beta[1:114,])
plot.ts(tau2[5:468])

colMeans(beta[100:1100,])
summary(mod)
accept.count/1000}

# priors
PC = prior.control(beta.prior = "flat", 
                              sigmasq.prior = "reciprocal", 
                              phi.prior = "squared.reciprocal", 
                              tausq.rel.prior = "reciprocal", 
                              tausq.rel.discrete = seq(.2, 9, length.out = 40), 
                              phi.discrete = seq(.01, 8, length.out = 40))

# model
MC = model.control(trend.d = '2nd', trend.l = '2nd', cov.model = 'matern', kappa = 2.5)

# output
OC = output.control(n.post = 1000, moments = T)

# data
geodata = list(y, loc)
names(geodata) = c('data', 'coords')

# defining the grid of prediction locations:
grid <- as.matrix(expand.grid(seq(-90, -50, l=30), seq(10.5, 45, l=30)))

#make grid of data only on land!!!!
library('maptools')
data(wrld_simpl)
pts <- SpatialPoints(grid, proj4string=CRS(proj4string(wrld_simpl)))
## Find which points fall over land
land <- pts[!is.na(over(pts, wrld_simpl)$FIPS)]@coords

bayes.fit = krige.bayes(geodata, model=MC, locations = land,
                        prior = PC, output = OC, borders = NULL)
bayes.fit2 = krige.bayes(geodata, model=MC, locations = loc,
                        prior = PC, output = OC, borders = NULL)

estimates = lapply(bayes.fit$posterior$sample, function(x) quantile(x, c(0.025, 0.5, 0.975)))

par(mar=c(4,2,1,1))
par(mfrow = c(3,3))
hist(bayes.fit)

par(mfrow = c(1,2))
plot(bayes.fit, col=2:3)

predict = bayes.fit$predictive$mean

# plot predictions
color.gradient <- function(x, colors=c("blue","purple","red"), colsteps=10) {
  return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
}
newmap <- getMap(resolution = "low")
{layout.matrix <- matrix(c(1, 2), nrow = 1, ncol = 2)
  
  layout(mat = layout.matrix,
         widths = c(5, 2))} 
plot(newmap, xlim = c(-90, -85), ylim = c(5, 45), asp = 1, main = 'Predicted Weighted Albedo Measurements')
points(land[,1], land[,2], col = color.gradient(predict), pch = 16, cex = 1.65)
#legend
{plot(x = rep(max(predict), 100),
      y = seq(min(predict), max(predict), length = 100),
      pch = 15, col = color.gradient(seq(0, 1, length = 100)),
      cex = 1.5, axes=F, xlab = '', ylab = '')
  
  lines(rep(max(predict)+.2, 2), range(predict))
  for (i in 1:6){
    lines(max(predict) + c(.2,.3), rep(min(predict)+(i-1)*diff(range(predict)/5), 2))
    #   text(max(loc[,1])+1.15, min(loc[,2])+(i-1)*diff(range(loc[,2])/5),
    #       round(quantile(y, (i-1)/5), 3), pos=4)
    text(max(predict)+.2, min(predict)+(i-1)*diff(range(predict)/5),
         round(seq(min(predict), max(predict), length = 6)[i], 2), pos=4, cex = 0.7)
  }
}

# Model check

pred.int = apply(bayes.fit$predictive$simulations, 1, quantile, c(.025,.975))
ind=1:25
matplot(rbind(ind,ind), pred.int[,1:25], type="l", lty=1,col=1, xlab="Observation",ylab="log weighted albedo", main = '95% Intervals for Weighted Albedo')
points(ind, y[1:25], pch=19)

1-(sum(y > pred.int[2,] | y < pred.int[1,])/500)
# 71% accuracy


# non weighted albedo 
pred.int = apply(bayes.fit2$predictive$simulations, 1, quantile, c(.025,.975))
matplot(rbind(ind,ind), pred.int[,1:25], type="l", lty=1,col=1, xlab="Observation",ylab="log albedo", main = '95% Intervals for Albedo')
points(ind, y[1:25], pch=19)

1-(sum(y > pred.int[2,] | y < pred.int[1,])/500)
