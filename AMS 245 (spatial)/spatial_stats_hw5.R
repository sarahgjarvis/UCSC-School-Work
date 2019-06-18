library(ncdf4)
library('maptools')
library(rworldmap)
library(gtools)
library(maps)
library(geoR)
library(spBayes)
library(fields)
library(maptools)

# get data
nc = nc_open('/Users/SarJar/Downloads/GSA_AlbedoProd_GOES_075_VIS02_2000_181.nc')
#nc = nc_open('/Users/SarJar/Downloads/GSA_AlbedoProd_GOES_135_VIS02_2000_181.nc')
lon = as.vector(ncvar_get(nc, 'longitude'))
lat = as.vector(ncvar_get(nc, 'latitude'))
albedo = as.vector(ncvar_get(nc, 'BHRiso'))
temp = data.frame(cbind(lon, lat, albedo))
temp = temp[complete.cases(temp),]

# subset my region
temp = temp[temp$lon >= -126 & temp$lon <= -90 & temp$lat >= 10.5 & temp$lat <= 34, ]

## Find which points fall over land
data(wrld_simpl)
pts <- SpatialPoints(temp, proj4string=CRS(proj4string(wrld_simpl)))
temp <- pts[!is.na(over(pts, wrld_simpl)$FIPS)]@coords

# plot
quilt.plot(x=temp[,1], y=temp[,2],
           z=temp[,3], main = 'Albedo')
map("world", add = TRUE)
quilt.plot(x=goes75[,1], y=goes75[,2],
           z=goes75[,3], main = 'Albedo')
map("world", add = TRUE)


quilt.plot(x=temp[,1], y=temp[,2],
           z=temp[,3], main = 'Albedo')
map("world", add = TRUE)
quilt.plot(x=goes135[,1], y=goes135[,2],
           z=goes135[,3], main = 'Albedo')
map("world", add = TRUE)


# testing set
test = goes75[sample(nrow(goes75), size = 1000),]

goes75 = goes75[sample(nrow(goes75), size = 1000),]
goes135 = goes135[sample(nrow(goes135), size = 1000),]

goes75 = data.frame(goes75)
goes135 = data.frame(goes135)

# knots (where to predict)
grid <- as.matrix(expand.grid(seq(-126, -90, l=10), seq(10.5, 34, l=10)))

# land only
data(wrld_simpl)
pts <- SpatialPoints(grid, proj4string=CRS(proj4string(wrld_simpl)))
knots <- pts[!is.na(over(pts, wrld_simpl)$FIPS)]@coords


#---------------------------------------------------------------------------
# explore data

hist(goes75$albedo)
hist(goes135$albedo)
# skewed - try logit because albedo is in (0,1)
hist(log(goes75$albedo))
hist(log(goes135$albedo))

plot(density(log(goes75$albedo)))
plot(density(log(goes135$albedo)))

qqnorm(log(goes75$albedo))
qqline(log(goes75$albedo), col = 'red')

qqnorm(log(goes135$albedo))
qqline(log(goes135$albedo), col = 'red')

y1 = log(goes75$albedo)
y2 = log(goes135$albedo)

plot(goes75)
plot(goes135)


quilt.plot(x=goes75[,1], y=goes75[,2],
           z=y1, main = 'GOES 075 log Albedo')
points(knots[,1], knots[,2], pch = 16)
map("world", add = TRUE)

quilt.plot(x=goes135[,1], y=goes135[,2],
           z=y2, main = 'GOES 135 log Albedo', xlim = c(-118, -89.8), ylim = c(15, 34.2))
points(knots[,1], knots[,2], pch = 16)
map("world", add = TRUE)

#------------------------ trend ------------------------

plot(goes75$lon, y1)
plot(goes75$lat, y1)
plot(goes135$lon, y2)
plot(goes135$lat, y2)

loc1 = goes75[,1:2]
loc2 = goes135[,1:2]

geodata1 = list(y1, loc1)
geodata2 = list(y2, loc2)
names(geodata2) = c('data', 'coords')

mod = step(lm(y2 ~ . + .^2 + I(loc2^2), data = data.frame(loc2)), scope = list("lower" = lm(y2 ~ 1), "upper" = lm(y2 ~ . + .^2 + I(loc2^2), data = data.frame(loc2))), direction = "both")
summary(mod)
AIC(mod)

mod1 = lm(y1 ~ lon + lat + I(lon^2) + I(lat^2) + lon*lat, data = loc1)
summary(mod1)
AIC(mod1) # 1348.827

mod2 = lm(y2 ~ lon + lat + I(lon^2) + I(lat^2) + lon*lat, data = loc2)
summary(mod2)
AIC(mod2) # 1193.431

mle1 = likfit(geodata1, trend = '2nd', cov.model = 'matern', kappa = 2.5, ini.cov.pars = c(0.2,1))
summary(mle1)
# 0.5, AIC = 656.5, BIC = 700.7, nugget = 0.0176, phi = 1.218, sigsq = 0.2677

mle2 = likfit(geodata2, trend = '2nd', cov.model = 'matern', kappa = .5, ini.cov.pars = c(0.2,1))
summary(mle2)
# 0.5, AIC = 537.6, BIC = 581.8, nugget = 0.0295, phi = 1.387, sigsq = 0.2179


#---------------------- predictive process GOES 075 ----------------
X = cbind(loc1$lon, loc1$lat, I(loc1$lon^2), I(loc1$lat^2), loc1$lon*loc1$lat)


priors = list("beta.flat", "sigma.sq.ig" = c(1, 0.5), "tau.sq.ig" = c(1, 0.01),
              "phi.Unif" = c(1e-6, 5), "nu.Unif" = c(0.1, 3.5))
starting = list('beta' = mod1$coefficients, "sigma.sq" = 0.2677, "tau.sq" = 0.017,
                "phi" = 1/1.2, "nu" = 0.5)
tuning = list("sigma.sq" = 0.05, "tau.sq" = 0.05, "phi" = 0.05, "nu" = 0)


pp1 = spLM(y1 ~ X, coords = X[,1:2], knots = knots, n.samples = 5000, cov.model = 'matern', modified.pp = FALSE, priors = priors, starting = starting, tuning = tuning)

# modified GP
pp1.mod = spLM(y1 ~ X, coords = X[,1:2], knots = knots, n.samples = 5000, cov.model = 'matern', modified.pp = TRUE, priors = priors, starting = starting, tuning = tuning)

traceplot(pp1$p.theta.sample)
colMeans(pp1$p.theta.sample[500:5000,])

traceplot(pp1.mod$p.theta.sample)
colMeans(pp1.mod$p.theta.sample[500:5000,])

#---------------------- predictive process GOES 135 ----------------
X2 = cbind(loc2$lon, loc2$lat, I(loc2$lon^2), I(loc2$lat^2), loc2$lon*loc2$lat)


priors = list("beta.flat", "sigma.sq.ig" = c(1, 0.5), "tau.sq.ig" = c(1, 0.01),
              "phi.Unif" = c(1e-6, 5), "nu.Unif" = c(0.1, 3.5))
starting = list('beta' = mod2$coefficients, "sigma.sq" = 0.2179, "tau.sq" = 0.017,
                "phi" = 1/1.387, "nu" = 0.5)
tuning = list("sigma.sq" = 0.05, "tau.sq" = 0.05, "phi" = 0.05, "nu" = 0)


pp2 = spLM(y2 ~ X2, coords = X2[,1:2], knots = knots, n.samples = 5000, cov.model = 'matern', modified.pp = FALSE, priors = priors, starting = starting, tuning = tuning)

# modified GP
pp2.mod = spLM(y2 ~ X2, coords = X2[,1:2], knots = knots, n.samples = 5000, cov.model = 'matern', modified.pp = TRUE, priors = priors, starting = starting, tuning = tuning)

traceplot(pp1$p.theta.sample)
colMeans(pp1$p.theta.sample[500:5000,])

traceplot(pp1.mod$p.theta.sample)
colMeans(pp1.mod$p.theta.sample[500:5000,])

traceplot(pp2$p.theta.sample)

#----------------------- predict ---------------------

X.test = cbind(1, test[,1], test[,2], I(test[,1]^2), I(test[,2]^2), test[,1]*test[,2])
pp1.samples = spPredict(pp1, pred.coords = test[,1:2], pred.covars = X.test)
pp1.mod.samples = spPredict(pp1.mod, pred.coords = test[,1:2], pred.covars = X.test)

pred.pp1 = apply(pp1.samples$p.y.predictive.samples, 1, mean)
pred.pp1.mod = apply(pp1.mod.samples$p.y.predictive.samples, 1, mean)

pp2.samples = spPredict(pp2, pred.coords = test[,1:2], pred.covars = X.test)
pp2.mod.samples = spPredict(pp2.mod, pred.coords = test[,1:2], pred.covars = X.test)

pred.pp2 = apply(pp2.samples$p.y.predictive.samples, 1, mean)
pred.pp2.mod = apply(pp2.mod.samples$p.y.predictive.samples, 1, mean)

#-------------------- accuracy ------------------
 
pred.int = apply(pp1.samples$p.y.predictive.samples, 1, function(x) quantile(x, c(0.05, 0.95)))
pred.int.mod = apply(pp1.mod.samples$p.y.predictive.samples, 1, function(x) quantile(x, c(0.05, 0.95)))

pred.int2 = apply(pp2.samples$p.y.predictive.samples, 1, function(x) quantile(x, c(0.05, 0.95)))
pred.int.mod2 = apply(pp2.mod.samples$p.y.predictive.samples, 1, function(x) quantile(x, c(0.05, 0.95)))

# proportion of observatins that fall inside the prediction interval
1-(sum(log(test[,3]) > pred.int[2,] | log(test[,3]) < pred.int[1,])/nrow(test))
1-(sum(log(test[,3]) > pred.int.mod[2,] | log(test[,3]) < pred.int.mod[1,])/nrow(test))

1-(sum(log(test[,3]) > pred.int2[2,] | log(test[,3]) < pred.int2[1,])/nrow(test))
1-(sum(log(test[,3]) > pred.int.mod2[2,] | log(test[,3]) < pred.int.mod2[1,])/nrow(test))

plot(density(pp1.samples$p.y.predictive.samples[1,]))
abline(v = log(test[1,3]))

plot(density(pp1.samples$p.y.predictive.samples[100,]))
abline(v = log(test[100,3]))

plot(density(pp1.samples$p.y.predictive.samples[500,]))
abline(v = log(test[500,3]))

ind=1:25
par(mfrow = c(2,2))

matplot(rbind(ind,ind), pred.int[,1:25], type="l", lty=1,col=1, xlab="Observation",ylab="log albedo", main = '90% Prediction Intervals for GOES 075')
points(ind, log(test[1:25,3]), pch=19)

matplot(rbind(ind,ind), pred.int.mod[,1:25], type="l", lty=1,col=1, xlab="Observation",ylab="log albedo", main = '90% Prediction Intervals for GOES 075 (Modified GP)')
points(ind, log(test[1:25,3]), pch=19)

matplot(rbind(ind,ind), pred.int2[,1:25], type="l", lty=1,col=1, xlab="Observation",ylab="log albedo", main = '90% Prediction Intervals for GOES 135')
points(ind, log(test[1:25,3]), pch=19)

matplot(rbind(ind,ind), pred.int.mod2[,1:25], type="l", lty=1,col=1, xlab="Observation",ylab="log albedo", main = '90% Prediction Intervals for GOES 135 (Modified GP)')
points(ind, log(test[1:25,3]), pch=19)


#---------------------------- model check plot ---------------------------
par(mfrow = c(2,3))

# GOES 075
quilt.plot(x=test[,1], y=test[,2],
           z=test[,3], main = 'Observed Albedo')
map("world", add = TRUE)

quilt.plot(x=test[,1], y=test[,2],
           z=exp(pred.pp1), main = 'Predicted Albedo')
map("world", add = TRUE)

quilt.plot(x=test[,1], y=test[,2],
           z=exp(pred.pp1.mod), main = 'Predicted Albedo (Modified GP)')
map("world", add = TRUE)

# GOES 135
quilt.plot(x=test[,1], y=test[,2],
           z=test[,3], main = 'Observed Albedo')
map("world", add = TRUE)

quilt.plot(x=test[,1], y=test[,2],
           z=exp(pred.pp2), main = 'Predicted Albedo')
map("world", add = TRUE)

quilt.plot(x=test[,1], y=test[,2],
           z=exp(pred.pp2.mod), main = 'Predicted Albedo (Modified GP)')
map("world", add = TRUE)

# -------------------------- predict over whole area ---------------------

grid <- as.matrix(expand.grid(seq(-126, -90, l=100), seq(10.5, 34, l=100)))

data(wrld_simpl)
pts <- SpatialPoints(grid, proj4string=CRS(proj4string(wrld_simpl)))
pred.grid <- pts[!is.na(over(pts, wrld_simpl)$FIPS)]@coords

# make sure the grid covers the whole area
plot(pred.grid[,1], pred.grid[,2], pch=16)
map("world", add = TRUE)

X.grid = cbind(1, pred.grid[,1], pred.grid[,2], I(pred.grid[,1]^2), I(pred.grid[,2]^2), pred.grid[,1]*pred.grid[,2])

# GOES 075
pp1.grid = spPredict(pp1, pred.coords = pred.grid[,1:2], pred.covars = X.grid)
pred.pp1.grid = apply(pp1.grid$p.y.predictive.samples, 1, mean)

# modified
pp1.grid.mod = spPredict(pp1.mod, pred.coords = pred.grid[,1:2], pred.covars = X.grid)
pred.pp1.grid.mod  = apply(pp1.grid$p.y.predictive.samples, 1, mean)

# GOES 135
pp2.grid = spPredict(pp2, pred.coords = pred.grid[,1:2], pred.covars = X.grid)
pred.pp2.grid = apply(pp2.grid$p.y.predictive.samples, 1, mean)

# modified
pp2.grid.mod = spPredict(pp2.mod, pred.coords = pred.grid[,1:2], pred.covars = X.grid)
pred.pp2.grid.mod  = apply(pp2.grid$p.y.predictive.samples, 1, mean)


# plot predictions
par(mfrow = c(2,2))
quilt.plot(x=pred.grid[,1], y=pred.grid[,2],
           z=exp(pred.pp1.grid), main = 'GOES 075 Predicted Albedo')
map("world", add = TRUE)

quilt.plot(x=pred.grid[,1], y=pred.grid[,2],
           z=exp(pred.pp1.grid.mod), main = 'GOES 075 Predicted Albedo (modified pp)')
map("world", add = TRUE)

quilt.plot(x=pred.grid[,1], y=pred.grid[,2],
           z=exp(pred.pp2.grid), main = 'GOES 135 Predicted Albedo')
map("world", add = TRUE)

quilt.plot(x=pred.grid[,1], y=pred.grid[,2],
           z=exp(pred.pp2.grid.mod), main = 'GOES 135 Predicted Albedo (modified pp)')
map("world", add = TRUE)


#----------------------- mean difference b/w satellites ------------------------

# Using the model that provides the best fit to the data, obtain 100 samples from the predictive distribution for each satellite and use them to compute the mean difference. Perform a descriptive analysis of the mean difference. 

predictions = cbind(pred.grid, pred.pp1.grid, pred.pp2.grid)

pred.samp = data.frame(predictions[sample(nrow(predictions), size = 100),])
test.diff = data.frame(predictions[sample(nrow(predictions), size = 100),])

pred.samp$diff = pred.samp[,3] - pred.samp[,4]
test.diff$diff = test.diff[,3] - test.diff[,4]

plot(density(pred.samp$diff), xlab = 'GOES 075-GOES 135', main = 'Difference')
abline(v = mean(diff), col = 'red', lty = 2)
text(1.7, 1, labels = 'mean difference = ')
text(2.7, 1, labels = round(mean(pred.samp$diff), 4))


# fit a gaussian process to the differences
X.diff = cbind(pred.samp[,1], pred.samp[,2], pred.samp[,1]*pred.samp[,2])

loc.diff = pred.samp[,1:2]
y.diff = pred.samp$diff
geodata.diff = list(y.diff, loc.diff)
names(loc.diff) = c('lon', 'lat')
names(geodata.diff) = c('data', 'coords')

mod.diff = step(lm(y.diff ~ . + .^2 + I(loc.diff^2), data = data.frame(loc.diff)), scope = list("lower" = lm(y.diff ~ 1), "upper" = lm(y.diff ~ . + .^2 + I(loc.diff^2), data = data.frame(loc.diff))), direction = "both")
summary(mod.diff)
AIC(mod.diff)

# covariates: lon, lat, lon:lat

# find starting values
mle.diff = likfit(geodata.diff, trend = '1st', cov.model = 'matern', kappa = 2.5, ini.cov.pars = c(0.2,1))
summary(mle.diff)
# 2.5: sigsq = 0.3808, phi = 2.089, nugget = 0.0016, AIC = -109.4


priors = list("beta.flat", "sigma.sq.ig" = c(1, 0.5), "tau.sq.ig" = c(1, 0.01),
              "phi.Unif" = c(1e-6, 10), "nu.Unif" = c(0.1, 3.5))
starting = list("sigma.sq" = 0.3808, "tau.sq" = 0.0016,
                "phi" = 1/2.089, "nu" = 2.5)
tuning = list("sigma.sq" = 0.01, "tau.sq" = 0.01, "phi" = 0.01, "nu" = 0)


pp.diff = spLM(geodata.diff$data ~ X.diff, coords = as.matrix(geodata.diff$coords), knots = knots, n.samples = 5000, cov.model = 'matern', modified.pp = FALSE, priors = priors, starting = starting, tuning = tuning)

traceplot(pp.diff$p.theta.samples)


# --------------------------------- test model --------------------------------
X.test.diff = cbind(1, test.diff[,1], test.diff[,2], test.diff[,1]*test.diff[,2])

test.diff.predict = spPredict(pp.diff, pred.coords = test.diff[,1:2], pred.covars = X.test.diff)


test.pred.means  = apply(test.diff.predict$p.y.predictive.samples, 1, mean)

# plot predictions
par(mfrow = c(1,2))

quilt.plot(x=test.diff[,1], y=test.diff[,2],
           z=test.diff$diff, main = 'Differences')
map("world", add = TRUE)

quilt.plot(x=test.diff[,1], y=test.diff[,2],
           z=test.pred.means, main = 'Predicted Differences')
map("world", add = TRUE)

# looks good!

pred.int.test = apply(test.diff.predict$p.y.predictive.samples, 1, function(x) quantile(x, c(0.05, 0.95)))

1-(sum(test.diff$diff > pred.int.test[2,] | test.diff$diff < pred.int.test[1,])/nrow(test.diff))

# 97% accuracy

plot(density(test.diff.predict$p.y.predictive.samples[1,]))
abline(v = test.diff$diff[1])

# ----------------------- predict over whole area --------------------

X.diff.grid = cbind(1, pred.grid[,1], pred.grid[,2], pred.grid[,1]*pred.grid[,2])

test.diff.predict = spPredict(pp.diff, pred.coords = pred.grid[,1:2], pred.covars = X.diff.grid)

diff = apply(test.diff.predict$p.y.predictive.samples, 1, mean)

quilt.plot(x=pred.grid[,1], y=pred.grid[,2],
           z=diff, main = 'Predicted Differences')
map("world", add = TRUE)

# just use longitude
