data = sarah_data
library(rworldmap)
library(geoR)

color.gradient <- function(x, colors=c("blue","purple","red"), colsteps=10) {
  return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
}

newmap <- getMap(resolution = "low")

# plot of all of my data
# layout
{layout.matrix <- matrix(c(1, 2), nrow = 1, ncol = 2)

layout(mat = layout.matrix,
       widths = c(5, 2))}
plot(newmap, xlim = c(-90, -85), ylim = c(5, 35), asp = 1, main = 'Albedo Measurements')
points(sarah_data$longitude, sarah_data$latitude, col = color.gradient(sarah_data$bhriso), cex = 0.3, pch=20)
#legend
{plot(x = rep(max(data$bhriso) + 0.67, 100),
      y = seq(min(data$bhriso), max(data$bhriso), length = 100),
      pch = 15, col = color.gradient(seq(0, 1, length = 100)),
      cex = 2.5, axes=F, xlab = '', ylab = '')
  
  lines(rep(max(data$bhriso) + 0.97, 2), range(data$bhriso))
  for (i in 1:6){
    lines(max(data$bhriso) + c(0.97, 1.07), rep(min(data$bhriso)+(i-1)*diff(range(data$bhriso)/5), 2))
    #   text(max(loc[,1])+1.15, min(loc[,2])+(i-1)*diff(range(loc[,2])/5),
    #       round(quantile(y, (i-1)/5), 3), pos=4)
    text(max(data$bhriso)+1.12, min(data$bhriso)+(i-1)*diff(range(data$bhriso)/5),
         round(seq(min(data$bhriso), max(data$bhriso), length = 6)[i], 3), pos=4)
  }
}

# random sample of 500
# data = data[sample(nrow(data), size = 500),]

# plot of 500 randomly selected points
# layout
{layout.matrix <- matrix(c(1, 2), nrow = 1, ncol = 2)

layout(mat = layout.matrix,
       widths = c(5, 2))} 
plot(newmap, xlim = c(-90, -85), ylim = c(5, 35), asp = 1, main = 'Albedo Measurements')
points(data$longitude, data$latitude, col = color.gradient(data$bhriso), pch = 20)
# legend
{plot(x = rep(max(data$bhriso) + 0.67, 100),
       y = seq(min(data$bhriso), max(data$bhriso), length = 100),
       pch = 15, col = color.gradient(seq(0, 1, length = 100)),
       cex = 2.5, axes=F, xlab = '', ylab = '')
  
lines(rep(max(data$bhriso) + 0.97, 2), range(data$bhriso))
for (i in 1:6){
  lines(max(data$bhriso) + c(0.97, 1.07), rep(min(data$bhriso)+(i-1)*diff(range(data$bhriso)/5), 2))
  #   text(max(loc[,1])+1.15, min(loc[,2])+(i-1)*diff(range(loc[,2])/5),
  #       round(quantile(y, (i-1)/5), 3), pos=4)
  text(max(data$bhriso)+1.12, min(data$bhriso)+(i-1)*diff(range(data$bhriso)/5),
       round(seq(min(data$bhriso), max(data$bhriso), length = 6)[i], 3), pos=4)
}
}

hist(data$bhriso, col = 'purple', main = 'Histogram of Albedo Measurement', xlab = 'BHRiso')

#log transformation of the albedo measurement
y = log(data$bhriso)

# plot transformed variable on the map
{layout.matrix <- matrix(c(1, 2), nrow = 1, ncol = 2)
  
  layout(mat = layout.matrix,
         widths = c(5, 2))} 
plot(newmap, xlim = c(-90, -85), ylim = c(5, 35), asp = 1, main = 'Transformed Albedo Measurements')
points(data$longitude, data$latitude, col = color.gradient(y), pch = 20)
#legend
{plot(x = rep(max(data$bhriso) + 0.67, 100),
      y = seq(min(data$bhriso), max(data$bhriso), length = 100),
      pch = 15, col = color.gradient(seq(0, 1, length = 100)),
      cex = 2.5, axes=F, xlab = '', ylab = '')
  
  lines(rep(max(data$bhriso) + 0.97, 2), range(data$bhriso))
  for (i in 1:6){
    lines(max(data$bhriso) + c(0.97, 1.07), rep(min(data$bhriso)+(i-1)*diff(range(data$bhriso)/5), 2))
    #   text(max(loc[,1])+1.15, min(loc[,2])+(i-1)*diff(range(loc[,2])/5),
    #       round(quantile(y, (i-1)/5), 3), pos=4)
    text(max(data$bhriso)+1.12, min(data$bhriso)+(i-1)*diff(range(data$bhriso)/5),
         round(seq(min(y), max(y), length = 6)[i], 3), pos=4, cex = 0.8)
  }
}

# show transformation is normal
par(mfrow = c(1, 2))
hist(y, xlab = 'log(BHRiso)', main = 'Transformed Albedo Meaurement', col = 'purple')
qqnorm(y)
qqline(y, col = 'red', lwd=2)

loc = cbind(data$longitude, data$latitude)

# Is there evidence of a first or second order trend function of location?
# explore trend
# fit regression with linear and quadratic terms in each direction. Select covariates with a stepwise AIC procedure.
mod = step(lm(y ~ . + .^2 + I(loc^2), data = data.frame(loc)), scope = list("lower" = lm(y ~ 1), "upper" = lm(y ~ . + .^2 + I(loc^2), data = data.frame(loc))), direction = "both")
summary(mod)
# y ~ longitude + latitude + (longitude + latitude)^2 + I(loc^2) AIR=-712.03
yhat = predict(mod)

par(mfrow = c(2,2))
plot(data$longitude, y, pch = 20, xlab = 'longitude', ylab = 'log(albedo)', main = 'Longitude')
points(data$longitude, yhat, col = 'slateblue1', pch = 20)
plot(data$latitude, y, pch = 20, xlab = 'latitude', ylab = 'log(albedo)', main = 'Latitude')
points(data$latitude, yhat, col = 'slateblue1', pch = 20)

plot(data$longitude, mod$residuals, pch = 20, ylab = 'residuals', xlab = 'longitude')
abline(h=0, col = 'red', lwd = 2)
plot(data$latitude, mod$residuals, pch = 20, ylab = 'residuals', xlab = 'latitude')
abline(h=0, col = 'red', lwd = 2)

# detrend
resid = y-yhat
plot(resid)

# Plot the variogram 

par(mfrow =c(1,2))
plot(variog(data = resid, coords = loc, op="cloud"), col = 'slateblue1', main='Semi-variogram Cloud for the Residuals')
plot(variog(data = resid, coords = loc, uvec = 16), col = 'slateblue1', main='Binned Semi-Variogram for the Residuals')

# Explore possible anisotropies using a directional variogram.
plot(variog4(data = resid, coords = loc), main = 'Directional Variogram for the Residuals', legend = F)
legend('topleft', legend = c(expression(paste("0",degree), paste("45", degree), paste("90",degree), paste("135",degree))), lty = c(1,2,3,4), col = c(1,2,3,4), cex = 0.8)

# Use least squares to fit the covariograms in the Matèrn family with smoothness equal to .5; 1; 1.5; 2.5. Plot the results. Use the plots and the values of the LSE to select the best fit.
var = variog(data = resid, coords = loc, messages = F)

smoothness = list(0.5, 1, 1.5, 2.5)
vf = lapply(smoothness, function(x) {variofit(var, kappa = x)})

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

# best fit is the second model with kappa = 1. the estimated nugget for this model is tau hat^2 = 0.00825979 and range (phi hat) = 0.9516736

# Plot the likelihood function for the sill and the range corresponding to each of the correlations in the previous point. If a nugget is needed, you can plug an estimated value.

X <- as.matrix(cbind(1,mod$model[,-1]))

n.grid = 30
# parameter grids
sig2.grid = seq(0.01, 15, length = n.grid)
phi.grid = seq(0.01, 15, length = n.grid)

# range likelihood over all 4 smoothness parameters
vals = list()
  vals[[1]] = matrix(NA, nrow = n.grid, ncol = n.grid)
    for (i in 1:n.grid){
      for (j in 1:n.grid) {
        vals[[1]][i,j] = loglikesillrange(s = loc, y, D = X, sig2.grid[j], phi.grid[i], kappa=smoothness[[1]], tau2 = vf[[1]]$nugget)
        print(i)
        print(j)}}
  vals[[2]] = matrix(NA, nrow = n.grid, ncol = n.grid)
  for (i in 1:n.grid){
    for (j in 1:n.grid) {
      vals[[2]][i,j] = loglikesillrange(s = loc, y, D = X, sig2.grid[j], phi.grid[i], kappa=smoothness[[2]], tau2 = vf[[2]]$nugget)
      print(i)
      print(j)}}
  vals[[3]] = matrix(NA, nrow = n.grid, ncol = n.grid)
  for (i in 1:n.grid){
    for (j in 1:n.grid) {
      vals[[3]][i,j] = loglikesillrange(s = loc, y, D = X, sig2.grid[j], phi.grid[i], kappa=smoothness[[3]], tau2 = vf[[3]]$nugget)
      print(i)
      print(j)}}
  vals[[4]] = matrix(NA, nrow = n.grid, ncol = n.grid)
  for (i in 1:n.grid){
    for (j in 1:n.grid) {
      vals[[4]][i,j] = loglikesillrange(s = loc, y, D = X, sig2.grid[j], phi.grid[i], kappa=smoothness[[4]], tau2 = vf[[3]]$nugget)
      print(i)
      print(j)}}


# plot likelihood
par(mfrow = c(2,2))
contour(phi.grid, sig2.grid, vals[[1]], levels = c(170,175,180,185,190,195, 197), main = 'kappa = 0.5', xlab = expression(phi), ylab = expression(sigma^2))
contour(phi.grid, sig2.grid, vals[[2]], levels = c(165, 170, 175, 180,185, 190, 192), main = 'kappa = 1', xlab = expression(phi), ylab = expression(sigma^2))
contour(phi.grid, sig2.grid, vals[[3]], levels = c(160, 170, 175, 180, 185, 190,194), main = 'kappa = 1.5', xlab = expression(phi), ylab = expression(sigma^2))
contour(phi.grid, sig2.grid, vals[[4]], levels = c(100,140, 158, 160, 170, 175,180 ,184), main = 'kappa = 2.5', xlab = expression(phi), ylab = expression(sigma^2))


# Plot the marginal likelihood for the range parameter for each of the examples above
range = list()
range[[1]] = rep(NA, length(phi.grid))
for (i in 1:n.grid){{
    range[[1]][i] = marglikerange(s = loc, y, D = X, phi.grid[i], kappa=smoothness[[1]], vf[[1]]$nugget/vf[[1]]$cov.pars[1])
    print(i)}}
range[[2]] = rep(NA, length(phi.grid))
for (i in 1:n.grid){{
  range[[2]][i] = marglikerange(s = loc, y, D = X, phi.grid[i], kappa=smoothness[[2]], vf[[2]]$nugget/vf[[2]]$cov.pars[1])
  print(i)}}
range[[3]] = rep(NA, length(phi.grid))
for (i in 1:n.grid){{
  range[[3]][i] = marglikerange(s = loc, y, D = X, phi.grid[i], kappa=smoothness[[3]], vf[[3]]$nugget/vf[[3]]$cov.pars[1])
  print(i)}}
range[[4]] = rep(NA, length(phi.grid))
for (i in 1:n.grid){{
  range[[4]][i] = marglikerange(s = loc, y, D = X, phi.grid[i], kappa=smoothness[[4]], vf[[4]]$nugget/vf[[4]]$cov.pars[1])
  print(i)}}


par(mfrow = c(2,2))
l = spline(phi.grid, range[[1]])
plot(phi.grid, range[[1]], ylim = c(min(l$y), max(l$y)), xlab = expression(phi), ylab ='log-likelihood', main = 'kappa = 0.5')
lines(l)
abline(v=l$x[which.max(l$y)], lty = 2, col = 'red')

l = spline(phi.grid, range[[2]])
plot(phi.grid, range[[2]], ylim = c(min(l$y), max(l$y)), xlab = expression(phi), ylab ='log-likelihood', main = 'kappa = 1')
lines(l)
abline(v=l$x[which.max(l$y)], lty = 2, col = 'red')

l = spline(phi.grid, range[[3]])
plot(phi.grid, range[[3]], ylim = c(min(l$y), max(l$y)), xlab = expression(phi), ylab ='log-likelihood', main = 'kappa = 1.5')
lines(l)
abline(v=l$x[which.max(l$y)], lty = 2, col = 'red')

l = spline(phi.grid, range[[4]])
plot(phi.grid, range[[4]], ylim = c(min(l$y), max(l$y)), xlab = expression(phi), ylab ='log-likelihood', main = 'kappa = 2.5')
lines(l)
abline(v=l$x[which.max(l$y)], lty = 2, col = 'red')


