#Spatial Statistics HW1

#3 Plot all the covariograms and variograms in the tables of the second set of slides. Take the variance to be 1, and take the range parameter to be such that the correlation is .05 at a distance of one unit

#list of covariance functions with parameters variance, distance, range
cov_fn = list(spherical = function(sig2, tau, phi, nu = 0){
              ifelse(phi < tau, 0, sig2*(1 - 1.5*(tau/phi) + 0.5*(tau/phi)^3))
            },
            
            powered.exp = function(sig2, tau, phi, nu){
              stopifnot(nu > 0 && nu < 2)
              sig2 * exp(-abs(tau/phi)^nu)
            },
            
            rational.quad = function(sig2, tau, phi, nu = 0){
              sig2*(1 - tau^2/(phi^2 + tau^2))
            },
            
            wave = function(sig2, tau, phi, nu = 0){
              sig2*sin(tau/phi)/(tau/phi)
            },
            
            matern = function(sig2, tau, phi, nu){
              sig2/(2^(nu-1)*gamma(nu)) * (tau/phi)^nu * besselK(tau/phi, nu)
            })
semi.variogram = function(cov_fn, d.0){
  function(sig2, tau, phi, nu){
    cov_fn(sig2, d.0, phi, nu) - cov_fn(sig2, tau, phi, nu)
  }
}

types = c("Spherical", 'Powered Exponential', 'Rational Quadratic', 'Wave', 'Matern')

#find the phi so that the correlations is r at distance tau, given nu
find_phi = function(cov_fn, r, tau, nu){
  phi.solve = function(phi)
    cov_fn(sig2 = 1, tau, phi, nu) - r
  uniroot(phi.solve, c(1E-10, 10))
}

#phi for which the correlation is 0.05, variance is 1, distance is 1, nu is 1
phi <- lapply(cov_fn, function(x) find_phi(x, r = 0.05, tau = 1, nu = 1)$root)

#plot covariance
par(mfrow=c(3,2), oma = c(0, 0, 2, 0))
for (i in 1:length(phi)) {
  plot(0, bty = 'n', type = 'n', ylim = c(-0.2, 1), xlim = c(0, 5), main = types[i], xlab = 'distance', ylab = 'covariance')
  plot.fn = function(x) cov_fn[[i]](x, sig2 = 1, phi[[i]], nu = 1)
  curve(plot.fn, add = T, lwd = 2)
  abline(h=.05, v=1, lty=2, col='grey')
}
mtext('Covariance with respect to Distance', outer = TRUE)

#plot semi-variogram
par(mfrow=c(3,2), oma = c(0, 0, 2, 0))
for (i in 1:length(phi)) {
  plot(0, bty = 'n', type = 'n', ylim = c(-0, 1.2), xlim = c(0, 5), main = types[i], xlab = 'distance', ylab = 'semi-variogram')
  plot.fn = function(x) semi.variogram(cov_fn[[i]], d.0 = 1E-10)(x, sig2 = 1, phi[[i]], nu = 1)
  curve(plot.fn, add = T, lwd = 2)
  abline(h=1, v=1, lty=2, col='grey')
}
mtext('Semi-variogram with respect to Distance', outer = TRUE)
dev.off()

#4 Assume that the correlation functions in the previous point correspond to one dimensional Gaussian processes. Simulate one 100-points realization of the process corresponding to each of the plotted functions.

library(geoR)
type = c("spherical", "powered.exponential", "cauchy", "wave", "matern")
leg.loc = c('bottomright', 'topright', 'bottomright', 'bottomleft', 'bottomleft')

#length of Gaussian process vector
n = 100
x <- seq(-5,5,len=n)
sig2 = 1
par(mfrow=c(3,2), oma = c(0, 0, 2, 0))
for (i in 1:length(type)) {
  y = grf(grid = cbind(x, 0), cov.model = type[i], cov.pars = c(sig2, phi[[i]]))
  plot(x, y$data, type = 'l', main = types[i], lwd = 2, ylab = 'process realization', xlab = 'distance')
}
mtext('Gaussian Processes', outer = TRUE)

