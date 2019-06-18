# Eigen approximation to exponential correlation

eigen.exp.approx = function(tau, J, L){
  #index
  j = 1:J
  
  # spectrum for exponential correlation at k
  f = function(k) { 1/(k^2 + 1) }
  
  lambda = f(j*pi/(2*L))
  
  psi = function(tau) {
    return(exp(1i*j*pi*tau / (2*L)))
  }
  res = Re(sum(lambda *  psi(tau)) / sum(lambda * psi(0)) )
  return(list(res=res, psi=psi, lambda=lambda))
}

L = 2
L = 10
actual <- cov_fn$powered.exp(tau = d, phi=1, sig2=1, nu=1)
eig.approx <- sapply(d, function(s) eigen.exp.approx(s, J=2, L=L))
eig.approx2 <- sapply(d, function(s) eigen.exp.approx(s, J=10, L=L))
eig.approx3 <- sapply(d, function(s) eigen.exp.approx(s, J=20, L=L))
eig.approx4 <- sapply(d, function(s) eigen.exp.approx(s, J=100, L=L))


plot(d, actual, type = 'l', xlab = 'distance' ~ (tau), ylab = 'correlation', main = 'Eigen Approximate of Exponential Correlation', ylim = c(-0.2, 1), lwd = 2)
lines(d, eig.approx, lwd = 3, lty = 2, col = 5)
lines(d, eig.approx2, lwd = 3, lty = 2, col = 2)
lines(d, eig.approx3, lwd = 3, lty = 2, col = 4)
lines(d, eig.approx4, lwd = 3, lty = 2, col = 3)
legend('bottomleft', legend = c('actual','J=2', 'J=10', 'J=20', 'J=100'), bty = 'n', lwd=2, col = c(1, 5, 2, 4, 3), cex=0.9)
text(3.5, 0.9, labels = 'L=10', cex = 2)

