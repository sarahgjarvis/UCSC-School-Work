loglikesillrange = function(s, y, D, sig2, phi, kappa, tau2){
  n = length(y)
  k = ncol(D)
  I_n = diag(n)
  V = geoR::matern(as.matrix(dist(s)), phi=phi, kappa=kappa) + tau2 * I_n
  Vinv = solve(V)
  DVinvD = t(D) %*% Vinv %*% D
  beta.hat = solve(DVinvD) %*% (t(D) %*% Vinv %*% y)
  resid = y - D %*% beta.hat
  S2 = t(resid) %*% Vinv %*% resid
  
  ldV = determinant(V, logarithm = T)$modulus
  ldDVinvD = determinant(DVinvD, logarithm = T)$modulus
  
  return((-n/2) * log(sig2) - (1/2) * ldV - S2/(2*sig2) - (k/2) * log(sig2) - (1/2)*ldDVinvD)
}


marglikerange <- function(s, y, D, phi, kappa, tau2OverSig2){
  n = length(y)
  k = ncol(D)
  I_n = diag(n)
  V = geoR::matern(as.matrix(dist(s)), phi=phi, kappa=kappa) + tau2OverSig2 * I_n
  Vinv = solve(V)
  DVinvD = t(D) %*% Vinv %*% D
  beta.hat = solve(DVinvD) %*% (t(D) %*% Vinv %*% y)
  resid = y - D %*% beta.hat
  S2 = t(resid) %*% Vinv %*% resid
  
  ldV = determinant(V, logarithm = T)$modulus
  ldDVinvD = determinant(DVinvD, logarithm = T)$modulus
  
  return( - (1/2) * ldV - (n-k)/2*log(S2) - (1/2)*ldDVinvD)
}

