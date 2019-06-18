semiv = function(u, phi, nu, sig2, kappa)
  sig2 * (1-matern(u, phi, nu)) + kappa
make.f = function(nu){
  f.opt = function(p){
    phi = p[1]
    sig2 = p[2]
    kappa = p[3]
    sum( (semiv(var$u, phi, nu, sig2, kappa) - var$v)^2)
  }
  return (f.opt)
}


