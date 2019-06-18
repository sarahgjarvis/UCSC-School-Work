tausq <- sigsq <- rep(NA, 11000)
tausq[1] <- sigsq[1] <- 1

sigisq <- mui <- matrix(0, 11000, 86)
sigisq[1,] <- mui[1,] <- 1

beta <- matrix(0, 11000, 6)
beta[1,] <- 1

X = model.matrix(temp ~ lat + lon +Type, sst)
n <- nrow(sst)
n.samples <- 2000
yi <- sst$temp
ni <- sst$N
alpha <- 1
a <- 1
b <- 1
for (i in 2:11000) 
  {
  v <- (ni/sigi[i-1,])+(1/tau[i-1])^(-1)
  w <- (yi/(sigi[i-1,]/ni))+((t(X)*beta[i-1,])/tau[i-1])
  mui[i,] <- rnorm(86, v*w, v)
  
  tausq[i] <- 1/rgamma(1, 42, 0.5*sum((mui[i,]-t(X)*beta[i-1,])^2))
  
  beta[i,] <- mvrnorm(1, mat %*% mui[i,], tausq[i] * xtxinv)
  
  sigsq[i] <- rgamma(1, 86*alpha+86+a, b+alpha*sum(1/sigisq[i-1,]))
  
  sigisq[i,] <- 1/rgamma(1, alpha+(3/2), alpha*sigsq[i] + 0.5*ni*(yi-mui[i,])^2)
}


dim(xtxinv %*% t(X) %*% mui[i,])
dim(xtxinv %*% t(X))
dim(mui[1,])
rep(1, 86)
lon


as.matrix(X, nrow=86, ncol=3)
