#GLM homework 1
#problem 3b
y = c(-0.774, 0.597, 7.575, 0.397, -0.865, -0.318, -0.125, 0.961, 1.039)

theta0 <- 0

t.prime <- function(y, theta0) {
  return(sum((2*y-2*theta0)/((1+(y-theta0)^2))))
}

t.dprime <- function(y, theta0) {
  numerator <- -1*(1+(y-theta0)^2)+(2*(y-theta0)^2)
  denom <- (1+(y-theta0)^2)^2
  return(2*sum(numerator/denom))
}
  
theta1 <- theta0 - t.prime(y, theta0)/t.dprime(y, theta0)

#Newton-Raphson
k <- 50
theta <- NA
theta[1] <- theta0
for (i in 1:k)
  {theta[i+1] <- theta[i] - t.prime(y, theta[i])/t.dprime(y, theta[i])}
#theta converges to 0.1791862 at the 3rd iteration

#Fisher's Scoring
fish.info <- length(y)/2
k <-  50
theta <- NA
theta[1] <- theta0
for (i in 1:k)
  {theta[i+1] <- theta[i] + t.prime(y, theta[i])/fish.info}
#theta converges to 0.1791862 at the 14th iteration


likelihood <- function(y, theta){
  n <- length(y)
  l <- 1
  for (i in 1:length(y))
    l <- l * (pi*(1+(y[i]-theta)^2))
  return(1/l)
}

#plot the likelihood
theta.grid = seq(-3, 3, length = 500)
plot(theta.grid, likelihood(y, theta.grid), type = 'l', main = 'Likelihood Plot')

abline(v=theta[10], col = 'red')

#Problem 3c
y <- c(0, 5, 9)

theta0 = -1
#Newton-Raphson: does not converge
#Fisher's Scoring: converges to 0.36233636 at the 11th iteration

theta0 = 4.67
#Newton-Raphson: converges to 5.047230 at the 4th iteration
#Fisher's Scoring: converges to 5.047230 at the 9th iteration

theta0 = 10
#Newton-Raphson: does not converge
#Fisher's Scoring: converges to 8.545466 at the 12th iteration


#Problem 4
#Each row represents a year from 1984-1988
y = c(1,6,16,23,27,39,31,30,43,51,63,70,88,97,91,104,110,113,149,159)

x <- log(1:20)
  
j <- function(x, beta) {
  J11 <- sum(exp(beta[1] + beta[2]*x))
  J12 <- sum(x*exp(beta[1] + beta[2] * x))
  J22 <- sum(x^2*exp(beta[1] + beta[2]*x))
  j <- matrix(c(J11, J12, J12, J22), nrow = 2)
  return(solve(j))
}

U <- function(x, y, beta){
  U1 = sum(y - exp(beta[1] + beta[2]*x))
  U2 = sum(x * (y - exp(beta[1] + beta[2]*x) ))
  return(matrix(c(U1, U2), nrow = 2, ncol = 1))
}

beta <- c(1, 1)

update <- function(x, y, beta){
  beta_nplus1 <- beta + (j(x, beta) %*% U(x, y, beta))
  return(beta_nplus1)
}



iter.1 <- update(x, y, beta)
iter.2 <- update(x, y, iter.1)
iter.3 <- update(x, y, iter.2)
iter.4 <- update(x, y, iter.3)
iter.5 <- update(x, y, iter.4)
iter.6 <- update(x, y, iter.5)

iterations <- cbind(iter.1, iter.2, iter.3, iter.4, iter.5, iter.6)
#converges to beta1=0.995998, beta2=1.326610 on the 5th iteration


#4b
model = glm(y ~ x, family = 'poisson')
summary(model)

#results confirmed with glm function
