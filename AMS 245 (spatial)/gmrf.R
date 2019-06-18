### Gibbs Sampler (generic) ###
gibbs <- function(init, update, B, burn, print_every=0) {
  out <- as.list(1:B)
  out[[1]] <- init
  
  for (i in 1:(B+burn)) {
    if (i - burn < 2) {
      out[[1]] <- update(out[[1]])
    } else {
      out[[i-burn]] <- update(out[[i-burn-1]])
    }
    if (print_every > 0 && i %% print_every == 0) {
      cat("\rProgress: ", i, " / ", B+burn)
    }
  }
  cat("\n")
  
  out
}

### genW ###
genW <- function(D, eps=1E-4) {
  # Generate unscaled precision matrix (W)
  # from distance matrix (D)
  
  #' @param eps: in computing distance between knot points, to determine if distance is minimal,
  #              a small precision error needs to be tolerated. Default is 1E-3.
  min_d <- min(c(D)[c(D) > 0])
  W <- matrix(0, nrow(D), ncol(D))
  
  for (i in 1:nrow(W)) {
    neighbors <- which(abs(D[i,] - min_d) < eps)
    # number of neighbors
    W[i,i] <- length(neighbors)
    W[i,neighbors] <- -1
  }
  
  W
}
# Test genW(D): ###
#grid <- expand.grid(1:3, 1:3)
#D <- as.matrix(dist(grid))
#W <- genW(D)
#solve(W + diag(1e-3,9))


### Convolution Process with GMRF prior for weights ###
gp_conv_gmrf_fit <- function(y, X, u, s,
                             B, burn, print_freq=0,
                             init=list(NULL), 
                             a_tau=2, b_tau=1,
                             a_lamz=.01, b_lamz=.01,
                             a_v=2, b_v=1, cs_v=1,
                             eps=1e-4) {
  
  #' @param eps: parameter in genW(D)
  
  stopifnot(ncol(u) == ncol(s))
  m <- nrow(u)
  n <- nrow(s)
  p <- ncol(X)
  su <- rbind(s,u)
  D_all <- as.matrix(dist(su))
  D_su <- D_all[1:n, -c(1:n)]
  D_u <- D_all[-c(1:n), -c(1:n)]
  W <- genW(D_u, eps)
  XXi <- solve(t(X) %*% X)
  
  I_n <- diag(n)
  I_m <- diag(m)
  
  ### TODO:
  update <- function(param) {
    ## update beta (conj)
    K <- kern_sph(D_su, param$v)
    m_beta <- XXi %*% t(X) %*% (y - K %*% param$z)
    new_beta <- mvrnorm(m_beta, XXi * param$tau2)
    
    ## update z    (conj)
    V_z <- solve((t(K) %*% K)/param$tau2 + W*param$lamz)
    m_z <- V_z %*% t(K) %*% (y - X %*% new_beta) / param$tau2
    new_z <- mvrnorm(m_z, V_z)
    
    ## update lamz (conj)
    new_lamz <- rgamma(1, a_lamz + m/2, b_lamz + t(new_z)%*% W %*% new_z/2)
    
    ## update tau2 (conj)
    resid <- y - X %*% new_beta - K %*% new_z
    new_tau2 <- 1 / rgamma(1, a_tau + n/2, b_tau + sum(resid^2)/2)
    
    ## update v    (metropolis)
    ll_v <- function(log_v) {
      v <- exp(log_v)
      sum(dnorm(y, X%*%new_beta + kern_sph(D_su,v)%*%new_z, sqrt(new_tau2), log=TRUE))
    }
    lp_v <- function(log_v) lp_log_invgamma(log_v, a_v, b_v)
    new_v <- exp(mh(log(param$v), ll_v, lp_v, cs_v))
    
    new_state <- list(beta=new_beta, z=new_z, v=new_v,
                      tau2=new_tau2, lamz=new_lamz)
    new_state
  }
  
  INIT <- if (is.null(init[[1]])) {
    list(beta=rep(0,p),
         z=rep(0,m),
         v=b_v / (a_v-1),
         tau2=b_tau / (a_tau-1),
         lamz=a_lamz / b_lamz)
  } else init
  
  gibbs(INIT, update, B, burn, print_freq)
}


gp_conv_gmrf_pred <- function(y, X, s, X_new, s_new, u, post) {
  stopifnot(nrow(X) == nrow(s)) 
  stopifnot(nrow(X_new) == nrow(s_new)) 
  
  n <- nrow(X)
  m <- nrow(u)
  n_new <- nrow(X_new)
  #D_all <- as.matrix(dist(rbind(s_new, s, u)))
  #D_sns_u <- D_all[1:(n+n_new), -c(1:(n+n_new))]
  D_sn_u <- as.matrix(dist(rbind(s_new, u)))[1:n_new, -c(1:n_new)]
  
  one_pred <- function(p) {
    rnorm(n_new, X_new %*% p$b + kern_sph(D_sn_u,p$v) %*% p$z, sqrt(p$tau2))
  }
  
  sapply(post, one_pred)
}
