# Sherical Kernel
kern_sph <- function(d, v) {
  abs_d <- abs(d/phi)
  ifelse(abs_d < 1, (1 - abs_d^2) ^ v, 0)
}

# Convolution Process with Spherical Kernel
gp_conv_fit <- function(y, X, u, s,
                        B, burn, print_freq=0,
                        init=list(NULL), 
                        a_sig=2, b_sig=1,
                        a_tau=2, b_tau=1,
                        a_v=2, b_v=1, cs_v=1) {
  
  stopifnot(ncol(u) == ncol(s))
  m <- nrow(u)
  n <- nrow(s)
  p <- ncol(X)
  su <- rbind(s,u)
  D_all <- as.matrix(dist(su))
  D_su <- D_all[1:n, -c(1:n)]
  XXi <- solve(t(X) %*% X)
  
  I_n <- diag(n)
  I_m <- diag(m)
  
  ### TODO:
  update <- function(param) {
    ## update beta (conj)
    K <- kern_sph(D_su, param$v)
    z_beta <- y - K %*% param$w
    m_beta <- XXi %*% t(X) %*% z_beta
    new_beta <- mvrnorm(m_beta, XXi * param$tau2)
    
    ## update w    (conj)
    V_w <- solve(t(K) %*% K / param$tau2 + I_m / param$sig2)
    m_w <- V_w %*% t(K) %*% (y - X %*% new_beta) / param$tau2
    new_w <- mvrnorm(m_w, V_w)
    
    ## update sig2 (conj)
    new_sig2 <- 1 / rgamma(1, m/2 + a_sig, sum(new_w^2)/2 + b_sig)
    
    ## update tau2 (conj)
    resid_tau <- y - X %*% new_beta - K %*% new_w
    new_tau2 <- 1 / rgamma(1, n/2 + a_tau, sum(resid_tau^2)/2+ b_tau)
    
    ## update v    (metropolis)
    ll_v <- function(log_v) {
      v <- exp(log_v)
      sum(dnorm(y, X%*%new_beta + kern_sph(D_su,v)%*%new_w, sqrt(new_tau2), log=TRUE))
    }
    lp_v <- function(log_v) lp_log_invgamma(log_v, a_v, b_v)
    new_v <- exp(mh(log(param$v), ll_v, lp_v, cs_v))
    
    new_state <- list(beta=new_beta, w=new_w, v=new_v,
                      sig2=new_sig2, tau2=new_tau2)
    new_state
  }
  
  INIT <- if (is.null(init[[1]])) {
    list(beta=rep(0,p),
         w=rep(0,m),
         v=b_v / (a_v-1),
         sig2=b_sig / (a_sig-1),
         tau2=b_tau / (a_tau-1))
  } else init
  
  gibbs(INIT, update, B, burn, print_freq)
}


gp_conv_pred <- function(y, X, s, X_new, s_new, u, post) {
  stopifnot(nrow(X) == nrow(s)) 
  stopifnot(nrow(X_new) == nrow(s_new)) 
  
  n <- nrow(X)
  m <- nrow(u)
  n_new <- nrow(X_new)
  #D_all <- as.matrix(dist(rbind(s_new, s, u)))
  #D_sns_u <- D_all[1:(n+n_new), -c(1:(n+n_new))]
  D_sn_u <- as.matrix(dist(rbind(s_new, u)))[1:n_new, -c(1:n_new)]
  
  one_pred <- function(p) {
    rnorm(n_new, X_new %*% p$b + kern_sph(D_sn_u,p$v) %*% p$w, sqrt(p$tau2))
  }
  
  sapply(post, one_pred)
}