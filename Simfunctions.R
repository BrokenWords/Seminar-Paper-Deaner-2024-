


# treatment model
mu0 <- function(alpha) {
  return(alpha + alpha^2)
}
mu1 <- function(alpha) {
  return(2*alpha + alpha^2 + 1)
}
px <- function(alpha) {
  alpha <- alpha - 0.5
  return(exp(alpha+alpha^2)/(1 + exp(alpha+alpha^2)))
}


# 1 is not used anymore.
hdfun1 <- function(u,v) {
  out <- sin(pi*(u+v))
  return(out)
}

hdfun2 <- function(u,v) {
  out <- (u-v)^2
  return(out)
}

hdfun3 <- function(u, v, sig=.1) {
  out <- (1/sqrt(2*pi)/sig) * exp(-(u-v)^2/sig^2)
  return(out)
}

hdfun4 <- function(u, v, sig=.1) {
  out <- exp(-abs(u-v)/sig)
  return(out)
}

funlist <- list(hdfun1=hdfun1, hdfun2=hdfun2, hdfun3=hdfun3, hdfun4=hdfun4)


# DGP generator
# sigma: s.d. for error in pretreatment model
dgp <- function(n, T0, hdmodel, sigma=0.5) {
  # latent variable
  alpha <- runif(n, 0, 1)
  
  # outcome
  Y0 <- mu0(alpha) + rnorm(n)
  Y1 <- mu1(alpha) + rnorm(n)
  
  # logistic assignment
  prob  <- px(alpha)
  runis <- runif(n,0,1)
  D     <- ifelse(runis < prob,1,0)
  
  y <- D*Y1+(1-D)*Y0
  
  # HD covariates
  lambda <- runif(T0, 0, 1)
  pre <- outer(alpha, lambda, FUN = hdmodel) + matrix(rnorm(n*T0, 0, sigma), n, T0)   # n by p matrix
  
  return(list(y=y, y0=Y0, y1=Y1, pre=pre, alpha=alpha, lambda=lambda, d=D,prob=prob))
}



partition_folds <- function(data, K) {
  #to get loo, set K=n
  n <- length(data$y)
  folds <- sample(rep(1:K, length.out=n))
  
  # add上 fold assignment
  data$folds <- folds
  
  return(data)
}



# compute pseudo distance matrix d_hat (n x n)
compute_distance_matrix <- function(data) {
  pre   <- data$pre   # n x T0 matrix
  folds <- data$folds #  the corresponding fold
  n  <- nrow(pre)
  T0 <- ncol(pre)
  
  d_hat <- matrix(NA, n, n)
  
  for (i in 1:n) {
    for (j in 1:n) {
      # exclude this same individual or in the same fold
      if (i == j || folds[i] == folds[j]) next 
      
      dif <- pre[i, ] - pre[j, ]
      inner_prods <- pre[-c(i, j), , drop=FALSE] %*% dif
      d_hat[i, j] <- max(abs(inner_prods)) / T0
    }
  }
  return(d_hat)
}

make_hlist <- function(d_hat, m=10, lower=0.1, upper=0.9) {
  # exlude NA
  vals <- as.vector(d_hat)
  vals <- vals[!is.na(vals)]
  
  # quantile range
  h_min <- quantile(vals, lower, na.rm=TRUE)
  h_max <- quantile(vals, upper, na.rm=TRUE)
  
  # get list
  h.list <- seq(h_min, h_max, length.out=m)
  return(h.list)
}

# Kernel Function (Epanechnikov)
Kfun_Epanechnikov <- function(u) ifelse(abs(u) <= 1, 0.75*(1-u^2), 0)

# giving h，calculate mu0, mu1, p
estimate_mu_p <- function(data, d_hat, h, Kfun=Kfun_Epanechnikov) {
  n <- length(data$y)
  mu0_hat <- rep(NA, n)
  mu1_hat <- rep(NA, n)
  p_hat   <- rep(NA, n)
  
  for (i in 1:n) {
    train_idx <- which(data$folds != data$folds[i])   # train sample of individual i
    dists     <- d_hat[i, train_idx]
    weights   <- Kfun(dists / h)
    
    if (sum(weights) == 0) next
    weights <- weights / sum(weights)
    
    # treatment / control partition
    treated_idx <- train_idx[data$d[train_idx] == 1]
    control_idx <- train_idx[data$d[train_idx] == 0]
    
    if (length(treated_idx) > 0) {
      mu1_hat[i] <- sum(weights[data$d[train_idx] == 1] * data$y[treated_idx]) /
        sum(weights[data$d[train_idx] == 1])
    }
    if (length(control_idx) > 0) {
      mu0_hat[i] <- sum(weights[data$d[train_idx] == 0] * data$y[control_idx]) /
        sum(weights[data$d[train_idx] == 0])
    }
    
    p_hat[i] <- sum(weights * data$d[train_idx])
  }
  
  return(list(mu0=mu0_hat, mu1=mu1_hat, p=p_hat))
}

# calculate cv error
cv_error <- function(data, d_hat, h) {
  est <- estimate_mu_p(data, d_hat, h)
  mu0 <- est$mu0
  mu1 <- est$mu1
  
  treated_loss <- (data$y[data$d==1] - mu1[data$d==1])^2
  control_loss <- (data$y[data$d==0] - mu0[data$d==0])^2
  
  return(mean(c(treated_loss, control_loss), na.rm=TRUE))
}

# choose optimal bandwith h
select_bandwidth <- function(data, d_hat, hlist) {
  errors <- sapply(hlist, function(h) cv_error(data, d_hat, h))
  h.opt  <- hlist[which.min(errors)]
  return(list(h.opt=h.opt, errors=errors))
}


# estimate
estimate_ATT_and_CF <- function(data, h.opt, alpha=0.05) {
  n  <- length(data$y)
  N1 <- sum(data$d)
  
  # --- Step A: calculate nuisance params based on optimal bandwidth
  d_hat    <- compute_distance_matrix(data)
  nuisance <- estimate_mu_p(data, d_hat, h.opt)
  mu0_hat <- nuisance$mu0
  mu1_hat <- nuisance$mu1
  p_hat   <- nuisance$p
  
  # estimate ATT
  sum_term <- 0
  for (i in 1:n) {
    term <- data$y[i] * data$d[i] -
      (1 - data$d[i]) * data$y[i] * p_hat[i] / (1 - p_hat[i]) +
      (data$d[i] - p_hat[i]) * mu0_hat[i] / (1 - p_hat[i])
    sum_term <- sum_term + term
  }
  att_hat <- sum_term / N1
  
  # estimate variance
  sum_var <- 0
  for (i in 1:n) {
    term <- data$y[i] * data$d[i] -
      (1 - data$d[i]) * data$y[i] * p_hat[i] / (1 - p_hat[i]) +
      (data$d[i] - p_hat[i]) * mu0_hat[i] / (1 - p_hat[i]) -
      (N1/n) * att_hat
    sum_var <- sum_var + term^2
  }
  var_hat <- (n / (N1^2)) * sum_var
  
  # estimate conterfactual mean
  Ey1_treated <- mean(data$y[data$d == 1])
  mu0_cf_hat  <- Ey1_treated - att_hat
  
  # ci: conterfactual expectation
  se <- sqrt(var_hat / n)
  ci_length <- 2 * qnorm(1 - alpha/2) * se
  
  return(list(att=att_hat,
              mu0_cf=mu0_cf_hat,
              var=var_hat,
              ci_length=ci_length))
}

# The whole algorithm

deaner_algorithm1 <- function(data, h.list=NULL, K=5, alpha=0.05, m=10) {
  # Step 0: partition folds
  data <- partition_folds(data, K)
  
  # Step 1: compute pseudo distance
  d_hat <- compute_distance_matrix(data)
  
  # Step 2: make hlist if not provided
  if (is.null(h.list)) {
    h.list <- make_hlist(d_hat, m=m)
  }
  
  # Step 3: choose optimal bandwidth
  bw_res <- select_bandwidth(data, d_hat, h.list)
  h.opt  <- bw_res$h.opt
  
  # Step 4: estimate based on optimal bandwidth
  res <- estimate_ATT_and_CF(data, h.opt, alpha)
  
  return(list(result=res, h.opt=h.opt, hlist=h.list, cv.errors=bw_res$errors))
}

# simulation function

sim <- function(n, T0, hdmodel,theta0, sigma=0.5, K=5, alpha=0.05, m=10) {
  # Step 1: generate data
  data <- dgp(n, T0, hdmodel, sigma)
  
  # Step 2: generate hlist
  result <- deaner_algorithm1(data, h.list=NULL, K=K, alpha=alpha, m=m)
  # Step 3: add deviation & rejection
  mu0_cf_hat <- result$result$mu0_cf
  bias       <- mu0_cf_hat - theta0
  
  # rej
  half_len <- result$result$ci_length / 2
  ci_lower <- mu0_cf_hat - half_len
  ci_upper <- mu0_cf_hat + half_len
  rej      <- ifelse(theta0 < ci_lower | theta0 > ci_upper, 1, 0)
  
  result$result$bias <- bias
  result$result$rej  <- rej
  
  
  return(result$result)
}


