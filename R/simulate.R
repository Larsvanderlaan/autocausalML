

generate_data_linear <- function(n, const = 1.5) {

  W1 <- runif(n, -1 , 1)
  W2 <- runif(n, -1 , 1)
  W3 <- runif(n, -1 , 1)
  W4 <- runif(n, -1 , 1)
  X <- cbind(W1, W2, W3, W4)
  colnames(X) <- c("W1", "W2", "W3", "W4")
  pi <- plogis(const*(W1+W2 + W3 + W4))
  A <- rbinom(n, 1, pi)
  g1 <- 1 + 1*(1 + W1 + W2 + W3 + W4 ) + W1 + W2 + W3 + W4
  g0 <- 1 + W1 + W2 + W3 + W4
  Y <- rnorm(n, A*g1 + (1-A)*g0, 1.5)
  return(list(X=X, A=A, Y=Y, data= as.data.frame(cbind(X,A,Y)), g1  = g1, g0 = g0, pi=pi,  ATE =  1))
}

generate_data_nonlinear <- function(n, const = 1.5) {

  W1 <- runif(n, -1 , 1)
  W2 <- runif(n, -1 , 1)
  W3 <- runif(n, -1 , 1)
  W4 <- runif(n, -1 , 1)
  X <- cbind(W1, W2, W3, W4)
  colnames(X) <- c("W1", "W2", "W3", "W4")
  A <- rbinom(n, 1, plogis(const*(W1+W2 + W3 + W4)))
  Y <- rnorm(n, 1 + A + A*(abs(W1) + sin(4*W2) + abs(W3) + cos(4*W4)  ) + sin(4*W1) + abs(W2) + sin(4*W3) + exp(W4) , 1.5)
  return(list(X=X, A=A, Y=Y, data= as.data.frame(cbind(X,A,Y))))
}

get_data_generator_linear <- function(const) {
  out <- function(n) {
    generate_data_linear(n, const)
  }
  out
}

get_data_generator_nonlinear <- function(const) {
  out <- function(n) {
    generate_data_nonlinear(n, const)
  }
  out
}


bootstrap_datasets <- function(data, nboot = 100) {
  n <- nrow(data[[1]])
  lapply(1:nboot, function(i) {
    indices <- sample(1:n, n, replace = T)
    data_boot <- list(X = data$X[indices, ], A = data$A[indices], Y = data$Y[indices])
    return(data_boot)
  })
}

get_HAL_fit <- function(data, smoothness_orders = smoothness_orders, max_degree = max_degree, num_knots = num_knots, ...) {
  A <- data$A
  hal_fit_A1 <- fit_hal(data$X[A==1,], data$Y[A==1], smoothness_orders = smoothness_orders,  num_knots = num_knots, max_degree = max_degree, family = "gaussian", ...)
  g1 <- predict(hal_fit_A1, new_data = data$X)
  hal_fit_A1$basis_list <- c(list(list(cols = 1, cutoffs = min(data$X)-1, orders = 0)),  hal_fit_A1$basis_list)
  reduced_blist_A1 <- hal_fit_A1$basis_list[hal_fit_A1$coefs!=0]

  hal_fit_A0 <- fit_hal(data$X[A==0,], data$Y[A==0], smoothness_orders = smoothness_orders,  num_knots = num_knots, max_degree = max_degree, family = "gaussian", ...)
  g0 <- predict(hal_fit_A0, new_data = data$X)
  hal_fit_A0$basis_list <- c(list(list(cols = 1, cutoffs = min(data$X)-1, orders = 0)),  hal_fit_A0$basis_list)
  reduced_blist_A0 <- hal_fit_A0$basis_list[hal_fit_A0$coefs!=0]

  sieve_list <- list(sieveA0 = reduced_blist_A0,
                     sieveA1 = reduced_blist_A1, g1 = g1, g0 = g0)
}





compute_oracle_gapprox<- function(sieve_list,  n= 25000, data_generator){
  data <- data_generator(n)
  sieveA0 <- sieve_list$sieveA0
  sieveA1 <- sieve_list$sieveA1
  xbasis_A0 <- as.matrix(hal9001::make_design_matrix(data$X, sieveA0))
  xbasis_A1 <- as.matrix(hal9001::make_design_matrix(data$X, sieveA1))
  Y <- data$Y
  A <- data$A
  V <- xbasis_A0[A==0,]

  g0 <- xbasis_A0 %*% solve(t(V)%*%V, t(V) %*% Y[A==0])
  V <- xbasis_A1[A==1,]
  g1 <- xbasis_A1 %*% solve(t(V)%*%V, t(V) %*% Y[A==1])
  return(list(g1 = g1, g0 = g0))
}




do_one_experiment <- function(n = 200, data_generator,  smoothness_orders, num_knots, max_degree, lrnr_pi, ...) {
  data <- data_generator(n = n)
  A <- data$A
  Y <- data$Y
  sieve_list <- get_HAL_fit_linear(data)
  sieve_list_sample_split <-  get_HAL_fit_linear(data_generator(n = n), smoothness_orders = smoothness_orders, max_degree = max_degree, num_knots = num_knots, ... )
  alpha_n_list <- compute_alpha_n(data, sieve_list_sample_split)
  boot_data <- bootstrap_datasets(data, 250)
  estimate <- compute_Reisz_estimate(data, alpha_n_list)
  g1 <- sieve_list$g1
  g0 <- sieve_list$g0
  pi <- compute_pi(data$data, lrnr_pi = lrnr_pi)
  g <- ifelse(data$A==1, g1, g0)
  IF <- alpha_n_list$alpha * (Y - g)  + g1 - g0 - estimate
  se <- sd(IF) /sqrt(n)
  CI_IF <- estimate + 1.96*se*c( -1, 1)
  boot_estimates <- na.omit(sapply(boot_data, function(data) {
    tryCatch({
      alpha_n_list <- compute_alpha_n(data, sieve_list_sample_split)
      compute_Reisz_estimate(data, alpha_n_list)
    },
    error = function(cond) {
      return(NA)
    })
  }))
  CI_boot <- estimate + 1.96*sd(boot_estimates)*c( -1, 1)
  se_boot <- sd(boot_estimates)
  approx_g <- compute_oracle_gapprox(sieve_list_sample_split, data_generator = data_generator)
  oracle_est <- mean(as.vector(approx_g$g1 - approx_g$g0))

  est_TMLE <- compute_TMLE(data$data, pi, g1, g0)
  est_AIPW <- compute_AIPW(data$data, pi, g1, g0)
  est_CRWO <- list(estimate = estimate, se_IF = se, se_boot = se_boot, CI_IF = CI_IF, CI_boot = CI_boot)
  output <- list(est_CRWO = est_CRWO, est_TMLE = est_TMLE, est_AIPW = est_AIPW , true_approx = oracle_est, true = ATE_linear)
  return(output)
}
