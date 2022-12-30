


library(autocausalML)
library(glmnet)
library(hal9001)
multiplier <-1
sim_beta <- function(n, d, k) {
  pi <- 3.14
  X <- as.matrix(replicate(d, runif(n, -1 , 1)))
  colnames(X) <- paste0("X", 1:ncol(X))
  zeta <- rowSums(X[,1:k, drop = F]) / sqrt(k)
  eta <- sign(zeta) * sqrt(abs(zeta))
  alpha <- pmax(0.0001, pmin(0.9999, 1 / (1+ exp(- 3*eta))))
  quantile(alpha)
  W <- rbeta(n, alpha, 1 - alpha)
  mu <- eta + 0.2* (alpha - 0.5)
  tau <- -0.2
  theta <- mu + W * multiplier*tau  + 3
  if(any(theta <0 )) {
    print(range(theta))
    warning("Theta has negative values")
  }
  Y <- rgamma(n, shape = theta^2, rate = abs(theta))
  tau <-  abs(mu + 1 * multiplier*tau  + 3) -   abs(mu + 0 * multiplier*tau  + 3)
  data <- as.data.frame(cbind(X,W,Y, theta = theta, tau = tau, ATE = mean(tau)))

  return(data)
}

sim_pois <- function(n, d, k) {
  pi <- 3.14
  X <- as.matrix(replicate(d, runif(n, -1 , 1)))
  colnames(X) <- paste0("X", 1:ncol(X))
  tau <- rowMeans(cos(X[, 1:k, drop = F]*pi /3))
  lambda <- 0.2 + tau^2 + rowMeans(pmax(X,0))
  mu <- 4 * rowMeans(X) + 2* lambda
  W <- rpois(n, lambda)
  theta <- mu + W * multiplier*tau + 4
  if(any(theta <0 )) {
    warning("Theta has negative values")
  }
  Y <- rgamma(n, shape = theta^2, rate = abs(theta))
  tau <-  abs(mu + 1 * multiplier*tau  + 4) -   abs(mu + 0 * multiplier*tau  + 4)
  return(as.data.frame(cbind(X,W,Y, theta = theta, tau = tau, ATE = mean(tau))))
}

sim_norm <- function(n, d, k) {
  pi <- 3.14
  X <- as.matrix(replicate(d, runif(n, -1 , 1)))
  colnames(X) <- paste0("X", 1:ncol(X))
  eta <- sign(X[, 1:k, drop = F]) * sqrt(abs(X[, 1:k, drop = F]))
  mu <-2^(k-1) * apply(eta, 1 , prod)
  #mu <-   sign(eta) * sqrt(abs(eta))
  lambda <- 0.1 * sign(mu) + mu + rowMeans(X)
  tau <- (pmax(X[,1],0) + pmax(X[,2], 0)) / 2
  W <- rnorm(n, mean = lambda, sd = 0.3 + abs(lambda)/3)
  theta <- mu + W * multiplier*tau + 10
  if(any(theta <0 )) {
    warning("Theta has negative values")
  }
  Y <- rgamma(n, shape = theta^2, rate = abs(theta))
  tau <-  (mu + 1 * multiplier*tau  + 10) -   (mu + 0 * multiplier*tau  + 10)
  return(as.data.frame(cbind(X,W,Y, theta = theta, tau = tau, ATE = mean(tau))))
}

sim_log <- function(n, d, k) {
  pi <- 3.14
  X <- as.matrix(replicate(d, runif(n, -1 , 1)))
  colnames(X) <- paste0("X", 1:ncol(X))
  zeta <- rowSums(X[,1:k, drop = F]) / sqrt(k)
  mu <- pmax(0, 2*zeta) + rowMeans(X)
  lambda <- 1 / (1 + exp(-sign(zeta) * zeta^2))
  tau <- sin(2 * pi * X[,1])
  W <- exp(rnorm(n, mean = lambda, sd = 1/3))
  theta <- mu + W * multiplier*tau  + 6
  if(any(theta <0 )) {
    warning("Theta has negative values")
  }
  Y <- rgamma(n, shape = theta^2, rate = abs(theta))
  tau <-  abs(mu + 1 * multiplier*tau  + 6) -   abs(mu + 0 * multiplier*tau  + 6)
  return(as.data.frame(cbind(X,W,Y, theta = theta, tau = tau, ATE = mean(tau))))
}

sim_bin <- function(n, d, k) {
  pi <- 3.14
  X <- as.matrix(replicate(d, runif(n, -1 , 1)))
  colnames(X) <- paste0("X", 1:ncol(X))
  zeta <- 2*rowSums(X[,1:k, drop = F]) / sqrt(k)
  mu <- pmax(0, 2*zeta) + rowMeans(X)
  lambda <- 1 / (1 + exp(-sign(zeta) * zeta^2))
  tau <- X[,1] + X[,2]
  W <- rbinom(n, 1, lambda)
  theta <- mu + W * multiplier*tau  + 2
  if(any(theta <0 )) {
    warning("Theta has negative values")
  }
  Y <- rgamma(n, shape = theta^2, rate = abs(theta))
  tau <-  abs(mu + 1 * multiplier*tau  + 2) -   abs(mu + 0 * multiplier*tau  + 2)
  return(as.data.frame(cbind(X,W,Y, theta = theta, tau = tau, ATE = mean(tau))))
}





get_projection <- function(sim_generator, d, k, nknots, nknots_tau, s) {
  data_big <- sim_generator(n= 25000,d,k)
  X <- as.matrix(data_big[,paste0("X", 1:d)])
  A <- data_big$W
  Y <-  data_big$Y
  theta <- data_big$theta

  basis_list_mu <- hal9001::enumerate_basis((X), max_degree = 2, num_knots = nknots, smoothness_orders =s)
  basis_list_tau <- hal9001::enumerate_basis((X), max_degree = 1, num_knots = nknots_tau, smoothness_orders = 1)
  x_basis_mu <- cbind(1, as.matrix(make_design_matrix(X, basis_list_mu, 0.9)))
  x_basis_tau <- cbind(1, as.matrix(make_design_matrix(X, basis_list_tau, 0.9)))
  x_basis <- cbind(  x_basis_mu, A * x_basis_tau)
  penalty.factor <-  rep(1, ncol(x_basis))
  #fit_A <- fit_hal(cbind(X,A), theta,   family = "gaussian", smoothness_orders = 0, max_degree =3, num_knots=  c(50,50, 10), fit_control = list(parallel = TRUE))
  #mean(predict(fit_A, new_data = cbind(X,1)) -  predict(fit_A, new_data=cbind(X,0)))
  lasso_fit <- glmnet::glmnet(x_basis, theta, family = "gaussian", lambda = 1e-8,  intercept = F,   standardize = FALSE, lambda.min.ratio = 1e-4)
  beta <- as.vector(coef(lasso_fit, s = "lambda.min"))[-1]

 # lasso_fit <- speedglm::speedlm.fit(theta,x_basis,  intercept = F )

  # beta <- coef(lasso_fit)
  beta[is.na(beta)] <- 0
  x_basis_lasso_A1 <- cbind(x_basis_mu, 1 * x_basis_tau)
  x_basis_lasso_A0 <- cbind(x_basis_mu, 0 * x_basis_tau)

  ATE_approx <- mean((x_basis_lasso_A1 - x_basis_lasso_A0) %*% beta)
  ATE_approx
}





run_sims <- function(n, d, k, sim_generator, nsims, nknots, nknots_tau, s) {
  data_big <- sim_generator(n= 500000,d,k)

  trueATE <- mean(data_big$ATE)
  approxATE <- get_projection(sim_generator, d, k, nknots, nknots_tau, s)
  sim_list <- list()
  for(iter in 1:nsims) {
    print(iter)
    try({
      data <- sim_generator(n,d,k)
      X <- as.matrix(data[,paste0("X", 1:d)])
      A <- data$W
      Y <- data$Y

      print(trueATE)
      print(approxATE)
      #fit_hal_g_params <- list(formula <- ~ h(.) + h(.,A), max_degree = 2, num_knots = c(20,20), smoothness_orders = 1)
      basis_list_mu <- hal9001::enumerate_basis((X), max_degree = 2, num_knots = nknots, smoothness_orders =s)
      basis_list_tau <- hal9001::enumerate_basis((X), max_degree = 1, num_knots = nknots_tau, smoothness_orders = 1)

      x_basis_mu <- cbind(1,  as.matrix(make_design_matrix(X, basis_list_mu, 0.9)))
      x_basis_tau <- cbind(1, as.matrix(make_design_matrix(X, basis_list_tau, 0.9)))

      # add intercept
      x_basis <- cbind(   x_basis_mu,  A * x_basis_tau)
      penalty.factor <-  rep(1, ncol(x_basis))
      penalty.factor[1] <- 0
      penalty.factor[1 + ncol(x_basis_mu) ] <-0

      # Compute lasso fit
      lasso_fit <- glmnet::cv.glmnet(x_basis, Y, family = "gaussian", intercept = F, penalty.factor = penalty.factor, relax = F, gamma = 0, parallel = T, standardize = FALSE, lambda.min.ratio = 1e-5)
      beta <- as.vector(coef(lasso_fit, s = "lambda.min"))[-1]
      keep_cols_1 <- which(abs(beta) > 1e-8)

      x_basis_lasso <- x_basis[, keep_cols_1, drop = F]


      # Compute relaxed fit
      #relaxed_fit <- glmnet::cv.glmnet( x_basis_lasso, Y, family = "gaussian" ,intercept = F) # screen hyper correlates variables just in case
      relaxed_fit <- glm.fit( x_basis_lasso, Y, family = gaussian() ,intercept = F)
      beta_relaxed <- coef(relaxed_fit)#[-1]
      keep_cols_2 <-  which(!is.na(beta_relaxed) & abs(beta_relaxed) > 1e-8 )
      x_basis_lasso <- x_basis_lasso[,keep_cols_2, drop = F]
      x_basis_lasso_A1 <- cbind(x_basis_mu, 1 * x_basis_tau)[, keep_cols_1, drop = F][, keep_cols_2, drop = F]
      x_basis_lasso_A0 <- cbind(x_basis_mu, 0 * x_basis_tau)[, keep_cols_1, drop = F][, keep_cols_2, drop = F]
      relaxed_fit <- glm.fit(x_basis_lasso, Y, family = gaussian(), intercept = F) # redo fit with uncorrelated variables
      beta <- coef(relaxed_fit) # get final beta estimate for AdaLSE

      theta_n <- x_basis_lasso %*% beta
      theta_n_A1 <- x_basis_lasso_A1 %*% beta
      theta_n_A0 <- x_basis_lasso_A0 %*% beta


      solve(t(x_basis_lasso)  %*% x_basis_lasso / n, t(x_basis_lasso) %*% Y)
      # Compute empirical Riesz representor
      theta_n_basis_row <- colMeans(x_basis_lasso_A1 - x_basis_lasso_A0)
      alpha <- x_basis_lasso %*% solve(t(x_basis_lasso)  %*% x_basis_lasso / n, theta_n_basis_row)

      IF <-  alpha * (Y - as.vector(x_basis_lasso %*% beta)) + as.vector(theta_n_basis_row %*% beta) - mean(theta_n_basis_row %*% beta)
      se <- sd(IF)/sqrt(n) * sqrt(n) / sqrt(n - ncol(x_basis_lasso))
      # Compute estimate of ATE
      est <- mean(theta_n_basis_row %*% beta)
      CI <- est + c(-1, 1) * abs(qnorm(1-0.975)) * se
      cover_approx <- as.numeric(approxATE >= CI[1] & approxATE <= CI[2])
      cover_true <- as.numeric(trueATE >= CI[1] & trueATE <= CI[2])
      sim_list$estimate <- c(sim_list$estimate , est)
      sim_list$CI <- c(sim_list$CI , CI)
      sim_list$cover_true <- c(sim_list$cover_true , cover_true)
      sim_list$cover_approx <- c(sim_list$cover_approx , cover_approx)

      sim_list$se <- c(sim_list$se, se)
      sim_list$trueATE <- trueATE
      sim_list$approxATE <- approxATE

      print("AdaLSE")
      print(est)
      print(se)
      print("AdaLSE CI")
      print(mean(sim_list$cover_true ))
      print(mean(sim_list$cover_approx))
      print(as.vector(table(sim_list$cover_true )))
      print(as.vector(table(sim_list$cover_approx )))

      # DR Competitor
      print("Fitting mean A")
      fit_A <- fit_hal(X, A,   family = "gaussian", smoothness_orders = 0, max_degree =2, num_knots=  nknots, fit_control = list(parallel = TRUE))
      EA <- predict(fit_A, new_data = X)
      print("Fitting variance A")
      if(!(all(A %in% c(0,1)))) {
        fit_varA <- fit_hal(X, (A- EA)^2,  family = "gaussian", smoothness_orders = 0, max_degree =2, num_knots=  nknots, fit_control = list(parallel = TRUE))
        VarA <- pmax(predict(fit_varA, new_data = X), 1e-3)
      } else {
        VarA <- pmax(EA * (1- EA), 1e-3)
      }

      DR <- mean(theta_n_A1 - theta_n_A0 + (A  - EA) / VarA * (Y - theta_n))
      DR_se <- sd(theta_n_A1 - theta_n_A0 + (A  - EA) / VarA * (Y - theta_n)) / sqrt(n)
      DR_CI <- DR + c(-1, 1) * abs(qnorm(1-0.975)) * DR_se

      DR_cover_approx <- as.numeric(approxATE >= DR_CI[1] & approxATE <= DR_CI[2])
      DR_cover_true<- as.numeric(trueATE >= DR_CI[1] & trueATE <= DR_CI[2])

      cover_true <- as.numeric(trueATE >= CI[1] & trueATE <= CI[2])
      sim_list$DR_estimate <- c(sim_list$DR_estimate,DR)
      sim_list$DR_se <- c(sim_list$DR_se, DR_se)
      sim_list$DR_CI <- c(sim_list$DR_CI,DR_CI)
      sim_list$DR_cover_true <- c(sim_list$DR_cover_true, DR_cover_true)
      sim_list$DR_cover_approx <- c(sim_list$DR_cover_approx, DR_cover_approx)

      sim_list$iter <- c(sim_list$iter, iter)
      sim_list$n <- c(sim_list$n, n)
      print("DR")
      print(DR)
      print(DR_se)
      print("DR CI")
      print(mean(sim_list$DR_cover_true ))
      print(mean(sim_list$DR_cover_approx))
      print(as.vector(table(sim_list$DR_cover_true )))
      print(as.vector(table(sim_list$DR_cover_approx )))
      print("se comparison")
      print(sd(sim_list$estimate))
      print(sd(sim_list$DR_estimate))
    }) # end try





  }
  return(sim_list)

}
