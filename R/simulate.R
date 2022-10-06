


generate_data_univariate <- function(n, const = 1.5) {

  W1 <- runif(n, -1 , 1)
  X <- as.matrix(W1)
  colnames(X) <- c("W1")
  pi <- plogis(const*(2*W1))
  A <- rbinom(n, 1, pi)
  CATE <- 1 + W1*cos(5*W1)
  g0 <- sin(5*W1)
  g1 <- g0 + CATE
  g <- A*g1 + (1-A)*g0
  Y <- rnorm(n, g, 0.2)
  return(list(X=X, A=A, Y=Y, data= as.data.frame(cbind(X,A,Y)), g1  = g1, g0 = g0, pi=pi,  ATE =  1))
}




generate_data_linear_fluct <- function(n, const = 1.5) {

  W1 <- runif(n, -1 , 1)
  W2 <- runif(n, -1 , 1)
  W3 <- runif(n, -1 , 1)
  W4 <- runif(n, -1 , 1)
  X <- cbind(W1, W2, W3, W4)
  colnames(X) <- c("W1", "W2", "W3", "W4")
  pi <- plogis(const*(W1+W2 + W3 + W4))
  A <- rbinom(n, 1, pi)
  g1 <-   2 * ( (1 )*(1 + W1 + W2 + W3 + W4)+ 3/sqrt(n) * (sin(5*W1) +  sin(5*W2) + sin(5*W3) + sin(5*W3)))
  g0 <-  (1 )*(1 + W1 + W2 + W3 + W4) + 3/sqrt(n) * (sin(5*W1) +  sin(5*W2) + sin(5*W3) + sin(5*W3))
  Y <- rnorm(n, A*g1 + (1-A)*g0, 1.5)
  return(list(X=X, A=A, Y=Y, data= as.data.frame(cbind(X,A,Y)), g1  = g1, g0 = g0, pi=pi,  ATE =  1))
}


generate_data_binary <- function(n, const = 1.5) {

  W1 <- runif(n, -1 , 1)
  W2 <- runif(n, -1 , 1)
  W3 <- runif(n, -1 , 1)
  X <- cbind(W1, W2, W3)
  colnames(X) <- c("W1", "W2", "W3")
  pi <- plogis(const*(W1+W2 + W3 ))
  A <- rbinom(n, 1, pi)
  CATE <- (0.3 + (W1 + W2 + W3 )/15)
  g0 <- 0.1 + 0.3*plogis((W1 + W2 + W3 ))
  g1 <- g0 +  CATE
  g1 <- pmax(g1, 0)
  g1 <- pmin(g1,1)

  Y <- rbinom(n, size = 1 , A*g1 + (1-A)*g0)
  return(list(X=X, A=A, Y=Y, data= as.data.frame(cbind(X,A,Y)), g1  = g1, g0 = g0, pi=pi,  ATE =  0.3))
}


generate_data_nonlinear_CATE <- function(n, const = 1.5) {

  W1 <- runif(n, -1 , 1)
  W2 <- runif(n, -1 , 1)
  W3 <- runif(n, -1 , 1)
  W4 <- runif(n, -1 , 1)
  X <- cbind(W1, W2, W3, W4)
  colnames(X) <- c("W1", "W2", "W3", "W4")
  pi <- plogis(const*(W1+W2 + W3 + W4))
  A <- rbinom(n, 1, pi)
  CATE <- 1 + W1 + W2 + W3 + W4
  g0 <- W1*sin(5*W1) + cos(5*W2) + sin(5*W3) + W4*cos(5*W4)
  g1 <- g0 + A * CATE
  Y <- rnorm(n,  A*g1 + (1-A)*g0  , 1.5)
  return(list(X=X, A=A, Y=Y, data= as.data.frame(cbind(X,A,Y)), g1  = g1, g0 = g0, pi=pi,  ATE =   1))
}


generate_data_linear_smallsample <- function(n, const = 1.5) {

  W1 <- runif(n, -1 , 1)
  W2 <- runif(n, -1 , 1)
  W3 <- runif(n, -1 , 1)
  W4 <- runif(n, -1 , 1)
  X <- cbind(W1, W2, W3, W4)
  colnames(X) <- c("W1", "W2", "W3", "W4")
  pi <- plogis(const*(W1+W2 + W3 + W4))
  A <- rbinom(n, 1, pi)
  g1 <- 1 + 1*(1 + W3   + W4 ) + W1 + W2
  g0 <- 1 + W1 + W2
  Y <- rnorm(n, A*g1 + (1-A)*g0, 0.3)
  return(list(X=X, A=A, Y=Y, data= as.data.frame(cbind(X,A,Y)), g1  = g1, g0 = g0, pi=pi,  ATE =  1))
}



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



generate_data_linear_lasso <- function(n, const = 1.5) {

  W1 <- runif(n, -1 , 1)
  W2 <- runif(n, -1 , 1)
  W3 <- runif(n, -1 , 1)
  W4 <- runif(n, -1 , 1)
  W5 <- runif(n, -1 , 1)
  W6 <- runif(n, -1 , 1)
  W7 <- runif(n, -1 , 1)
  W8 <- runif(n, -1 , 1)
  W9 <- runif(n, -1 , 1)
  W10 <- runif(n, -1 , 1)

  X <- cbind(W1, W2, W3, W4, W5, W6, W7, W8, W9, W10)
  colnames(X) <- c("W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9", "W10")
  pi <- plogis(const*(W1  + W2 + W3 + W5 + W7 + W9)/2)
  A <- rbinom(n, 1, pi)
  g1 <- 1 + 1*(1 + W2 + W3 + W6   )  + W2 + W3 + W6 + W8 + W10
  g0 <- 1 + W2 + W3 + W6 + W8 + W10
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
  pi <- plogis(const*(W1+W2 + W3 + W4))
  A <- rbinom(n, 1, pi)
  g1 <- 1 + 1 + 1*(abs(W1) + sin(5*W2) + abs(W3) + cos(5*W4)  ) + sin(5*W1) + abs(W2) + sin(5*W3) + 1/(1.25 + W4)
  g0 <- 1 + 0 + 0*(abs(W1) + sin(5*W2) + abs(W3) + cos(5*W4)  ) + sin(5*W1) + abs(W2) + sin(5*W3) + 1/(1.25 + W4)
  Y <- rnorm(n,  A*g1 + (1-A)*g0  , 1.5)
  return(list(X=X, A=A, Y=Y, data= as.data.frame(cbind(X,A,Y)), g1  = g1, g0 = g0, pi=pi,  ATE =   1.806832))
}




get_data_generator_univariate <- function(const) {
  out <- function(n) {
    generate_data_univariate(n, const)
  }
  out
}

get_data_generator_linear_fluct <- function(const) {
  out <- function(n) {
    generate_data_linear_fluct(n, const)
  }
  out
}


get_data_generator_binary <- function(const) {
  out <- function(n) {
    generate_data_binary(n, const)
  }
  out
}

get_data_generator_linear <- function(const) {
  out <- function(n) {
    generate_data_linear(n, const)
  }
  out
}

get_data_generator_linear_smallsample <- function(const) {
  out <- function(n) {
    generate_data_linear_smallsample(n, const)
  }
  out
}


get_data_generator_linear_lasso <- function(const) {
  out <- function(n) {
    generate_data_linear_lasso(n, const)
  }
  out
}

get_data_generator_nonlinear <- function(const) {
  out <- function(n) {
    generate_data_nonlinear(n, const)
  }
  out
}
get_data_generator_nonlinear_CATE <- function(const) {
  out <- function(n) {
    generate_data_nonlinear_CATE(n, const)
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


library(data.table)


run_sims <- function(const, n, nsims,  formula_hal = ~ h(.) + h(.,A), num_knots = c(1,1), screen_basis = F, gen_fun, lrnr_pi = Lrnr_glmnet$new(), lrnr_g = Lrnr_glmnet$new(formula = ~ . + A * .),nboots = 500) {

  out_list <- list()
  out <- lapply(1:nsims, function(iter) {
    try({
      print("Current settings:")
      print(n)
      print(const)
      print(iter)
      datam_list <- gen_fun(const)(n)

      X <- datam_list$X
      A <- datam_list$A
      Y <- datam_list$Y

      g_basis_gen <-make_g_basis_generator_HAL(X,A,Y, formula_hal =formula_hal, smoothness_orders=1, num_knots = num_knots, max_degree = 2, screen_basis = screen_basis)
      # datam_list_oracle <- gen_fun(const)(n=5000)
      # g_basis_oracle <- g_basis_gen(X=datam_list_oracle$X, A=datam_list_oracle$A)
      # print(dim(g_basis_oracle))
      # beta <- coef(glm.fit(g_basis_oracle, datam_list_oracle$Y, family = gaussian()))
      # ATE <- mean((g_basis_gen(X=datam_list_oracle$X, A= 1)-g_basis_gen(X=datam_list_oracle$X, A= 0) ) %*%beta)
      ATE <- datam_list$ATE

      causal_sieve <- causalsieve$new(X, A, Y, g_basis_gen, nboots = nboots)
      causal_sieve$add_target_parameter(g(A=1,X=X) - g(A=0,X=X) ~ 1, name = "ATE")
      # causal_sieve$add_target_parameter(g(A=1,X=X) - g(A=0,X=X) ~ 1 + W1)
      causal_sieve$estimate()
      #\causal_sieve$summary()
      name <- unlist(sapply(causal_sieve$estimates, `[[`, "name"))

      estimates <- unlist(sapply(causal_sieve$estimates, `[[`, "estimate"))
      CI_IF_df <- do.call(rbind, lapply(causal_sieve$estimates, `[[`, "CI"))
      causal_sieve$confint(include_se_df_correction = FALSE)
      CI_IF <- do.call(rbind, lapply(causal_sieve$estimates, `[[`, "CI"))
      CI_boot <- do.call(rbind, lapply(causal_sieve$estimates, `[[`, "CI_boot"))
      out <- cbind(t(as.data.table(c(iter, name))), t(as.data.table(as.numeric(c(estimates, CI_IF, CI_IF_df, CI_boot)))))
      print(t(as.data.table(as.numeric(c(estimates, CI_IF, CI_IF_df, CI_boot)))))

      colnames(out) <- c("iter", "name", "estimate", "CI_left", "CI_right", "CI_df_left", "CI_df_right", "CI_boot_left", "CI_boot_right")
      #

      data <- as.data.frame(cbind(X,A,Y))
      g_ests <- compute_g(data, lrnr_g = lrnr_g)
      g1 <- g_ests$g1
      g0 <- g_ests$g0
      pi <-compute_pi(as.data.frame(cbind(X,A,Y)), lrnr_pi = lrnr_pi)
      tmle <- compute_TMLE (data, pi, g1, g0,level = 0.05)
      aipw <- compute_AIPW (data, pi, g1, g0,level = 0.05)
      lm <- compute_glm (data,level = 0.05)

      comp <- as.numeric(unlist(c(tmle[-2], aipw[-2], lm[-2])))
      names(comp) <- c(paste0(c("estimate", "CI_left", "CI_right"), "_tmle"),
                       paste0(c("estimate", "CI_left", "CI_right"), "_aipw"),
                       paste0(c("estimate", "CI_left", "CI_right"), "_lm"))
      #
      out <- cbind(out, t(as.data.table(comp)))
      colnames(out) <- c(
        c("iter", "name", "estimate", "CI_left", "CI_right", "CI_df_left", "CI_df_right", "CI_boot_left", "CI_boot_right"),
        c(paste0(c("estimate", "CI_left", "CI_right"), "_tmle"),
          paste0(c("estimate", "CI_left", "CI_right"), "_aipw"),
          paste0(c("estimate", "CI_left", "CI_right"), "_lm"))
      )
      out_list[[iter]] <<- out
      out_full <- as.data.table(do.call(rbind, out_list))

      print("sieve IF")

      print(out_full[,mean(ATE >= CI_left & ATE <= CI_right), by = "name"][[2]])
      print(out_full[, mean(as.numeric(CI_right) - as.numeric(CI_left))])
      print("sieve IF - df adjusted")
      print(out_full[,mean(ATE >= CI_df_left & ATE <= CI_df_right), by = "name"][[2]])
      print(out_full[, mean(as.numeric(CI_df_right) - as.numeric(CI_df_left))])
      print("sieve IF - boot")
      print(    out_full[,mean(ATE >= CI_boot_left & ATE <= CI_boot_right), by = "name"][[2]])
      print(out_full[, mean(as.numeric(CI_boot_right) - as.numeric(CI_boot_left))])


      print("tmle")
      print(    out_full[,mean(ATE >= CI_left_tmle & ATE <= CI_right_tmle), by = "name"][[2]]
      )
      print(out_full[, mean(as.numeric(CI_right_tmle) - as.numeric(CI_left_tmle))])
      print("lm")
      print(    out_full[,mean(ATE >= CI_left_lm & ATE <= CI_right_lm), by = "name"][[2]]
      )
      print(out_full[, mean(as.numeric(CI_right_lm) - as.numeric(CI_left_lm))])

      return(out)
    })
  })

  out <- as.data.frame(do.call(rbind, out_list))
  out$const <- const
  out$n <- n
  return(out)

}

