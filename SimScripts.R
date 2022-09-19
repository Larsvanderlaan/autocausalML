

library(data.table)
library(sl3)

run_sims <- function(const, n, nsims,  formula_hal = ~ h(.) + h(.,A), num_knots = c(1,1), screen_basis = F, gen_fun, lrnr_pi = Lrnr_glmnet$new(), lrnr_g = Lrnr_glmnet$new(formula = ~ . + A * .),nboots = 500) {
  print("HERE")
  datam_list <- gen_fun(1)(1000)

  X <- datam_list$X
  A <- datam_list$A
  Y <- datam_list$Y

  g_basis_gen <-make_g_basis_generator_HAL(X,A,Y, formula_hal =formula_hal, smoothness_orders=1, num_knots = num_knots, max_degree = 2, screen_basis = screen_basis)
  print("basis dimension:")
  print(dim(g_basis_gen(X,A)))

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


### HIGH DIM

out_list <- list()
outs <- lapply(c(  3,5, 8), function(const) {
  out_list[[as.character(const)]] <<-  list()
  lapply(rev(c(   250, 500, 1000  )) ,function(n) {

    out <- run_sims(const,n,1000,  formula_hal = ~ h(.) + h(.,A), num_knots = c(1,1), screen_basis = TRUE, gen_fun = get_data_generator_linear_lasso, lrnr_pi = Lrnr_glmnet$new(), lrnr_g = Lrnr_hal9001$new(formula = ~h(.)  , smoothness_orders = 1, max_degree =1, num_knots = c(1)))

    out_list[[as.character(const)]][[as.character(n)]] <<- out
    out2 <- rbindlist(unlist(out_list, recursive = F))
    fwrite(out2, file = "LassoHighDim.csv")
    return(out)

  })

})
outs <- rbindlist(unlist(outs, recursive = F))
fwrite(outs, file = "LassoHighDim.csv")

### SIMPLE

out_list <- list()
outs <- lapply(c(  3,5, 8), function(const) {
  out_list[[as.character(const)]] <<-  list()
  lapply(rev(c(   500, 1000,  2500 ,5000 )) ,function(n) {

    out <- run_sims(const,n,1000,  formula_hal = ~ h(.) + h(.,A), num_knots = c(20,20), screen_basis = TRUE, gen_fun = get_data_generator_linear, lrnr_pi = Lrnr_gam$new(), lrnr_g = Lrnr_hal9001$new(formula = ~h(.)  , smoothness_orders = 1, max_degree =2, num_knots = c(20)), nboots=500)

    out_list[[as.character(const)]][[as.character(n)]] <<- out
    out2 <- rbindlist(unlist(out_list, recursive = F))
    fwrite(out2, file = "SimpleParametricHAL2.csv")
    return(out)

  })

})



### COMPLEX


library(sl3)
out_list <- list()
outs <- lapply(c( 3,5,8), function(const) {
  out_list[[as.character(const)]] <<-  list()
  lapply(rev(c(   500, 1000,  2500 ,4000 )) ,function(n) {
    if(n == 4000){
      nknots <- 100
    } else if(n == 2500){
      nknots <- 75
    } else if(n == 1000){
      nknots <- 50
    } else if(n == 500){
      nknots <- 30
    }
    print(nknots)
    out <- run_sims(const,n,1000,  formula_hal = ~ h(.) + h(.,A), num_knots = c(nknots,nknots), screen_basis = TRUE, gen_fun = get_data_generator_nonlinear, lrnr_pi = Lrnr_gam$new(), lrnr_g = Lrnr_hal9001$new(formula = ~h(.)  , smoothness_orders = 1, max_degree =2, num_knots = c(nknots)), nboots=2)

    out_list[[as.character(const)]][[as.character(n)]] <<- out
    out2 <- rbindlist(unlist(out_list, recursive = F))
    fwrite(out2, file = "ComplexParametricHAL4_3.csv")
    return(NULL)

  })

})



### SMALL SAMPLE

out_list <- list()
outs <- lapply(c(  1, 4, 7), function(const) {
  out_list[[as.character(const)]] <<-  list()
  lapply(rev(c(  50, 100, 150, 200,  300 , 500)) ,function(n) {

    out <- run_sims(const,n,1000,  formula_hal = ~ h(.) + h(.,A), num_knots = c(1,1), screen_basis = TRUE, gen_fun = get_data_generator_linear_smallsample, lrnr_pi = Lrnr_glmnet$new(), lrnr_g =Lrnr_glmnet$new(), nboots=5000)

    out_list[[as.character(const)]][[as.character(n)]] <<- out
    out2 <- rbindlist(unlist(out_list, recursive = F))
    fwrite(out2, file = "SimpleParametricHALSmallSamples.csv")
    return(out)

  })

})

