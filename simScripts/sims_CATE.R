
library(autocausalML)
library(data.table)
library(sl3)
library(doParallel)
library(hal9001)
registerDoParallel(11)


run_sims_CATE <- function(const, n, nsims,   nboots = 2) {
  screen_basis = TRUE
  gen_fun <- get_data_generator_nonlinear_CATE
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
      screen_basis <- T
      g_basis_gen <-make_g_basis_generator_HAL(X,A,Y, formula_hal = ~ h(.) + h(.,A), smoothness_orders=1, num_knots = c(20, 20), max_degree = 2, screen_basis = screen_basis)
      print(dim(g_basis_gen(X=X,A=A)))
      ATE <- datam_list$ATE
      causal_sieve <- causalsieve$new(X, A, Y, g_basis_gen, nboots = nboots)
      causal_sieve$add_target_parameter(g(A=1,X=X) - g(A=0,X=X) ~ 1 + W1 + W2 + W3 + W4, name = "CATE")
      causal_sieve$estimate()
      name <- unlist(sapply(causal_sieve$estimates, `[[`, "name"))


      estimates <- unlist(sapply(causal_sieve$estimates, `[[`, "estimate"))
      CI_IF_df <- do.call(rbind, lapply(causal_sieve$estimates, `[[`, "CI"))
      causal_sieve$confint(include_se_df_correction = FALSE)
      CI_IF <- do.call(rbind, lapply(causal_sieve$estimates, `[[`, "CI"))
      CI_boot <- do.call(rbind, lapply(causal_sieve$estimates, `[[`, "CI_boot"))
      out1 <- cbind((as.data.table(cbind(iter, name))), (as.data.table( (cbind(estimates, CI_IF, CI_IF_df, CI_boot)))))
      colnames(out1) <- c("iter", "name", "estimate", "CI_left", "CI_right", "CI_df_left", "CI_df_right", "CI_boot_left", "CI_boot_right")



      g_basis_gen <-make_g_basis_generator_HAL(X,A,Y, formula_hal = ~ h(.) + h(.,A , k =1), smoothness_orders=1,  num_knots = c(20, 1), max_degree = 2, screen_basis = screen_basis)
      print(dim(g_basis_gen(X=X,A=A)))
      ATE <- datam_list$ATE
      causal_sieve <- causalsieve$new(X, A, Y, g_basis_gen, nboots = nboots)
      causal_sieve$add_target_parameter(g(A=1,X=X) - g(A=0,X=X) ~ 1 + W1 + W2 + W3 + W4, name = "CATE")
      causal_sieve$estimate()
      name <- unlist(sapply(causal_sieve$estimates, `[[`, "name"))


      estimates <- unlist(sapply(causal_sieve$estimates, `[[`, "estimate"))
      CI_IF_df <- do.call(rbind, lapply(causal_sieve$estimates, `[[`, "CI"))
      causal_sieve$confint(include_se_df_correction = FALSE)
      CI_IF <- do.call(rbind, lapply(causal_sieve$estimates, `[[`, "CI"))
      CI_boot <- do.call(rbind, lapply(causal_sieve$estimates, `[[`, "CI_boot"))
      out2 <- cbind((as.data.table(cbind(iter, name))), (as.data.table( (cbind(estimates, CI_IF, CI_IF_df, CI_boot)))))
      colnames(out2) <- c("iter", "name", "estimate_sp", "CI_left_sp", "CI_right_sp", "CI_df_left_sp", "CI_df_right_sp", "CI_boot_left_sp", "CI_boot_right_sp")


      lrnr_pi <- Lrnr_glmnet$new()
      lrnr_g <-  Lrnr_hal9001$new(formula = ~ h(.) , family = "gaussian", smoothness_orders = 1, max_degree = 1, num_knots = c(20, 20))
      g_ests <- compute_g(as.data.frame(cbind(X,A,Y)), lrnr_g = lrnr_g)
      g1 <- g_ests$g1
      g0 <- g_ests$g0
      pi <-compute_pi(as.data.frame(cbind(X,A,Y)), lrnr_pi = lrnr_pi)
      npglm <- compute_TMLE_CATE(data = as.data.frame(cbind(X,A,Y)), V = cbind(1,X), g1 = g1, g0 =  g0, pi = pi, level = 0.05)
      out_np <- cbind(npglm[[1]], npglm[[3]][,1] , npglm[[3]][,2])
     # npglm <- causalglm::npglm(~ 1 + W1 + W2 + W3 + W4, data = as.data.frame(cbind(X,A,Y)), W = colnames(X), A = "A", Y = "Y", estimand = "CATE", sl3_Learner_A = Lrnr_gam$new(), sl3_Learner_Y = Lrnr_hal9001$new(formula = ~ h(.) + h(.,A), family = "gaussian", smoothness_orders = 1, max_degree = 2, num_knots = c(20, 20)), verbose = T)
      colnames(out_np) <- c("estimate_npglm", "CI_left_npglm", "CI_right_npglm")

      spglm <- compute_TMLE_CATE_sp(data = as.data.frame(cbind(X,A,Y)), V = cbind(1,X), g1 = g1, g0 =  g0, pi = pi, level = 0.05)
      out_sp <- cbind(spglm[[1]], spglm[[3]][,1] , spglm[[3]][,2])
      #spglm <- causalglm::spglm(~ 1 + W1 + W2 + W3 + W4, data = as.data.frame(cbind(X,A,Y)), W = colnames(X), A = "A", Y = "Y", estimand = "CATE", sl3_Learner_A = Lrnr_glmnet$new(), sl3_Learner_Y =  Lrnr_hal9001$new(formula = ~ h(.), smoothness_orders = 1, max_degree = 1, num_knots = c(20, 1)), verbose = FALSE, append_interaction_matrix = F)
      colnames(out_sp) <- c("estimate_spglm", "CI_left_spglm", "CI_right_spglm")
      print(out_sp)
       out <- cbind(out1, out2[,-c(1,2)], out_np, out_sp)


      out_list[[iter]] <<- out
      out_full <- as.data.table(do.call(rbind, out_list))


       print("sieve IF - df adjusted")
      print(out_full[,mean(ATE >= CI_df_left & ATE <= CI_df_right), by = "name"][[2]])
      print(out_full[,mean(ATE >= CI_df_left_sp & ATE <= CI_df_right_sp), by = "name"][[2]])
      print(out_full[,mean(ATE >= CI_left_npglm & ATE <= CI_right_npglm), by = "name"][[2]])
      print(out_full[,mean(ATE >= CI_left_spglm & ATE <= CI_right_spglm), by = "name"][[2]])
      print(out_full[, sd(estimate), by = "name" ][[2]])
      print(out_full[, sd(estimate_sp), by = "name" ][[2]])
      print(out_full[, sd(estimate_npglm), by = "name" ][[2]])
      print(out_full[, sd(estimate_spglm), by = "name" ][[2]])

      return(out)
    })
  })

  out <- as.data.frame(do.call(rbind, out_list))
  out$const <- const
  out$n <- n
  return(out)

}






library(sl3)
out_list <- list()
outs <- lapply(c(   4, 7,1), function(const) {
  out_list[[as.character(const)]] <<-  list()
  lapply(rev(c(   500, 1000,  2500, 5000   )) ,function(n) {

    out <- run_sims_CATE(const,n,5000)

    out_list[[as.character(const)]][[as.character(n)]] <<- out
    out2 <- rbindlist(unlist(out_list, recursive = F))
    fwrite(out2, file = "SimsHALCATE.csv")
    return(out)

  })

})



