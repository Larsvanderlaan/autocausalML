library(data.table)
library(doMC)
doMC::registerDoMC(cores = 11)
ATE_all <- fread(paste0("benchmark_datasets/ACICDataChallenge2019/", "lowDim_trueATE.csv"))



file_names <- paste0("low", 1:3200)
file_name <- file_names[1]
ns <- sapply(file_names, function(file_name) {
  data <- fread(paste0("benchmark_datasets/low_dimensional_datasets/", file_name, ".csv"))
  (nrow(data))
})

file_names <- rev(file_names[ns<=1000])
ATE_all <- ATE_all[ATE_all$filename %in% file_names,]
dgpid_list <- sort(unique(ATE_all$DGPid))

for(DGPid in dgpid_list) {




  keep <- ATE_all$DGPid == DGPid
  file_names <- ATE_all[keep,]$filename


  out_list <- list()
  out <- lapply(file_names, function(file_name){
    try({
      print(file_name)

      data <- fread(paste0("benchmark_datasets/low_dimensional_datasets/", file_name, ".csv"))
      print(dim(data))
      ATE <- ATE_all[ATE_all$filename==file_name,trueATE]

      X <- data[,grep("^V", colnames(data)), with = F]
      A <- data$A
      Y <- data$Y
      print(quantile(Y))

      tskA1 <- sl3_Task$new(data, covariates = grep("^V", colnames(data), value = T), "Y" )[A==1]
      tskA0 <- sl3_Task$new(data, covariates = grep("^V", colnames(data), value = T), "Y" )[A==0]

      lrnr <- Lrnr_screener_coefs$new(Lrnr_glmnet$new(family = gaussian()))
      lrnr <- Lrnr_xgboost$new(max_depth = 4, nrounds = 20)
      lrnr1 <- lrnr$train(tskA1)
      lrnr0 <- lrnr$train(tskA0)
      #keep <- union(lrnr1$fit_object$selected, lrnr0$fit_object$selected)

      imp1 <- lrnr1$importance()[,c(1,2)]
      imp0 <- lrnr0$importance()[,c(1,2)]
      imp <- rbind(imp1, imp0)
      imp <- unique(imp[, Gain := max(Gain), by = "Feature"])
      keep <-
        imp[[1]][order(-imp[[2]])[1:min(20,nrow(imp))]]


      print("original num cov")
      print(ncol(X))
      print("Screened num cov")
      print(length(keep))
      X <- X[,keep, with = F]
      data <- data[,c("A", "Y", keep), with = F]



      fit_hal_g_params <- list(fit_control = list(parallel = TRUE, cv_select = T), formula = ~ h(.) + h(.,.)  , smoothness_orders = 0, max_degree = 2, num_knots = c(30,10))
      g_basis_gen <- make_g_basis_generator_HAL(X,A,Y,  fit_hal_g_params = fit_hal_g_params, relaxed_fit = F)
      causal_sieve <- causalsieve$new(X, A, Y, g_basis_gen, nboots = 1000, weights = NULL)
      causal_sieve$add_target_parameter(g(A=1,X=X) - g(A=0,X=X) ~ 1, name = "ATE")
      causal_sieve$estimate()
      est <- unlist(sapply(causal_sieve$estimates, `[[`, "estimate"))
      CI <- unlist(do.call(rbind, lapply(causal_sieve$estimates, `[[`, "CI")))
      CI_boot <- unlist(do.call(rbind, lapply(causal_sieve$estimates, `[[`, "CI_boot")))
      bias <- est - ATE
      CI_IF_cov <- CI[1] <= ATE & CI[2] >= ATE
      CI_boot_cov <- CI_boot[1] <= ATE & CI_boot[2] >= ATE
      CI_width_IF <- diff(CI)
      CI_width_boot <- diff(CI_boot)

      print(CI_IF_cov)
      print(CI_boot_cov)


      g1 <- causal_sieve$regression_fit$g_n1
      g0 <- causal_sieve$regression_fit$g_n0
      print(head(data))
      pi <- compute_pi(data, lrnr_pi = Lrnr_hal9001$new(max_degree = 1, num_knots=30, fit_control = list(parallel = TRUE), smoothness_orders = 0, family = "binomial"))
      out_tmle <- compute_TMLE(data, pi, g1, g0)
      out_aipw <- compute_AIPW(data, pi, g1, g0)

      print(quantile(pi))

      print( out_tmle$CI[1] <= ATE & out_tmle$CI[2] >= ATE)
      print( out_aipw$CI[1] <= ATE & out_aipw$CI[2] >= ATE)

      print(c(ATE, est, out_tmle$est, out_aipw$est))
      dt <- data.table(file_name = file_name, bias=  bias ,
                       CI_boot_cov=CI_boot_cov,
                       CI_IF_cov = CI_IF_cov,
                       CI_width_IF = CI_width_IF, CI_width_boot = CI_width_boot,
                       bias_tmle = out_tmle$est - ATE,
                       CI_tmle = out_tmle$CI,
                       cov_tmle = out_tmle$CI[1] <= ATE & out_tmle$CI[2] >= ATE,
                       bias_aipw = out_aipw$est - ATE,
                       CI_aipw = out_aipw$CI,
                       cov_aipw = out_aipw$CI[1] <= ATE & out_aipw$CI[2] >= ATE

      )
      out_list[[file_name]] <<- dt
      tmp <- rbindlist(out_list, fill = TRUE)
      print("bias")
      print(c(mean(tmp$bias), mean(tmp$bias_tmle), mean(tmp$bias_aipw)))
      print("se")
      print(c(sd(tmp$bias), sd(tmp$bias_tmle), sd(tmp$bias_aipw)))

      print("coverage")
      print(
        colMeans(as.matrix(tmp[, c("CI_boot_cov", "CI_IF_cov", "cov_tmle", "cov_aipw")], with = F)
                 , na.rm = T))
      return(dt)
    })
  })
  out_list <- out_list[sapply(out_list, function(item){
    return(is.data.table(item))
  })]
  out_list <- rbindlist(out_list, fill = TRUE)

  fwrite(out_list, paste0("benchmark_datasets/results/Sims_HAL1s","_DGPid=",DGPid,  ".csv"))
}
