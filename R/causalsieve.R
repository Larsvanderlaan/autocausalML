



#' @export
causalsieve <- R6::R6Class(
  classname = "causalsieve",
  portable = TRUE, class = TRUE,
  active = list(
    args = function(){
      private$.args
    },
    target_parameters = function(){
      return(private$.target_parameters)
    },
    estimates = function() {
      return(private$.estimates)
    },
    weights = function(){
      private$.weights
    },
    regression_fit = function() {
      return(private$.regression_fit)
    }),
  public = list(
    initialize = function(X, A, Y, g_basis_generator, nboots = 1000, bootstrap_indices = NULL, weights = NULL, adjust_for_variance = FALSE) {
      X <- as.matrix(X)
      A <- as.vector(A)
      Y <- as.vector(Y)
      if(length(colnames(X)) == 0) {
        colnames(X) <- paste0("X", 1:ncol(X))
      }
      args <- sl3:::args_to_list()
      private$.args <- args
      if(adjust_for_variance & is.null(weights)) {
        private$.fit_outcome_regression(compute_bootstrap_MLE = FALSE)
        g <- self$regression_fit$g_n
        if(all(Y %in% c(0,1))) {
          g <- pmax(g, 1e-3)
          g <- pmin(g, 1-1e-3)
          weights <- 1 / (g * (1 - g))
          weights <- pmin(weights, 10 * mean(weights))
        } else {

          data <- as.data.frame(cbind(X,A,  (Y - g)^2))

          colnames(data) <- c(colnames(X), "A", "res")
          task <- sl3_Task$new(data, covariates = setdiff(colnames(data), "res"), outcome = "res", outcome_type = "continuous")
          lrnr_cv <- Lrnr_sl$new(Stack$new( Lrnr_gam$new(family = gaussian()),
                                            Lrnr_glmnet$new(family = "gaussian"),
                                            Lrnr_hal9001$new(smoothness_orders = 0, max_degree = 1, num_knots = c(30), family = "gaussian" )
          ),
          meta_learner = Lrnr_cv_selector$new(loss_squared_error))
          lrnr_cv <- lrnr_cv$train(task)
          weights <- lrnr_cv$predict(task)
          weights <- pmax(weights, 1e-5)
          weights <- 1/weights
          weights <- pmin(weights, 20 * mean(weights))


        }
      }
      if(is.null(weights)) {
        weights <- rep(1,length(A))
      }

      private$.weights <- weights
      self$bootstrap(nboots = nboots, bootstrap_indices = bootstrap_indices)
      private$.fit_outcome_regression(compute_bootstrap_MLE = TRUE)
    },

    add_target_parameter = function(formula =  g(A=1,X=X) - g(A=0,X=X) ~ 1 ,  name = NULL, bootstrap_se = TRUE  ) {
      n <- length(self$args$A)
      m_functional_string <- as.character(formula)[2]
      msm_model_string <- as.character(formula)[3]
      formula_msm <-  as.formula(paste0("~", msm_model_string))


      if(is.null(name)) {

        name_key <- paste0(m_functional_string, " ~ " ,msm_model_string)
        terms <- attr(terms(formula_msm, data=as.data.frame(self$args$X)), "term.labels")
        if( attr(terms(formula_msm,  data=as.data.frame(self$args$X)), "intercept")) {terms <- c("intercept", terms)}
        name <- paste0(name_key, "_",terms)
      } else {
        if(length(name) > 1) {
          name_key <- paste0(name, collapse = "_")
        } else {

          name_key <- name
          var_terms <- attr(terms(formula_msm, data=as.data.frame(self$args$X)), "term.labels")
          if( attr(terms(formula_msm,  data=as.data.frame(self$args$X)), "intercept")) {var_terms <- c("intercept", var_terms)}
          name <- paste0(name, "_",var_terms)

        }


      }



      is_intercept_msm <- length(attr(terms(formula_msm, data = data.frame( (self$args$X))), "term.labels")) ==0

      theta_functional_msm = function(X,A, g) {
        #X <- self$args$X
        #A <- self$args$A
        m_A_X_g <- eval(parse(text=m_functional_string))
        if(is_intercept_msm) {
          return(colMeans(as.matrix(m_A_X_g )))
        }
        V <- model.matrix(formula_msm, as.data.frame(X))
        msm.fit <- fastglm::fastglm(V, m_A_X_g, family = gaussian(), intercept = FALSE)
        beta <- coef(msm.fit)
        return(beta)
      }

      m_functional_msm = function(X,A, g) {


        m_A_X_g <- eval(parse(text=m_functional_string))
        if(is_intercept_msm) {

          return(m_A_X_g - mean(m_A_X_g) )
        }

        V <- model.matrix(formula_msm, as.data.frame(X))

        msm.fit <- fastglm::fastglm(V, m_A_X_g, family = gaussian(), intercept = FALSE)
        beta <- coef(msm.fit)

        m_A_V_g <-  as.vector(V %*% beta)

        scale <- (t(V) %*% V) /n
        scaleinv <- solve(scale)


        m_X_A_g_msm <-  (V* (m_A_X_g - m_A_V_g)  ) %*% scaleinv
        return(m_X_A_g_msm)
      }




      parameter_specification <- list(name = name, formula = formula, theta_functional = theta_functional_msm, m_functional = m_functional_msm)
      private$.target_parameters[[name_key]] <- parameter_specification
      # return(self$add_target_parameter_custom(theta_functional = theta_functional_msm, m_functional = m_functional_msm, name = name, bootstrap_se = bootstrap_se   ) )

      return(invisible())

    },
    add_target_parameter_custom = function(theta_functional = NULL,  m_functional = NULL, name = NULL  ) {
      name_key <- paste0(name_key, collapse = "_")
      parameter_specification <- list(name = name, theta_functional = theta_functional, m_functional = m_functional)
      private$.target_parameters[[name_key]] <- parameter_specification
      return(invisible())
    },
    estimate = function(bootstrap_se = TRUE , include_se_df_correction = TRUE ) {
      target_parameters <- self$target_parameters
      X <- self$args$X
      A <- self$args$A
      Y <- self$args$Y
      n <- length(A)
      weights <- self$weights
      for(name_key in names(target_parameters)) {
        param <- target_parameters[[name_key]]
        name <- param$name
        theta_functional <- param$theta_functional
        m_functional <- param$m_functional
        name <- param$name


        g_basis_gen <- self$args$g_basis_generator
        g_basis <- g_basis_gen(X=X,A=A)
        # Compute theta_n functional at each basis function
        # If theta_n functional is vector-valued then a matrix of dimension (len(theta_n) x nbasis)
        cache_list <- list()
        theta_n_basis <- do.call(cbind, lapply(seq_len(ncol(g_basis)), function(j) {


          g_fun <- function(X,A) {
            uuid <- digest::digest(as.data.frame(cbind(X,A)))

            g_basis <- cache_list[[uuid]]
            if(is.null(g_basis)) {

              g_basis <- g_basis_gen(X=X,A=A)
              cache_list[[uuid]] <<- g_basis
            }
            g_basis[,j]
          }

          return(theta_functional(X=X, A=A, g = g_fun))
        }))

        sieve_fit <- self$regression_fit
        coef_sieve <- sieve_fit$coef
        g_n <- sieve_fit$g_n

        theta_n <- theta_n_basis %*% coef_sieve
        estimate <- theta_n
        alpha_n_mat <- apply(theta_n_basis, 1, function(theta_n_basis_row) {
          alpha_n_coef <- NULL

          tryCatch({alpha_n_coef <- solve(t(g_basis) %*% diag(self$weights) %*% g_basis / n, theta_n_basis_row)}, error = function(cond){print(cond)})






          alpha_n <- self$weights * (g_basis %*% alpha_n_coef)
          if(!all(abs(colMeans(as.vector(alpha_n)*g_basis) - theta_n_basis_row) <= 1e-6)) {
            print(colMeans(as.vector(alpha_n)*g_basis) - theta_n_basis_row)
            warning("Riez-rep Sieve-MLE did not solve all scores. g_basis might be singular.")
          }
          return(alpha_n)
        })


        if( any(abs(theta_n - colMeans(alpha_n_mat * g_n)) > 1e-5)){
          print(theta_n - colMeans(alpha_n_mat * g_n))
          print("Estimates Riesz-based and sieve-plugin-based estimates dont match. Contact developer.")
        }

        m_X_A_g_basis <-  lapply(seq_len(ncol(g_basis)), function(j) {
          g_fun <- function(X,A) {
            uuid <- digest::digest(as.data.frame(cbind(X,A)))

            g_basis <- cache_list[[uuid]]
            if(is.null(g_basis)) {

              g_basis <- g_basis_gen(X=X,A=A)
              cache_list[[uuid]] <<- g_basis
            }
            g_basis[,j]
          }

          return(coef_sieve[j] * as.matrix(m_functional(X=X,A=A, g=g_fun)))
        })

        m_X_A_g_n <- as.matrix(purrr::reduce(m_X_A_g_basis, `+`))


        m_X_A_g_mean_mat <- do.call(cbind,lapply(1:nrow(m_X_A_g_n), function(k){
          colMeans(as.matrix(m_X_A_g_n))
        } ))
        if(is.vector(m_X_A_g_mean_mat)) {
          m_X_A_g_mean_mat <- as.matrix(m_X_A_g_mean_mat)
        }
        if(nrow(m_X_A_g_mean_mat)!= nrow(m_X_A_g_n)) {
          m_X_A_g_mean_mat <- t(m_X_A_g_mean_mat)
        }

        IF <- m_X_A_g_n - m_X_A_g_mean_mat  + alpha_n_mat * as.vector(Y - g_n)
        df <- 0
        se <- sqrt(diag(as.matrix(var(IF))) * (n/ (n - df)))/sqrt(n)
        parameter_est <- list(name = name, estimate = estimate, se = se, riesz_rep = alpha_n_mat , IF = IF,  theta_functional = theta_functional, m_functional = m_functional )


        private$.estimates[[name_key]] <- parameter_est

        if(bootstrap_se) {
          private$.bootstrap_parameter(name_key)
        }

        self$confint( include_se_df_correction = include_se_df_correction)
      }
      return(invisible(private$.estimates))
    },
    summary = function(ndigits = 5) {
      out <- as.data.frame(do.call(rbind,lapply(self$estimates, function(item){
        item <- item[c("name", "estimate", "se", "CI", "se_boot", "CI_boot")]
        names(item) <- c("Param", "Estimate", "se", "CI", "se_boot", "CI_boot")
        item$CI <- signif(item$CI, ndigits)
        item$CI_boot <- signif(item$CI_boot, ndigits)
        item$se <- signif(item$se, ndigits)
        item$se_boot <- signif(item$se_boot, ndigits)
        item$Estimate <- signif(item$Estimate, ndigits)
        item$CI <- paste0("(", item$CI[, 1], ",", item$CI[,2], ")")
        item$CI_boot <- paste0("(", item$CI_boot[,1], ",", item$CI_boot[,2], ")")
        data.frame(Param = item$Param,Estimate = item$Estimate, se = item$se, CI = item$CI, se_boot = item$se_boot, CI_boot = item$CI_boot )
      })))
      rownames(out) <- NULL
      return(out)
    },
    confint = function(level = 0.95, include_se_df_correction = TRUE) {
      estimates <- self$estimates
      param_names <- names(estimates)
      n <- length(self$args$A)
      for(name in param_names) {
        param <- estimates[[name]]
        estimate <- param$estimate
        se <- param$se
        se_boot <- param$se_boot

        if(include_se_df_correction) {

          df <- ncol(self$regression_fit$g_basis)
        } else {
          df <- 0
        }
        se <- se * sqrt(n / (n-df))
        CI <- cbind(estimate - qnorm( 1 - (1-level)/2) * se , estimate + qnorm( 1 - (1-level)/2) * se )
        CI_boot <-  cbind(estimate - qnorm( 1 - (1-level)/2) * se_boot , estimate + qnorm( 1 - (1-level)/2) * se_boot )
        param$CI <- CI
        param$CI_boot <- CI_boot
        estimates[[name]] <- param
      }

      private$.estimates <- estimates
      return(estimates)

    },
    bootstrap = function(nboots = 5000, bootstrap_indices = NULL) {
      n <- length(private$.args$A)
      if(!is.null(bootstrap_indices)) {
        nboots <- length(bootstrap_indices )
      }
      if(is.null(bootstrap_indices)) {
        bootstrap_indices <- lapply(1:nboots, function(i) {
          sample(seq_len(n), n, replace = TRUE)})
      }
      private$.args$nboots <- nboots
      private$.args$bootstrap_indices <- bootstrap_indices
      private$.fit_outcome_regression(g_basis_generator = NULL, compute_bootstrap_MLE = TRUE)
    }
  ),
  private = list(
    .regression_fit = NULL,
    .sieve_bootstrap_fits = NULL,
    .target_parameters = list(),
    .estimates = list(),
    .args = list(),
    .weights = NULL,
    .fit_outcome_regression = function(g_basis_generator = NULL, compute_bootstrap_MLE = TRUE) {
      if(!is.null(g_basis_generator)) {
        self$args$g_basis_generator <- g_basis_generator
      } else {
        g_basis_generator <- self$args$g_basis_generator
      }
      X <- self$args$X
      A <- self$args$A
      Y <- self$args$Y

      g_basis <- g_basis_generator(X=X,A=A)

      g_sieve_fit <- glm.fit(g_basis, Y, family = gaussian(), weights = self$weights)


      coef <- coef(g_sieve_fit)
      remove_basis <- is.na(coef)
      if(any(remove_basis)) {
        g_basis_generator_new <- function(A,X) {
          g_basis_generator(X=X,A=A)[,!remove_basis, drop = F]
        }
        private$.args$g_basis_generator <- g_basis_generator_new
      }

      g_basis <- g_basis[, !remove_basis, drop = F]

      g_sieve_fit <- glm.fit(g_basis, Y, family = gaussian(), weights = self$weights)
      coef <- coef(g_sieve_fit)


      g_basis_generator2 <- private$.args$g_basis_generator
      g_n <- as.vector(g_basis %*% coef)
      #print(quantile(g_basis_generator2(X=X,A=1)%*% coef))
      #print(quantile(g_basis_generator2(X=X,A=0)%*% coef))
      if(!all(abs(colMeans( g_basis * (Y- g_n)) )<= 1e-8)) {
        warning("Sieve-MLE did not solve all scores. g_basis might be singular.")
      }

      private$.regression_fit <- list(coef = coef, g_n = g_n, g_basis= g_basis, remove_basis = remove_basis )

      if(compute_bootstrap_MLE) {
        estimates_boot <-  lapply(1:self$args$nboots, function(iter) {
          tryCatch({
            boot_indices <- self$args$bootstrap_indices[[iter]]
            X_boot <- X[boot_indices, , drop = F]
            A_boot <- A[boot_indices]
            Y_boot <- Y[boot_indices]
            g_basis_boot <- g_basis[boot_indices, , drop = F]
            sieve_MLE_boot <- fastglm::fastglm(g_basis_boot, Y_boot, family = gaussian(), weights = self$weights)
            beta_n_boot <- coef(sieve_MLE_boot)
            beta_n_boot[is.na(beta_n_boot)] <- 0
            #g_n_boot <- g_basis_boot %*% beta_n_boot
            return(list(    coef = beta_n_boot))
          }, error = function(cond) {  })
          return(NULL)
        })
        private$.sieve_bootstrap_fits <- estimates_boot
      }
      return(self$regression_fit)
    },
    .bootstrap_parameter = function(name) {

      sieve_bootstrap_fits <- private$.sieve_bootstrap_fits

      bootstrap_estimates <- sapply(seq_along(sieve_bootstrap_fits), function(i) {try({
        indices <- self$args$bootstrap_indices[[i]]


        g <- self$args$g_basis_generator
        X <- self$args$X[indices , , drop=F]
        A <- self$args$A[indices]
        n <- length(A)

        theta_functional <- private$.estimates[[name]]$theta_functional
        g_basis_gen <- self$args$g_basis_generator


        coef_boot <- sieve_bootstrap_fits[[i]]$coef
        g_fun <- function(X,A) {
          g_basis_gen(X=X,A=A) %*% coef_boot
        }
        bootstrap_estimate <- theta_functional(X=X,A=A,g= g_fun)

        return(bootstrap_estimate)
      })
        return(NA)
      })
      bootstrap_estimates <- na.omit(as.matrix(bootstrap_estimates))

      if(nrow(bootstrap_estimates)!=length(sieve_bootstrap_fits)){
        bootstrap_estimates <- t(bootstrap_estimates)
      }



      se_boot <- sqrt(diag(as.matrix(var(bootstrap_estimates))))

      private$.estimates[[name]]$se_boot <- se_boot
      #private$.estimates[[name]]$CI_boot <- CI_boot
      return(list(bootstrap_estimates = bootstrap_estimates, se = se_boot))
    }
  )
)

#' @export
make_g_basis_generator_HAL <- function(X, A, Y, fit_hal_g_params = list(formula = ~ h(.) + h(.,A), reduce_basis = 25/length(A), smoothness_orders = 1, max_degree =2, num_knots = c(100,100)), fit_hal_alpha_params = list(  smoothness_orders = 0, max_degree =1, num_knots = c(25,5)), screen_basis = TRUE, relaxed_fit = TRUE, weight_screen_by_alpha = TRUE, ...) {

  V <- as.matrix(cbind(X,A))
  library(hal9001)
  if(all(A %in% c(0,1)) & weight_screen_by_alpha) {

   fit_hal_alpha_params$X <- X
   fit_hal_alpha_params$A <- A
   fit_pi <- sl3:::call_with_args(fit_propensity_hal, fit_hal_alpha_params)
  alpha1 <- estimate_alpha1_cv(fit_pi$x_basisA1, A, fit_pi$coefs_matA1)
  alpha0 <- estimate_alpha1_cv(fit_pi$x_basisA0, 1-A, fit_pi$coefs_matA0)
  weights <- ifelse(A==1, alpha1, alpha0)
  } else {
    weights <- rep(1, length(A))
  }
  weights <- pmin(10, weights)

  print(quantile(weights))

  if(!screen_basis) {
    fit_hal_g_params$fit_control$cv_select <- FALSE
    fit_hal_g_params$lambda <- 1e-15
  }
  fit_hal_g_params$fit_control$relax <- relaxed_fit
  fit_hal_g_params$fit_control$gamma <- 0

  fit_hal_g_params$fit_control$weights <- weights
  basis_formula <- formula_hal(fit_hal_g_params$formula, smoothness_orders = fit_hal_g_params$smoothness_orders, num_knots = fit_hal_g_params$num_knots, X = as.data.frame(V) )$basis_list
  fit_hal_g_params$basis_list <- basis_formula
  fit_hal_g_params$X <- V
  fit_hal_g_params$Y <- Y
  fit_hal_g_params$family = "gaussian"
  fit_hal_g_params$reduce_basis = 25/length(A)
  hal_fit <- sl3:::call_with_args( hal9001::fit_hal, fit_hal_g_params)

  if(relaxed_fit) {
    lambda.min <- hal_fit$lasso_fit$relaxed$lambda.min
    coefs <- as.vector(coef(hal_fit$lasso_fit$glmnet.fit, s = lambda.min))
   # print(hal_fit$lasso_fit$relaxed$lambda.min)
    #print(hal_fit$lambda_star)
  } else {
    coefs <- as.vector(hal_fit$coefs)
  #  print(hal_fit$lambda_star)
  }





  basis_list <- c(list(list(cols = 1, cutoffs = min(V) - 1, orders = 0)) , hal_fit$basis_list)

  # print("BASIS before")
  # print(length(basis_list))
  # #print(((basis_list)))
  basis_list_reduced <- basis_list[abs(coefs) > 0]
  # print("BASIS SELECTED")
  # print(length(basis_list_reduced))

  basis_list <- basis_list_reduced



  g_basis_generator <- function(A,X) {
    V <- as.matrix(cbind(X,A))
    g_basis <- as.matrix(make_design_matrix(V, basis_list_reduced, 0.9))

    return(g_basis)
  }
  return(g_basis_generator)

}


#' @export
make_g_basis_generator_LASSO <- function(X, A, Y, formula = ~., use_lambda_min = TRUE, ...) {

  V <- model.matrix(formula, as.data.frame(cbind(X,A)))

  cvglmnet.fit <- glmnet::cv.glmnet(V, Y, family = "gaussian", intercept = FALSE, ...)
  if(use_lambda_min) {
    beta <- coef(cvglmnet.fit$glmnet.fit, s = cvglmnet.fit$lambda.min )[-1]
  } else {
    beta <- coef(cvglmnet.fit$glmnet.fit, s = cvglmnet.fit$lambda.min )[-1]
  }
  if(length(beta) != ncol(V)) {
    stop("Something went wrong.")
  }
  beta <- as.vector(beta)
  keep <- (beta!=0)


  g_basis_generator <- function(A,X) {

    g_basis <- model.matrix(formula, as.data.frame(cbind(X,A)))[, keep, drop = FALSE]
    return(g_basis)
  }
  return(g_basis_generator)

}

