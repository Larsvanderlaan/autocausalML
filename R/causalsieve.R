



#' @import
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
    regression_fit = function() {
      return(private$.regression_fit)
    }),
  public = list(
    initialize = function(X, A, Y, g_basis_generator, nboots = 1000, bootstrap_indices = NULL) {
      X <- as.matrix(X)
      A <- as.vector(A)
      Y <- as.vector(Y)
      args <- sl3:::args_to_list()
      private$.args <- args
      self$bootstrap(nboots = nboots, bootstrap_indices = bootstrap_indices)
      private$.fit_outcome_regression()
    },

    add_target_parameter = function(formula =  g(A=1,X=X) - g(A=0,X=X) ~ 1 ,  name = NULL, bootstrap_se = TRUE  ) {

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
        X <- self$args$X
        A <- self$args$A
        m_A_X_g <- eval(parse(text=m_functional_string))
        if(is_intercept_msm) {
          return(colMeans(as.matrix(m_A_X_g )))
        }
        V <- model.matrix(formula_msm, as.data.frame(X))
        msm.fit <- glm.fit(V, m_A_X_g, family = gaussian(), intercept = FALSE)
        beta <- coef(msm.fit)
        return(beta)
      }

      m_functional_msm = function(X,A, g) {


        m_A_X_g <- eval(parse(text=m_functional_string))
        if(is_intercept_msm) {

          return(m_A_X_g - mean(m_A_X_g) )
        }

        V <- model.matrix(formula_msm, as.data.frame(X))

        msm.fit <- glm.fit(V, m_A_X_g, family = gaussian(), intercept = FALSE)
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

      return()

    },
    add_target_parameter_custom = function(theta_functional = NULL,  m_functional = NULL, name = NULL  ) {
      name_key <- paste0(name_key, collapse = "_")
      parameter_specification <- list(name = name, theta_functional = theta_functional, m_functional = m_functional)
      private$.target_parameters[[name_key]] <- parameter_specification
      return()
    },
    estimate = function(bootstrap_se = TRUE  ) {
      target_parameters <- self$target_parameters
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
        theta_n_basis <- do.call(cbind, lapply(seq_len(ncol(g_basis)), function(j) {
          g_fun <- function(X,A) {

            g_basis_gen(X=X,A=A)[,j]

          }

          return(theta_functional(X=X, A=A, g = g_fun))
        }))

        sieve_fit <- self$regression_fit
        coef_sieve <- sieve_fit$coef
        g_n <- sieve_fit$g_n

        theta_n <- theta_n_basis %*% coef_sieve
        estimate <- theta_n
        alpha_n_mat <- apply(theta_n_basis, 1, function(theta_n_basis_row) {
          alpha_n_coef <- solve(t(g_basis) %*% g_basis / n, theta_n_basis_row)
          alpha_n <- g_basis %*% alpha_n_coef
        })

        if( any(abs(theta_n - colMeans(alpha_n_mat * g_n)) > 1e-9)){
          print(theta_n - colMeans(alpha_n_mat * g_n))
          stop("Estimates dont match")
        }

        m_X_A_g_basis <-  lapply(seq_len(ncol(g_basis)), function(j) {
          g_fun <- function(X,A) {

            g_basis_gen(X=X,A=A)[,j]

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
        se <- sqrt(diag(as.matrix(var(IF))))/sqrt(n)
        parameter_est <- list(name = name, estimate = estimate, se = se, riesz_rep = alpha_n_mat , IF = IF,  theta_functional = theta_functional, m_functional = m_functional )


        private$.estimates[[name_key]] <- parameter_est

        if(bootstrap_se) {
          private$.bootstrap_parameter(name_key)
        }

        self$confint()
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
    confint = function(level = 0.95) {
      estimates <- self$estimates
      param_names <- names(estimates)
      for(name in param_names) {
        param <- estimates[[name]]
        estimate <- param$estimate
        se <- param$se
        se_boot <- param$se_boot

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
      g_sieve_fit <- glm.fit(g_basis, Y, family = gaussian())
      coef <- coef(g_sieve_fit)
      g_n <- as.vector(g_basis %*% coef)
      private$.regression_fit <- list(coef = coef, g_n = g_n, g_basis= g_basis )

      if(compute_bootstrap_MLE) {
        estimates_boot <-  lapply(1:self$args$nboots, function(iter) {
          tryCatch({
            boot_indices <- self$args$bootstrap_indices[[iter]]
            X_boot <- X[boot_indices, , drop = F]
            A_boot <- A[boot_indices]
            Y_boot <- Y[boot_indices]
            g_basis_boot <- g_basis[boot_indices, , drop = F]
            sieve_MLE_boot <- glm.fit(g_basis_boot, Y_boot, family = gaussian())
            beta_n_boot <- coef(sieve_MLE_boot)
            g_n_boot <- g_basis_boot %*% beta_n_boot
            return(list(X  = X_boot, A = A_boot, g_n = g_n_boot, coef = beta_n_boot))
          }, error = function(cond) {  })
          return(NULL)
        })
        private$.sieve_bootstrap_fits <- estimates_boot
      }
      return(self$regression_fit)
    },
    .bootstrap_parameter = function(name) {

      sieve_bootstrap_fits <- private$.sieve_bootstrap_fits

      bootstrap_estimates <- sapply(seq_along(sieve_bootstrap_fits), function(i) {
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
      bootstrap_estimates <- as.matrix(bootstrap_estimates)

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

#' @import
make_g_basis_generator_HAL <- function(X, A, Y, max_degree = 2, smoothness_orders = 1, num_knots =  c(30,10), formula_hal = NULL, screen_basis = TRUE, lambda= NULL, fit_control = list() ,...) {
  V <- as.matrix(cbind(X,A))
  library(hal9001)

  if(!screen_basis) {
    fit_control$cv_select <- FALSE
    lambda <- 1e-6
  }
  hal_fit <- hal9001::fit_hal(V, Y, formula = formula_hal,  family = "gaussian", lambda = lambda, fit_control = fit_control, max_degree = max_degree, smoothness_orders = smoothness_orders, num_knots =  num_knots, reduce_basis = 0.001)
  basis_list <- c(list(list(cols = 1, cutoffs = min(V) - 1, orders = 0)) , hal_fit$basis_list)

  basis_list_reduced <- basis_list[abs(hal_fit$coefs) > 0]




  g_basis_generator <- function(A,X) {
    V <- as.matrix(cbind(X,A))
    g_basis <- as.matrix(make_design_matrix(V, basis_list_reduced))
    return(g_basis)
  }
  return(g_basis_generator)

}


#' @import
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

