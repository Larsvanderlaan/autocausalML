


library(origami)
library(hal9001)
#
# dat <- generate_data_univariate(500, 0.5)
# X <- dat$X
# A <- dat$A
#  fit_pi <- fit_propensity_hal(X,A, num_knots = 25, smoothness_orders = 0)
#  alpha1 <- estimate_alpha1_cv(fit_pi$x_basisA1, A, fit_pi$coefs_matA1)
#  alpha0 <- estimate_alpha1_cv(fit_pi$x_basisA0, 1-A, fit_pi$coefs_matA0)
#  alpha <- ifelse(A==1, alpha1, alpha0)
#   truealpha <- 1 /ifelse(A==1, dat$pi, 1 - dat$pi )
#   quantile(alpha)
#

fit_propensity_hal <- function(X, A, smoothness_orders = 0,  num_knots = c(50,25), max_degree = 2, fit_control = list(), ...) {
  fit_control$cv_select <- FALSE
  X <- as.matrix(X)
  A <- as.vector(A)
  hal_fit <- fit_hal(X, A, family = "gaussian", smoothness_orders = smoothness_orders, num_knots = num_knots, max_degree =max_degree,   reduce_basis  = 25/length(A), fit_control = fit_control,...)
  basis_list <- c(list(list(cols = 1, orders = 0, cutoffs = min(X)-10)),  hal_fit$basis_list)
  keep_coefs <- ncol(hal_fit$coefs) > 1e-9
  coefs_mat <- as.matrix(hal_fit$coefs[,keep_coefs, drop = F])
  basis_list <- basis_list[keep_coefs]
  x_basis <- make_design_matrix(X, basis_list)
  reduce_basis_1 <- hal9001:::make_reduced_basis_map(x_basis,  25/length(A))
  reduce_basis_0 <- hal9001:::make_reduced_basis_map(as(1-x_basis, "sparseMatrix"),  25/length(A))
  reduce_basis <- intersect(reduce_basis_1, reduce_basis_0)
  reduce_basis <- sort(c(1,reduce_basis))
  x_basis <- as.matrix(x_basis)


  keep_basis <- reduce_basis
  tryCatch({
  tmp <- which(!is.na(coef(glm.fit(x_basis, A, family = gaussian()))))
  keep_basis <- intersect(keep_basis, tmp)
  }, error = function(cond){})


  coefs_matA1 <- coefs_mat[ (keep_basis), , drop = F]

  x_basisA1 <- x_basis[, keep_basis , drop = F]
  keep_basis <- reduce_basis
  tryCatch({
    tmp <- which(!is.na(coef(glm.fit(x_basis, A, family = gaussian()))))
    keep_basis <- intersect(keep_basis, tmp)
  }, error = function(cond){})


  coefs_matA0 <- coefs_mat[ (keep_basis), , drop = F]

  x_basisA0 <- x_basis[, keep_basis , drop = F]
  return(list(coefs_matA1 = coefs_matA1,  coefs_matA0 = coefs_matA0, x_basisA1 = x_basisA1, x_basisA0 = x_basisA0))
}


estimate_alpha1_cv <- function(x_basis, A, coefs_mat, folds =  origami::folds_vfold(length(A))){

  alpha_cv <- do.call(rbind,lapply(folds, estimate_alpha1_fold, x_basis = x_basis, A = A, coefs_mat = coefs_mat))
  val_index <- sapply(folds, function(fold){ validation()})
  cv_index <- which.min(colMeans(A[val_index]*alpha_cv^2 - 2*alpha_cv))
  #print(as.vector(colMeans(A[val_index]*alpha_cv^2 - 2*alpha_cv)))
  #print(cv_index)
  Aval <- A[val_index]
  oracle_risk <- colMeans( Aval*(alpha_cv - 1/pi[val_index])^2)
 # print(as.vector(oracle_risk))
  #print(which.min(oracle_risk))
  alpha_full <- estimate_alpha1(x_basis, A, coefs_mat[, cv_index, drop = F])
  alphan1 <-  (x_basis %*% alpha_full)



  return(alphan1)
}


estimate_alpha1_fold <- function(fold,x_basis, A, coefs_mat) {
  train_index <- training()
  validation_index <- validation()
  coefs_train <- estimate_alpha1(x_basis[train_index, , drop = F], A[train_index], coefs_mat)
  alpha_n_val <- x_basis[validation_index,, drop = F] %*% coefs_train
}

estimate_alpha1 <- function(x_basis, A, coefs_mat) {
  x_basis_full <- x_basis
  n <- nrow(x_basis_full)
  alpha_n1_lambda <-   do.call(cbind, lapply(1:ncol(coefs_mat), function(j) {
    alpha_coef <- rep(0, nrow(coefs_mat))
    keep <- coefs_mat[, j] != 0
    x_basis <- x_basis_full[,keep, drop = F]

    alpha_coef_part <- NULL
    tryCatch({

      alpha_coef_part <-  as.vector(solve(t((A)*x_basis) %*% x_basis, t(x_basis) %*%  rep(1,n)))
      alpha_coef[keep] <- alpha_coef_part
    }, error = function(cond){ })




    return(alpha_coef)
  }))

}
