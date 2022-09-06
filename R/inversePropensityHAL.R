#
# library(CVXR)
# fit_inverse_propensity_hal <- function(X, A, smoothness_orders = 1, num_knots = c(30,10), max_degree = 2, fit_control = list(), ...) {
#   fit_control$cv_select <- FALSE
#   X <- as.matrix(X)
#   A <- as.vector(A)
#   hal_fit <- fit_hal(X, A, family = "gaussian", smoothness_orders = 1, num_knots = num_knots, max_degree =max_degree, ..., fit_control = fit_control)
#   basis_list <- c(list(list(cols = 1, orders = 0, cutoffs = min(X)-10)),  hal_fit$basis_list)
#   basis_list <- basis_list[hal_fit$coefs[,ncol(hal_fit$coefs)] != 0]
#   x_basis <- as.matrix(make_design_matrix(X, basis_list))
#   coefs_mat <- hal_fit$coefs[hal_fit$coefs[,ncol(hal_fit$coefs)] != 0,, drop = F]
#   estimate_alpha_cv(x_basis, A, coefs_mat)
# }
#
# fit_outcome_regression_hal <- function(X, A, Y,  smoothness_orders = 1, num_knots = c(30,10), max_degree = 2, ...) {
#   hal_fit <- fit_hal(cbind(X,A), Y, family = "gaussian", smoothness_orders = smoothness_orders, num_knots = num_knots, max_degree =max_degree, ... )
#   g1 <- predict(hal_fit, new_data = cbind(X,1))
#   g0 <- predict(hal_fit, new_data = cbind(X,0))
#   return(list(g = ifelse(A==1, g1, g0), g1 = g1, g0 = g0))
# }
#
#
# estimate_alpha_cv <- function(x_basis, A, coefs_mat, folds =  origami::folds_vfold(length(A))){
#   folds <- origami::folds_vfold(n)
#   alpha_cv <- do.call(rbind,lapply(folds, estimate_alpha1_fold, x_basis = x_basis, A = A, coefs_mat = coefs_mat))
#   val_index <- sapply(folds, function(fold){ validation()})
#   cv_index <- which.min(colMeans(A[val_index]*alpha_cv^2 - 2*alpha_cv))
#   print((colMeans(A[val_index]*alpha_cv^2 - 2*alpha_cv)))
#   print(cv_index)
#   alpha_full <- estimate_alpha1(x_basis, A, coefs_mat)[, cv_index]
#   alphan1 <-  (x_basis %*% alpha_full)
#
#   folds <- origami::folds_vfold(n)
#   alpha_cv <- do.call(rbind,lapply(folds, estimate_alpha1_fold, x_basis = x_basis, A =1-A, coefs_mat = coefs_mat))
#   val_index <- sapply(folds, function(fold){ validation()})
#   cv_index <- which.min(colMeans((1-A)[val_index]*alpha_cv^2 - 2*alpha_cv))
#   print(cv_index)
#   alpha_full <- estimate_alpha1(x_basis, 1-A, coefs_mat)[, cv_index]
#   alphan0 <-  (x_basis %*% alpha_full)
#
#   return(list(inverse_pA1X = alphan1, inverse_pA0X = alphan0))
# }
#
#
# estimate_alpha1_fold <- function(fold,x_basis, A, coefs_mat) {
#   train_index <- training()
#   validation_index <- validation()
#   coefs_train <- estimate_alpha1(x_basis[train_index, , drop = F], A[train_index], coefs_mat)
#   alpha_n_val <- x_basis[validation_index,, drop = F] %*% coefs_train
# }
#
# estimate_alpha1 <- function(x_basis, A, coefs_mat) {
#   x_basis_full <- x_basis
#   n <- nrow(x_basis_full)
#   alpha_n1_lambda <-   do.call(cbind, lapply(1:ncol(coefs_mat), function(j) {
#     alpha_coef <- rep(0, nrow(coefs_mat))
#     keep <- coefs_mat[, j] != 0
#     x_basis <- x_basis_full[,keep, drop = F]
#
#     alpha_coef_part <- NULL
#     tryCatch({
#
#       alpha_coef_part <-  as.vector(solve(t((A)*x_basis) %*% x_basis, t(x_basis) %*%  rep(1,n)))
#       alpha_coef[keep] <- alpha_coef_part
#     }, error = function(cond){ })
#
#     if(is.null(alpha_coef_part)) {
#       try({
#       beta <- Variable(ncol(x_basis))
#       lasso <- function(beta, lambda) {
#         lambda * sum(abs(beta))
#       }
#       loss <- sum(((A)*(x_basis%*% beta)^2 - 2*x_basis %*% beta))
#       obj <- loss + lasso(beta, 0.00001)
#       prob <- Problem(Minimize(obj))
#       result <- solve(prob)
#       alpha_coef_part <- result$getValue(beta)
#       alpha_coef[keep] <- alpha_coef_part
#       })
#
#     }
#
#
#     return(alpha_coef)
#   }))
#
# }
#
#
# n <- 500
# d <- 1
# X <- replicate(d, runif(n, -1 , 1))
# X <- as.matrix(X)
# A <- rbinom(n, 1, plogis(2*rowMeans(X)))
# Y <- rnorm(n, rowMeans(X) + A + A*X[,1])
# library(hal9001)
# alpha_list <- fit_inverse_propensity_hal(X,A, num_knots = c(100,1))
# g_list <- fit_outcome_regression_hal(X,A,Y)
#
# alpha_ATE <- A*alpha_list$inverse_pA1X - (1-A)*alpha_list$inverse_pA0X
# eps <- coef(glm.fit(alpha_ATE, Y, family = gaussian(), offset = g_list$g))
#
# g1_tmle <- g_list$g1 + alpha_list$inverse_pA1X * eps
# g0_tmle <- g_list$g0 - alpha_list$inverse_pA0X * eps
# g_tmle <- g_list$g + alpha_ATE * eps
#
# mean(A*alpha_list$inverse_pA1X * Y) - mean((1-A)*alpha_list$inverse_pA0X * Y)
# mean(g1_tmle - g0_tmle)
#
#
# plot(1/plogis(2*rowMeans(X)),alpha_list$inverse_pA1X )
# plot(1/(1-plogis(2*rowMeans(X))),alpha_list$inverse_pA0X )
#
# plot(X,alpha_list$inverse_pA1X)
# plot(X,1/plogis(2*rowMeans(X)))
#
#
# plot(X,1/(1-plogis(2*rowMeans(X))))
# plot(X,alpha_list$inverse_pA0X)
