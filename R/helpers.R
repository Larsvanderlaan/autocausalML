library(sl3)

compute_pi <- function(data, lrnr_pi = Lrnr_glm$new(family = binomial())) {
  task <- sl3_Task$new(data, covariates = grep( "^[WV]", colnames(data), value = TRUE), outcome = "A")

  lrnr_pi_trained <- lrnr_pi$train(task)
  pi <- lrnr_pi_trained$predict(task)
  pi <- pmax(pi, 0.001)
  pi <- pmin(pi, 0.999)
  return(pi)
}

compute_g <- function(data, lrnr_g = Lrnr_glm$new(family = binomial())) {
  task <- sl3_Task$new(data, covariates = c(grep( "^[WV]", colnames(data), value = TRUE)), outcome = "Y")
  A <- data$A
  g1 <- lrnr_g$train(task[A==1])$predict(task)
  g0 <- lrnr_g$train(task[A==0])$predict(task)

  return(list(g1 = g1, g0 = g0))

}

compute_TMLE <- function(data, pi, g1, g0,level = 0.05) {
  g <- ifelse(data$A==1, g1, g0)
  A <- data$A
  Y <- data$Y
  eps <- coef(glm.fit(A - (1-A), data$Y, weights = 1/ifelse(data$A==1, pi, 1-pi),
                      offset = g))
  g1_tmle <- g1 + eps
  g0_tmle <- g0 - eps
  est <- mean(g1_tmle - g0_tmle)
  IF <- g1_tmle - g0_tmle - est + (A/pi - (1-A)/(1-pi)) * (Y - g)
  se <- sd(IF)/sqrt(nrow(data))
  CI = est + qnorm(1-level/2) * se *c(-1,1)
  return(list(estimate = est, se = se, CI = CI))
}

compute_AIPW <- function(data, pi, g1, g0, level = 0.05) {
  g <- ifelse(data$A==1, g1, g0)
  A <- data$A
  Y <- data$Y
  est <-  mean(g1 - g0 + (A/pi - (1-A)/(1-pi)) * (Y - g))
  IF <- g1 - g0 - est + (A/pi - (1-A)/(1-pi)) * (Y - g)
  se <- sd(IF) /sqrt(nrow(data))
  CI = est + qnorm(1-level/2) * se *c(-1,1)
  return(list(estimate = est, se = se, CI = CI))
}

compute_glm<- function(data, level = 0.05) {
  data <- as.data.frame(data)
  glm_fit <- glm(Y ~ A + . , family = gaussian(), data = as.data.frame(data))

  sumry <- summary(glm_fit)$coef
  sumry <- sumry[which(rownames(sumry)=="A"),]
  est <- sumry[1]
  library(sandwich)
  se <-  diag(vcovHC(glm_fit, type = "HC"))^0.5
  se <- se[names(se)=="A"]
  CI <- est + qnorm(1-0.05/2) * se * c(-1,1)

  return(list(estimate = est, se = se, CI = CI))
}



compute_TMLE_CATE <- function(data, V, pi, g1, g0,level = 0.05) {
  g <- ifelse(data$A==1, g1, g0)
  pi <- pmax(pi, 1e-8)
  pi <- pmin(pi, 1 - 1e-8)
  print(range(pi))
  A <- data$A
  Y <- data$Y
  H <- V*(A - (1-A))
  eps <- coef(glm.fit(H, data$Y, weights = 1/ifelse(data$A==1, pi, 1-pi),
                      offset = g))
  g1_tmle <- g1 + V%*% eps
  g0_tmle <- g0 -  V %*% eps

  beta <- coef(glm.fit(V, g1_tmle - g0_tmle, family = gaussian() ))
  CATE <- as.vector(V  %*% beta)

  n <- length(A)
  scale <- (t(V) %*% V) / n
  scaleinv <- solve(scale)
  H <- H * (2*A-1) / ifelse(data$A==1, pi, 1-pi)
  EIF_Y <-   (H %*% scaleinv) * as.vector(Y - g)
  EIF_WA <- apply(V, 2, function(v) {
     (v * (g1_tmle - g0_tmle - CATE) - mean(v * (g1_tmle - g0_tmle - CATE)))
  }) %*% scaleinv


  IF <-   EIF_Y + EIF_WA

  se <- sqrt(diag(var(IF)))/sqrt(n)

  CI <- cbind( beta - qnorm(1-level/2) * se ,  beta + qnorm(1-level/2) * se )

  return(list(estimate = beta, se = se, CI = CI))
}



compute_TMLE_CATE_sp <- function(data, V,  pi, g1, g0,level = 0.05) {

  g <- ifelse(data$A==1, g1, g0)
  pi <- pmax(pi, 1e-8)
  pi <- pmin(pi, 1 - 1e-8)
  pi1 <- pi
  pi0 <- 1-pi
  A <- data$A
  Y <- data$Y

  beta <- coef(glm.fit(A*V, Y, offset = g0, family = gaussian()))
  g1 <- as.vector(g0 +   V %*% beta)
  g <- ifelse(data$A==1, g1, g0)
  var_Y1 <- var(Y )
  var_Y0 <- var(Y )
  var_Y <- ifelse(A==1, var_Y1, var_Y0)

  gradM <- V
  num <- gradM * (pi1 / var_Y1)
  denom <- (pi0 / var_Y0 + pi1 / var_Y1)
  hstar <- -num / denom
  H <- as.matrix((A * gradM + hstar) / var_Y)
  H1 <- as.matrix(( gradM + hstar) / var_Y1)
  H0 <- as.matrix((0 * gradM + hstar) / var_Y0)

  eps <- coef(glm.fit(H, data$Y,
                      offset = g))
  g1_tmle <- g1 + H1 %*% eps
  g0_tmle <- g0 + H0 %*% eps
  g_tmle <- g + H %*% eps

  scale <- apply(V, 2, function(v) {
    apply( H * (A * v), 2, mean)
  })
  scaleinv <- solve(scale)

  EIF_Y <-   (H %*% scaleinv) * as.vector(Y - g)




  beta <- coef(glm.fit(V, g1_tmle - g0_tmle, family = gaussian() ))


  n <- length(A)
  IF <-   EIF_Y

  se <- sqrt(diag(var(IF)))/sqrt(n)

  CI <- cbind( beta - qnorm(1-level/2) * se ,  beta + qnorm(1-level/2) * se )

  return(list(estimate = beta, se = se, CI = CI))
}

