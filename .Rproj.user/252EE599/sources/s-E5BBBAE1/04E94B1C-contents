compute_pi <- function(data, lrnr_pi = Lrnr_glm$new(family = binomial())) {
  task <- sl3_Task$new(data, covariates = grep( "^V", colnames(data), value = TRUE), outcome = "A")

  lrnr_pi_trained <- lrnr_pi$train(task)
  pi <- lrnr_pi_trained$predict(task)
  pi <- pmax(pi, 0.001)
  pi <- pmin(pi, 0.999)
  return(pi)
}

compute_TMLE <- function(data, pi, g1, g0,level) {
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

compute_AIPW <- function(data, pi, g1, g0, level) {
  g <- ifelse(data$A==1, g1, g0)
  A <- data$A
  Y <- data$Y
  est <-  mean(g1 - g0 + (A/pi - (1-A)/(1-pi)) * (Y - g))
  IF <- g1 - g0 - est + (A/pi - (1-A)/(1-pi)) * (Y - g)
  se <- sd(IF) /sqrt(nrow(data))
  CI = est + qnorm(1-level/2) * se *c(-1,1)
  return(list(estimate = est, se = se, CI = CI))
}

compute_glm<- function(data, level) {
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
