
#name <- "ComplexParametricHAL"
#name <- "LassoHighDim"
name <- "SimpleParametricHAL"
library(data.table)
consts <- c(3,5,8)
ns <- c(500,1000,2500,5000)

#name <- "LassoHighDim"
#consts <- c(1,3,5)
#ns <- c(250,500,1000,2000)
outs <- rbindlist(lapply(ns, function(n) {
  items <- lapply(consts, function(const) {
   try({ fread(paste0("./simScripts/", name, "_",const,"_", n, ".csv"))
  })})
  items <- items[sapply(items, is.data.frame)]
  rbindlist(items)
  })
)


 link <- name

ATE <- 1
 #ATE <- 1.811417



outs$const[outs$const == 3] <- "overlap: 1e-05"
outs$const[outs$const == 5] <- "overlap: 1e-08"
outs$const[outs$const == 8] <- "overlap: 1e-15"

outs$const[outs$const == 1] <- "overlap: 0.025"

outs$const[outs$const == 4] <- "overlap: 1e-06"
outs$const[outs$const == 7] <- "overlap: 1e-09"


out_sieveIF <- outs[,c("const", "n", "iter", "estimate", "CI_left", "CI_right")]
out_sieveIF$method = "Sieve - IF"
col_names <- colnames(out_sieveIF)

#out_sieve <- out_sieveIF_df
#out_sieve$method = "HAL-AdaMLE "

out_sieveIF_df <- outs[,c("const", "n", "iter", "estimate", "CI_df_left", "CI_df_right")]
out_sieveIF_df$method = "HAL-AdaMLE "
colnames(out_sieveIF_df) <- col_names

out_sieveboot <- outs[,c("const", "n", "iter", "estimate", "CI_boot_left", "CI_boot_right")]
out_sieveboot$method = "Sieve - bootstrap se"
colnames(out_sieveboot) <- col_names

out_TMLE <- outs[,c("const", "n", "iter", "estimate_tmle", "CI_left_tmle", "CI_right_tmle")]
out_TMLE$method = "TMLE"
colnames(out_TMLE) <- col_names

out_AIPW <- outs[,c("const", "n", "iter", "estimate_aipw", "CI_left_aipw", "CI_right_aipw")]
out_AIPW$method = "AIPW"
colnames(out_AIPW) <- col_names

out_lm <- outs[,c("const", "n", "iter", "estimate_lm", "CI_left_lm", "CI_right_lm")]
out_lm$method = "Linear-model"
colnames(out_lm) <- col_names

out_sieveIF_relaxed_df <- outs[,c("const", "n", "iter", "estimate_relaxed", "CI_df_left_relaxed", "CI_df_right_relaxed")]
out_sieveIF_relaxed_df$method = "HAL-AdaMLE (CV)"
colnames(out_sieveIF_relaxed_df) <- col_names

out_sieve_oracle <-out_sieveIF
out_sieve_oracle$method <- "Sieve - oracle se"

out_sieve_oracle[, CI_left:= estimate - qnorm(1-0.05/2) * sd(estimate), by = c("const", "n")]
out_sieve_oracle[, CI_right:= estimate + qnorm(1-0.05/2) * sd(estimate), by = c("const", "n")]


results <- rbind( out_sieve_oracle,  out_sieveIF_df, out_sieveIF_relaxed_df,   out_TMLE, out_AIPW, out_lm)


coverage <- results[, mean(ATE >= CI_left & ATE <= CI_right, na.rm=T), by = c("const", "method", "n")]
CI_width <- results[, mean(CI_right - CI_left, na.rm=T), by = c("const", "method", "n")]
cov_oracle = results[, mean(abs(estimate-ATE) <=  1.959964 * sd(estimate))  ,  by = c("const", "n", "method")]
bias <- results[, mean(estimate - ATE), by = c("const", "method", "n")]
se <- results[, sd(estimate - ATE, na.rm=T), by = c("const", "method", "n")]

plot_data <- data.table(const = coverage$const, method = coverage$method, n = coverage$n,  cov_oracle=cov_oracle[[4]],
                        bias = abs(bias[[4]]),
                        se = se[[4]],
                        coverage = coverage[[4]],
                        CI_width = CI_width[[4]]
)

library(ggplot2)

ggplot(plot_data[method %in% c("Linear-model", "AIPW", "TMLE", "HAL-AdaMLE ", "HAL-AdaMLE (CV)")], aes(x=n, y = bias, group = method, color = method, linetype=method)) + geom_line() +  facet_wrap(~ const) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5))  + labs(x = "Sample Size (n)", y = "Empirical Bias", color = "Method", group = "Method", linetype = "Method")





ggsave(paste0(link, "_bias.pdf"), width = 7, height = 4)

ggplot(plot_data[method %in% c("Linear-model", "AIPW", "TMLE", "HAL-AdaMLE ", "HAL-AdaMLE (CV)")], aes(x=n, y = se, group = method, color = method, linetype=method)) + geom_line() +  facet_wrap(~ const) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5))  + labs(x = "Sample Size (n)", y = "Empirical Standard Error", color = "Method", group = "Method", linetype = "Method")


#ggsave("SimPlotHALsmallSample_se.pdf", width = 7, height = 4)

ggsave(paste0(link, "_se.pdf"), width = 7, height = 4)



ggplot(plot_data[method %in% c("AIPW", "TMLE","Linear-model", "HAL-AdaMLE ", "HAL-AdaMLE (CV)")] , aes(x=n, y = cov_oracle, group = method, color = method, linetype=method)) + geom_line() +  facet_wrap(~ const) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5))  + labs(x = "Sample Size (n)", y = "Oracle CI coverage", color = "Method", group = "Method", linetype = "Method") + scale_y_continuous( )

ggsave(paste0(link, "_coverage_oracle.pdf"), width = 7, height = 4)


plot_data <- plot_data[(method %in%   c("AIPW", "TMLE","Linear-model", "HAL-AdaMLE ", "HAL-AdaMLE (CV)"))]


ggplot(plot_data , aes(x=n, y = coverage, group = method, color = method, linetype=method)) + geom_line() +  facet_wrap(~ const) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5))  + labs(x = "Sample Size (n)", y = "Empirical CI coverage", color = "Method", group = "Method", linetype = "Method") + scale_y_continuous(breaks = c(0.95, 0.93, 0.9, 0.8, 0.7))

ggsave(paste0(link, "_coverage.pdf"), width = 7, height = 4)


ggplot(plot_data  , aes(x=n, y = CI_width, group = method, color = method, linetype=method)) + geom_line() +  facet_wrap(~ const) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5)) + labs(x = "Sample Size (n)", y = "Average CI width", color = "Method", group = "Method", linetype = "Method")

ggsave(paste0(link, "_CI_width.pdf"), width = 7, height = 4)


