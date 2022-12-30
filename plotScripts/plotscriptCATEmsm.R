

name <- "SimsHALCATE"

library(data.table)
consts <- c(1,4,7)
ns <- c(500, 1000, 2500, 5000)
ns <- c(500 ,1000 ,2000 ,3000, 4000 ,5000)
outs <- rbindlist(lapply(ns, function(n) {
  items <- lapply(consts, function(const) {
    fread(paste0("./simScripts/", name, "_",const,"_", n, ".csv"))
  })

  rbindlist(items)
})
)


link <- name


ATE <- 1
library(data.table)

# outs$const[outs$const == 3] <- "overlap: 1e-06"
# outs$const[outs$const == 5] <- "overlap: 1e-10"
# outs$const[outs$const == 8] <- "overlap: 1e-15"
outs$const[outs$const == 1] <- "overlap: 0.025"

outs$const[outs$const == 4] <- "overlap: 1e-06"
outs$const[outs$const == 7] <- "overlap: 1e-09"


out_sieveIF <- outs[,c("name","const", "n", "iter", "estimate", "CI_left", "CI_right")]
out_sieveIF$method = "Sieve - IF"
col_names <- colnames(out_sieveIF)

out_sieve <- out_sieveIF
out_sieve$method = "Sieve - plugin"

out_sieveIF_df <- outs[,c("name","const", "n", "iter", "estimate", "CI_df_left", "CI_df_right")]
out_sieveIF_df$method = "HAL-AdaMLE "
colnames(out_sieveIF_df) <- col_names

out_sieveboot <- outs[,c("name","const", "n", "iter", "estimate", "CI_boot_left", "CI_boot_right")]
out_sieveboot$method = "Sieve - bootstrap se"
colnames(out_sieveboot) <- col_names

out_npglm <- outs[,c("name","const", "n", "iter", "estimate_npglm", "CI_left_npglm", "CI_right_npglm")]
out_npglm$method = "npWorkingTMLE"
colnames(out_npglm) <- col_names

out_spglm  <- outs[,c("name", "const", "n", "iter", "estimate_spglm", "CI_left_spglm", "CI_right_spglm")]
out_spglm$method = "spTMLE"
colnames(out_spglm) <- col_names





out_sieve_oracle <-out_sieveIF
out_sieve_oracle$method <- "Sieve - oracle se"

out_sieve_oracle[, CI_left:= estimate - qnorm(1-0.05/2) * sd(estimate), by = c("const", "n", "name")]
out_sieve_oracle[, CI_right:= estimate + qnorm(1-0.05/2) * sd(estimate), by = c("const", "n", "name")]


results <- rbind( out_sieve_oracle, out_sieveIF_df,  out_npglm, out_spglm)

results[-grep( "W",results$name ), "name"] <- "intercept"
results[grep( "W",results$name ), "name"] <- "X"

coverage <- results[, mean(ATE >= CI_left & ATE <= CI_right), by = c("const", "method", "n", "name")]
CI_width <- results[, mean(CI_right - CI_left), by = c("const", "method", "n", "name")]
cov_oracle = results[, mean(abs(estimate-ATE) <=  1.959964 * sd(estimate))  ,  by = c("const", "n", "method", "name")]
bias <- results[, mean(estimate - ATE), by = c("const", "method", "n", "name")]
se <- results[, sd(estimate - ATE), by = c("const", "method", "n", "name")]

plot_data <- data.table(name = coverage$name, const = coverage$const, method = coverage$method, n = coverage$n,  cov_oracle=cov_oracle[[5]],
                        bias = abs(bias[[5]]),
                        se = se[[5]],
                        coverage = coverage[[5]],
                        CI_width = CI_width[[5]]
)

library(ggplot2)
plot_data_orig <- plot_data

plot_data <- plot_data_orig

plot_data <- plot_data[plot_data$name %in% c("intercept")]
ggplot(plot_data[method %in% c("npWorkingTMLE", "spTMLE", "HAL-AdaMLE ")], aes(x=n, y = bias,   color = method, linetype=method)) + geom_line() +  facet_wrap(~ const) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5))  + labs(x = "Sample Size (n)", y = "Empirical Bias", color = "Method", group = "Method", linetype = "Method")



ggsave("SimPlotHALCATE_biasinter.pdf", width = 7, height = 4)


ggplot(plot_data[method %in%  c("npWorkingTMLE", "spTMLE", "HAL-AdaMLE ")], aes(x=n, y = se,  color = method, linetype=method)) + geom_line() +  facet_wrap(~ const) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5))  + labs(x = "Sample Size (n)", y = "Empirical Standard Error", color = "Method", group = "Method", linetype = "Method")

ggsave("SimPlotHALCATE_seinter.pdf", width = 7, height = 4)



ggplot(plot_data[method %in%  c("npWorkingTMLE", "spTMLE", "HAL-AdaMLE ")] , aes(x=n, y = cov_oracle, color = method, linetype=method)) + geom_line() +  facet_wrap(~ const) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5))  + labs(x = "Sample Size (n)", y = "Oracle CI coverage", color = "Method", group = "Method", linetype = "Method")

ggsave("SimPlotHALCATE_coverage_oracleinter.pdf", width = 7, height = 4)



ggplot(plot_data[method %in%  c("npWorkingTMLE", "spTMLE", "HAL-AdaMLE ")] , aes(x=n, y = coverage,  color = method, linetype=method)) + geom_line() +  facet_wrap(~ const) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5))  + labs(x = "Sample Size (n)", y = "Empirical CI coverage", color = "Method", group = "Method", linetype = "Method") + scale_y_continuous(breaks = c(0.95, 0.93, 0.9, 0.8, 0.7, 0.6))

ggsave("SimPlotHALCATE_coverageinter.pdf", width = 7, height = 4)



ggplot(plot_data[  method %in% c("npWorkingTMLE", "spTMLE", "HAL-AdaMLE ")] , aes(x=n, y = CI_width,  color = method, linetype=method)) + geom_line() +  facet_wrap(~ const) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5)) + labs(x = "Sample Size (n)", y = "Average CI width", color = "Method", group = "Method", linetype = "Method")

ggsave("SimPlotHALCATE_CI_widthinter.pdf", width = 7, height = 4)







plot_data <- plot_data_orig
plot_data <- plot_data[plot_data$name %in% c("X")]
ggplot(plot_data[method %in% c("npWorkingTMLE", "spTMLE", "HAL-AdaMLE ")], aes(x=n, y = bias,   color = method, linetype=method)) + geom_line() +  facet_wrap(~ const) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5))  + labs(x = "Sample Size (n)", y = "Empirical Bias", color = "Method", group = "Method", linetype = "Method")



ggsave("SimPlotHALCATE_bias_X.pdf", width = 7, height = 4)


ggplot(plot_data[method %in%  c("npWorkingTMLE", "spTMLE", "HAL-AdaMLE ")], aes(x=n, y = se,  color = method, linetype=method)) + geom_line() +  facet_wrap(~ const) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5))  + labs(x = "Sample Size (n)", y = "Empirical Standard Error", color = "Method", group = "Method", linetype = "Method")

ggsave("SimPlotHALCATE_se_X.pdf", width = 7, height = 4)



ggplot(plot_data[method %in%  c("npWorkingTMLE", "spTMLE", "HAL-AdaMLE ")] , aes(x=n, y = cov_oracle, color = method, linetype=method)) + geom_line() +  facet_wrap(~ const) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5))  + labs(x = "Sample Size (n)", y = "Oracle CI coverage", color = "Method", group = "Method", linetype = "Method")

ggsave("SimPlotHALCATE_coverage_oracle_X.pdf", width = 7, height = 4)



ggplot(plot_data[method %in%  c("npWorkingTMLE", "spTMLE", "HAL-AdaMLE ")] , aes(x=n, y = coverage,  color = method, linetype=method)) + geom_line() +  facet_wrap(~ const) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5))  + labs(x = "Sample Size (n)", y = "Empirical CI coverage", color = "Method", group = "Method", linetype = "Method") + scale_y_continuous(breaks = c(0.95,0.93, 0.9, 0.8, 0.7, 0.6))

ggsave("SimPlotHALCATE_coverage_X.pdf", width = 7, height = 4)



ggplot(plot_data[  method %in% c("npWorkingTMLE", "spTMLE", "HAL-AdaMLE ")] , aes(x=n, y = CI_width,  color = method, linetype=method)) + geom_line() +  facet_wrap(~ const) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5)) + labs(x = "Sample Size (n)", y = "Average CI width", color = "Method", group = "Method", linetype = "Method")

ggsave("SimPlotHALCATE_CI_width_X.pdf", width = 7, height = 4)




