
nsims <- 1
source("simsWagerHelper.R")
d <- 5
k <- 3
settings <- list("beta" = sim_beta, "pois" = sim_pois, "norm" = sim_norm, "log" = sim_log)
sim_generator <- settings[[setting]]
if(setting == "norm" {
  d <- 3
  k <- 2
})
n <- as.numeric(n)
nsims <- as.numeric(nsims)

results <- run_sims(n, d, k, sim_generator, nsims)
results <- as.data.frame(do.call(cbind, results))
fwrite(results,file = paste0("simsWager_" ,setting, "_", n, ".csv" ) )

