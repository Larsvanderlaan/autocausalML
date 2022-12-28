#setwd("./wagerScripts")
#n <- 1000
#setting <- "pois"


nsims <- 500
source("simsWagerHelper.R")
d <- 2
k <- 2
settings <- list("beta" = sim_beta, "pois" = sim_pois, "norm" = sim_norm, "log" = sim_log, "bin" = sim_bin)
sim_generator <- settings[[setting]]
if(setting == "norm") {
  d <- 2
  k <- 2
}
n <- as.numeric(n)
nsims <- as.numeric(nsims)
doMC::registerDoMC(cores = 11)
results <- run_sims(n, d, k, sim_generator, nsims)
results <- as.data.frame(do.call(cbind, results))
fwrite(results,file = paste0("simsWager_" ,setting, "_", n, "d2k2.csv" ) )

