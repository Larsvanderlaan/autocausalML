#setwd("./wagerScripts")
#n <- 1000
#setting <- "norm"

set.seed(17592057)
nsims <- 500
source("simsWagerHelper.R")
d <- 2
k <- 2
settings <- list("beta" = sim_beta, "pois" = sim_pois, "norm" = sim_norm, "log" = sim_log, "bin" = sim_bin)
sim_generator <- settings[[setting]]

n <- as.numeric(n)
nsims <- as.numeric(nsims)
doMC::registerDoMC(cores = 11)
test_data <- sim_generator(n ,d, k)
nknots <- log(n) * c(sqrt(n), sqrt(n))
nknots_tau <- sqrt(n)
# if(n== 500){
#   nknots <- c(50,50)
#   nknots_tau <- c(10)
# } else if(n == 1000) {
#   nknots <- c(100,100)
#   nknots_tau <- c(20)
# } else if(n == 1500) {
#   nknots <- c(75,75)
#   nknots_tau <- c(30)
# } else if(n == 2000) {
#   nknots <- c(100,100)
#   nknots_tau <- c(50)
# } else if(n == 3000) {
#   nknots <- c(200,200)
#   nknots_tau <- c(50)
# }

results <- run_sims(n, d, k, sim_generator, nsims, nknots, nknots_tau)
results <- as.data.frame(do.call(cbind, results))
fwrite(results,file = paste0("simsWager_" ,setting, "_", n, "d2k2_rootn.csv" ) )

