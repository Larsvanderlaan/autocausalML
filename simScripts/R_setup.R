.libPaths( c( .libPaths(), "~/Rlibs2") )
setwd("~/sieve")
print(getwd())
nsims = 5000
library(autocausalML)
#library(future)
#plan(multisession, workers = 16)

