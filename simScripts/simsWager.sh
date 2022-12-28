#!/bin/usr/env bash
filename="simsCATE_xgboost"
nsims=2500
export R_LIBS=~/Rlibs2
export R_LIBS_USER=~/Rlibs2
for n in 500 1000 1500 2000 3000
do
for setting in "beta" "norm" "log" "pois" "bin"
do
sbatch  --export=n=$n,setting=$setting ~/sieveSims/simsWager.sbatch
done
done
