#!/bin/usr/env bash
filename="simsCATE_xgboost"
nsims=2500
export R_LIBS=~/Rlibs
export R_LIBS_USER=~/Rlibs
for n in 250 500 1000 2000
do
  for const in 1 3 5
  do
    sbatch  --export=n=$n,const=$const ~/sieveSims/simsHighDim.sbatch
  done
done
