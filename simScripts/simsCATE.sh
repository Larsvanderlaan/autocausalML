#!/bin/usr/env bash
filename="simsCATE_xgboost"
nsims=2500
export R_LIBS=~/Rlibs2
export R_LIBS_USER=~/Rlibs2
for n in 500 1000 2500 5000
do
  for const in 4 7 1
  do
    sbatch  --export=n=$n,const=$const ~/sieveSims/simsCATE.sbatch
  done
done
