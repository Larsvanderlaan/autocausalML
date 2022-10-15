#!/bin/usr/env bash
filename="simsCATE_xgboost"
nsims=2500
export R_LIBS=~/Rlibs2
export R_LIBS_USER=~/Rlibs2
for n in 50 100 150 200 300 500
do
  for const in 4 7 1
  do
    sbatch  --export=n=$n,const=$const ~/sieveSims/simsSmallSample.sbatch
  done
done
