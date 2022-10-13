#!/bin/usr/env bash
filename="simsCATE_xgboost"
nsims=2500
export R_LIBS=~/Rlibs
export R_LIBS_USER=~/Rlibs
for n in 1000 #500 1000 2500 5000
do
  for const in 3 #5 8
  do
    sbatch  --export=n=$n,const=$const simsCATE_Complex.sbatch
  done
done
