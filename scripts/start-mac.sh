#!/bin/sh
numcpu=4

for (( i=1; i<=${numcpu}; i++ ))
do
  seed=$((${i} * 1000000))
  ./nle-lepton ${seed} 12 5 70 65 > leptonout-${i}.txt 2>&1 &
done
