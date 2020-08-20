#!/bin/sh
#
# sample script to start multiple threads on OSX
#

# edit this for the number of threads you want to start
numcpu=4

for (( i=1; i<=${numcpu}; i++ ))
do
  seed=$((${i} * 1000000))
  ./nle-lepton -s ${seed} > leptonout-${i}.txt 2>&1 &
done
