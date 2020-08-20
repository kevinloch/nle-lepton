#!/bin/sh
#
# sample script to start multiple threads on linux
#

numcpu=`cat /proc/cpuinfo | grep -c processor`
for (( i=1; i<=${numcpu}; i++ ))
do
  seed=$((${i} * 1000000))
  ./nle-lepton -s ${seed} > leptonout-${i}.txt 2>&1 &
done

