#!/bin/sh
#
# This script executes gnuplot for each .csv file in ./csv subdirectory
#
for fnam in `ls ./csv/E*.csv | cut -f 3 -d "/" | cut -f 1 -d "."`
do
  echo "Plotting ${fnam}.png"
  gnuplot -e "filename='./csv/${fnam}.csv';fnam='${fnam}'" gnuplot.cfg > ${fnam}.png
done

#rm -f ./csv/*.csv
