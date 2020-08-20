#!/bin/sh
#
# This script splits data.txt into individual .csv files in ./csv subdirectory for each geometry 
#
rm -f tmp.txt
rm -rf ./csv
mkdir ./csv
rm -f fnam.tmp

echo "Extracting and formatting columns from data.txt"

grep Solved data.txt | cut -f 4,5,6,8,11,12,13,14,15 -d "," | sed 's/e\+/G+/g' | sed 's/e\-/G-/g'| sed 's/o\-//g' | sed 's/.://g' | sed 's/[a-z]//g' | sed 's/ //g' | sed 's/G+/e+/g' | sed 's/G-/e-/g' | sed 's/\/M/-M/g' | sed 's/M\//+M/g' > tmp.txt

echo "Generating filename list"

cat tmp.txt | cut -f 1-3 -d "," | sort | uniq > fnam.tmp

for f in `cat fnam.tmp`
do
  echo "Creating ${f}.csv"
  grep "${f}" tmp.txt >> ./csv/${f}.csv
done

rm -f tmp.txt
rm -f fnam.tmp
rm -f data.txt
