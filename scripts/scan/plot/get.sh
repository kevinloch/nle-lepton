#!/bin/sh
#
# This script gets raw data files from multiple sources and combines them into data.txt
#

rm -f data.txt

echo "Finding raw data files..."

# Grab from output files in parent directory
for f in `ls .. | grep leptonout-`
do
  if [ -f ../${f} ]
  then
    echo "Found ../${f}"
    cat ../${f} >> data.txt
  fi
done

# Grab from output files in run* subdirectories of parent directory
for d in `ls .. | grep run`
do
  if [ -d ../${d} ]
  then
    echo "Found directory ../${d}"
    for f in `ls ../${d} | grep leptonout-`
    do
      if [ -f ../${d}/${f} ]
      then
        echo "Found ../${d}/${f}"
        cat ../${d}/${f} >> data.txt
      fi
    done
  fi
done
