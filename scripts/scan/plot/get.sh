#!/bin/sh
#
# This script gets raw data files from multiple sources and combines them into data.txt
#

rm -f data.tmp
rm -f data.txt

# grab from output files in parent directory
cat ../leptonout*.txt >> data.tmp

# grab from output files in run* subdirectories of parent directory
for f in `ls .. | grep run`
do
  if [ -d ../${f} ]
  then
    cat ../${f}/leptonout*.txt >> data.tmp
  fi
done
 
cat data.tmp | sort | uniq >> data.txt
rm -f data.tmp
