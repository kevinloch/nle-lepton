#!/bin/sh
#
# This master script runs sub-scripts to get data files and generate graphs of two-term test vs reference mass
#

#cd ../run-aws
#./rsync.sh
#cd ../plot

./get.sh
./split.sh
./plot.sh
