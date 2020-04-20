#!/bin/sh

#
# sample script to extract results from httpd logs
# run from cron as often as desired
#
# edit these paths

results=/var/www/results.html
resultsbak=/var/www/results-noreload.html.bak
resultstmp=/var/www/resultstmp.txt
resultsnoreload=/var/www/results-noreload.html
twoterm=/var/www/twoterm.html
twotermbak=/var/www/twoterm.html.bak
twotermnoreload=/var/www/twoterm-noreload.html
logfile=/var/log/http/httpd-log

# extract results
# back up noreload version
sleep 15  # avoid race with http log rotate script
rm -f ${resultsbak}
mv ${resultsnoreload} ${resultsbak}
rm -f ${twotermbak}
mv ${twotermnoreload} ${twotermbak}
# print header
echo '<!DOCTYPE html>' > ${resultsnoreload}
echo '<html>' >> ${resultsnoreload}
echo '<head><title>Latest lepton polyform search data</title>' >> ${resultsnoreload}
echo '</head>' >> ${resultsnoreload}
echo '<body><pre>' >> ${resultsnoreload}
# get results from http log file
rm -f ${resultstmp}
grep "result," ${logfile} | cut -f 2 -d '"' | cut -f 2 -d " " | cut -f 3-99 -d "/" | gsort -g -k2 -t,  > ${resultstmp}
# dedup and output
lastline="_"
for l in `cat ${resultstmp}`
do
  thisline=`echo ${l} | cut -f 8 -d ","`
  if [ ${thisline} != ${lastline} ]
  then
    echo "${l}" | sed 's/_/ /g' >> ${resultsnoreload}
    lastline=${thisline}
  fi
done
echo '</pre></body>' >> ${resultsnoreload}
echo '</html>' >> ${resultsnoreload}

# create truncated version with auto reload header
:>${results}
echo '<!DOCTYPE html>' > ${results}
echo '<html>' >> ${results}
echo '<head><title>Latest lepton polyform search data</title>' >> ${results}
echo '<META http-equiv="refresh" CONTENT="30">' >> ${results}
echo '</head>' >> ${results}
echo '<body><pre>' >> ${results}
grep "result," ${resultsnoreload} | head -10000 >> ${results}
echo '</pre></body>' >> ${results}
echo '</html>' >> ${results}

# extreact interesting two_term_test
echo '<!DOCTYPE html>' > ${twotermnoreload}
echo '<html>' >> ${twotermnoreload}
echo '<head><title>Latest lepton polyform search data</title>' >> ${twotermnoreload}
echo '</head>' >> ${twotermnoreload}
echo '<body><pre>' >> ${twotermnoreload}
grep "two-term_test," ${logfile} | cut -f 2 -d '"' | cut -f 2 -d " " | cut -f 3-99 -d "/" | sed 's/_/ /g' | gsort -n  >> ${twotermnoreload}
echo '</pre></body>' >> ${twotermnoreload}
echo '</html>' >> ${twotermnoreload}

# create truncated version with auto reload header
:>${twoterm}
echo '<!DOCTYPE html>' > ${twoterm}
echo '<html>' >> ${twoterm}
echo '<head><title>Latest lepton polyform search data</title>' >> ${twoterm}
echo '<META http-equiv="refresh" CONTENT="30">' >> ${twoterm}
echo '</head>' >> ${twoterm}
echo '<body><pre>' >> ${twoterm}
grep "two-term test," ${twotermnoreload} | head -10000 >> ${twoterm}
echo '</pre></body>' >> ${twoterm}
echo '</html>' >> ${twoterm}

