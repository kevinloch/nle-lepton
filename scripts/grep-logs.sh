#!/bin/sh

#
# sample script to extract results from httpd logs
# run from cron as often as desired
#
# edit these paths

out=/var/www/status.html
bak=/var/www/status-noreload.html.bak
noreload=/var/www/status-noreload.html
logfile=/var/log/http/httpd-log



sleep 15  # avoid race with http log rotate script
:>${bak}
mv ${noreload} ${bak}

echo '<!DOCTYPE html>' >> ${noreload}
echo '<html>' >> ${noreload}
echo '<head><title>Latest lepton formula search data</title>' >> ${noreload}
echo '</head>' >> ${noreload}
echo '<body><pre>' >> ${noreload}
grep result ${logfile} | cut -f 2 -d '"' | cut -f 2 -d " " | cut -f 3-99 -d "/" | sed 's/_/ /g' | gsort -g -k2 -t,  >>${noreload}
echo '</pre></body>' >> ${noreload}
echo '</html>' >> ${noreload}

:>${out}
echo '<!DOCTYPE html>' >> ${out}
echo '<html>' >> ${out}
echo '<head><title>Latest lepton formula search data</title>' >> ${out}
echo '<META http-equiv="refresh" CONTENT="30">' >> ${out}
echo '</head>' >> ${out}
echo '<body><pre>' >> ${out}
grep result ${noreload} | head -10000 >> ${out}
echo '</pre></body>' >> ${out}
echo '</html>' >> ${out}
