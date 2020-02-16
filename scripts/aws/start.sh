#!/bin/sh
#
# change localhost to url you will fetch the executable from
#
cmd=nle-lepton
cfg=nle-lepton.cfg

:>${cfg}

while true
do
  curl -s -k https://localhost/lep/$cfg > ${cfg}.new
  if  [ `diff ${cfg} ${cfg}.new | grep -c .` -gt 0 ]
  then
    mv ${cfg}.new ${cfg}
    killall ${cmd}
    rm -f ${cmd}
    curl -s -k https://localhost/lep/${cmd} > ${cmd}
    chmod 755 ./${cmd}
    if [ `grep -c 'cloud_run=yes' ${cfg}` -eq 1 ]
    then
      numcpu=`cat /proc/cpuinfo | grep -c processor`
      for (( i=1; i<=${numcpu}; i++ ))
      do
        seed=$((${i} * 1000000))
        ./${cmd} -s ${seed} > leptonout-${i}.txt 2>&1 &
      done
    fi
  fi
  sleep 60
done
