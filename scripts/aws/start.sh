#!/bin/sh
#
# change localhost to url you will fetch the executable from
#
:>nle-lepton.cfg
while true
do
  curl -s -k https://localhost/lep/nle-lepton.cfg > nle-lepton.cfg.new
  if  [ `diff nle-lepton.cfg nle-lepton.cfg.new | grep -c .` -gt 0 ]
  then
    oldcmd=`grep cmd nle-lepton.cfg | cut -f 2 -d "="`
    mv nle-lepton.cfg.new nle-lepton.cfg
    cmd=`grep cmd nle-lepton.cfg | cut -f 2 -d "="`
    args=`grep args nle-lepton.cfg | cut -f 2 -d "="`
    killall ${oldcmd}
    rm -f ${oldcmd}
    curl -s -k https://localhost/lep/${cmd} > ${cmd}
    chmod 755 ./${cmd}
    numcpu=`cat /proc/cpuinfo | grep -c processor`
    for (( i=1; i<=${numcpu}; i++ ))
    do
      seed=$((${i} * 1000000))
      ./${cmd} ${seed} ${args} > leptonout-${i}.txt 2>&1 &
    done
  fi
  sleep 60
done
