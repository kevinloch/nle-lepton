#!/bin/sh
#
# change localhost to url you will fetch the executable from
#
:>lepton.cfg
while true
do
  curl -s -k https://localhost/lep/lepton.cfg > lepton.cfg.new
  if  [ `diff lepton.cfg lepton.cfg.new | grep -c .` -gt 0 ]
  then
    oldcmd=`grep cmd lepton.cfg | cut -f 2 -d "="`
    mv lepton.cfg.new lepton.cfg
    cmd=`grep cmd lepton.cfg | cut -f 2 -d "="`
    args=`grep args lepton.cfg | cut -f 2 -d "="`
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
