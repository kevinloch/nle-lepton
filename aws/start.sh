#!/bin/sh
#
# change localhost to url you will fetch the executable from
#
:>polylepton.cfg
while true
do
  curl -s -k https://localhost/lep/polylepton.cfg > polylepton.cfg.new
  if  [ `diff polylepton.cfg polylepton.cfg.new | grep -c .` -gt 0 ]
  then
    oldcmd=`grep cmd polylepton.cfg | cut -f 2 -d "="`
    mv polylepton.cfg.new polylepton.cfg
    cmd=`grep cmd polylepton.cfg | cut -f 2 -d "="`
    args=`grep args polylepton.cfg | cut -f 2 -d "="`
    killall ${oldcmd}
    rm -f ${oldcmd}
    curl -s -k https://localhost/lep/${cmd} > ${cmd}
    chmod 755 ./${cmd}
    numcpu=`cat /proc/cpuinfo | grep -c processor`
    for (( i=1; i<=${numcpu}; i++ ))
    do
      seed=$((${i} * 1000000))
      ./${cmd} ${seed} ${args} > polyleptonout-${i}.txt 2>&1 &
    done
  fi
  sleep 60
done
