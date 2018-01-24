#!/bin/bash

## This script will automatically submit mpich2.sh to the queue until you
## kill it.  To have it run in the background on the machine, either use
## nohup (i.e. nohup ./resubmit.sh &) or use screen.

function timer()
{
    if [[ $# -eq 0 ]]; then
        echo $(date '+%s')
    else
        local  stime=$1
        etime=$(date '+%s')

        if [[ -z "$stime" ]]; then stime=$etime; fi

        dt=$((etime - stime))
        ds=$((dt % 60))
        dm=$(((dt / 60) % 60))
        dh=$((dt / 3600))
        #printf '%d:%02d:%02d' $dh $dm $ds
        echo $dh
    fi
}



t=$(timer)

##first time

qsub -sync y mpich2.sh

while [ 1 ]
do

  t1=$(timer $t)
  t=$(timer)
  if [[ $t1 -lt 1 ]]; then
    echo $t1;
    exit
  fi

  qsub -sync y mpich2.sh

  #let's be nice.
  #sleep 10;

done
