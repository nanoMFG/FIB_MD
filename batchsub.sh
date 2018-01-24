#!/bin/bash

JOB=`qsub jobsub`

for a in `seq 9`
do
  JOB=`qsub -W depend=afterany:$JOB jobsub`
done
