#!/bin/bash

JOB=`sbatch job.mpi`

for a in `seq 5`
do
  JOB=`sbatch --dependency=afterany:$JOB job.mpi`
done
