#!/bin/bash
#$ -S /bin/bash
#
#$ -cwd
# Combine stdout and stderr into one file
#$ -j y
# Set the Parallel Environment and number of procs.
#$ -pe mpich2 64
# The queue, select shadowfax.q for the shadowfaxes
#$ -q smaug.q
# Job name
#$ -N si_ind35
# Hard runtime limit for this job, 1 hour
#$ -l h_rt=06:00:0
#$ -m e
#$ -M shubhro.kallol@gmail.com
#$ -o out.out

port=$((JOB_ID % 5000 + 20000))

#source /opt/intel/cce/10.1.018/bin/iccvars.sh
#source /opt/intel/fce/10.1.018/bin/ifortvars.sh

echo "Got $NSLOTS slots."

mpiexec -n $NSLOTS -machinefile $TMPDIR/machines -port $port -phrase asdfasdf mdrun2 >> log.out

exit 0

