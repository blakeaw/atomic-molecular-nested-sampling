#!/bin/bash

#SBATCH --job-name=Dimer.job
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#Print  detailed  event  logging to error file
#SBATCH -v
#SBATCH --share


# run job on single node as a single task
# but request that parameter 1 cpus be assigned
# for this task
cp1=${1}
cp2=${2}
cptotal=`echo "$((cp1 * cp2))" `
srun -s -n${cp1} -c${cp2} -l ./${3}
##srun -s -n ${cptotal} --hint=nomultithread -l ./${3}
