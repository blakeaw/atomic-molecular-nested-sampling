#!/bin/bash

#SBATCH --job-name=OMP-%j.dimers
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
##SBATCH --nodes=1
##SBATCH --ntasks-per-node=4
#Print  detailed  event  logging to error file
##SBATCH -v


#echo "Running slurm batch $2 with $1 tasks" 
##srun -n1 -c$1 -l ./$2
#srun -n$1 -l ./$2
# run job on single node as a single task
# but request that parameter 1 cpus be assigned
# for this task
cp1=${1}
cptotal=`echo "$((cp1 * 2))" `
srun -N1 -n1 -c${cptotal} -l ./${2}
