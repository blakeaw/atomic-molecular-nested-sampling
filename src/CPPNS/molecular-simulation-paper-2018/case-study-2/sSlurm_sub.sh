#!/bin/bash

#SBATCH --job-name=run-%j.k1rep
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err


# run job on single node as a single task
# but request that parameter 1 cpus be assigned
# for this task

srun -n1 -l ./${1}
