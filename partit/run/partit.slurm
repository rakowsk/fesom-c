#!/bin/bash

#SBATCH --account=kuesten.kuesten          # Your account
#SBATCH --time 02:30:00
#SBATCH -p smp
#SBATCH --tasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --mem=50G
#SBATCH --job-name=partit
#SBATCH --output=out_%x.%j
#SBATCH --qos=12h

# disable hyperthreading
#SBATCH --hint=nomultithread

# module load    xthi/1.0-openmpi4.1.3-gcc8.5.0   openmpi/4.1.3

## Uncomment the following line to enlarge the stacksize if needed,
##  e.g., if your code crashes with a spurious segmentation fault.
# ulimit -s unlimited

# This binds each thread to one core
export OMP_PROC_BIND=TRUE

# OpenMP and srun, both need to know the number of CPUs per task
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

srun ../bin/fvom_ini.x

