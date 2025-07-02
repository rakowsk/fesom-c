#!/bin/bash

#SBATCH --account=kuesten.kuesten          # Your account
#SBATCH --time 3:00:00
#SBATCH -p smp
#SBATCH -N 1
#SBATCH --tasks-per-node 1
#SBATCH --cpus-per-task 64
#SBATCH --hint=nomultithread
#SBATCH --job-name=fsmc_SNS
#SBATCH --output=out_%x.%j
#SBATCH --qos=12h

module purge
module load intel.compiler
module load intel.mpi
module load netcdf/4.6.0_gnu
module list

module purge
module load  intel-oneapi-mpi/2021.6.0 netcdf-fortran/4.5.4-oneapi2022.1.0

EXE=../../src/fesom_coast.exe
export OMP_NUM_THREADS=1

ulimit -s unlimited

date
srun ${EXE}
date
