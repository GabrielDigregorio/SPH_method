#!/bin/bash
#
#SBATCH --job-name=SPH_mpi
#SBATCH --output=SPH_Result_mpi.txt
#
#SBATCH --ntasks=1
#SBATCH --time=10:00
#SBATCH --mem-per-cpu=100

module load gcc
module load intel
module load cmake
module load slurm
module load openmpi
make
