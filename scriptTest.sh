#!/bin/bash
#
#SBATCH --job-name=Test_Wave
#SBATCH --mail-user=mail@student.ulg.ac.be
#SBATCH --mail-type=ALL
#SBATCH --output=TEST.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=100:00
#SBATCH --mem-per-cpu=400

Para="../Playgrounds/Wave_Para.kzr"
Geom="../Playgrounds/Wave_Geom.kzr"
TestName="Wave"


module load openmpi/1.6.4/gcc-4.9.2
module load cmake/3.5.2
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
mpirun sph  $Para $Geom $TestName
