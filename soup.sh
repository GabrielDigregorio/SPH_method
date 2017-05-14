#!/bin/bash
#
#SBATCH --job-name=soup
#SBATCH --mail-user=xavier.adriaens@student.ulg.ac.be
#SBATCH --mail-type=ALL
#SBATCH --output=soup.txt
#
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=16
#SBATCH --time=360:00
#SBATCH --mem-per-cpu=500

Para="../Playgrounds/soupPara.kzr"
Geom="../Playgrounds/soupGeom.kzr"
TestName="soup"

module load openmpi/1.6.4/gcc-4.9.2
module load cmake/3.5.2

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
mpirun sph  $Para $Geom $TestName
