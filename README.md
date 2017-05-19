# Smoothed-particle hydrodynamics
MATH0471 – Spring 2017
v.2 (19/05/2017)



## Assignment

This project consists in studying, implementing and validating a numerical scheme for the
solution of Navier-Stokes equations using the “smoothed-particle hydrodynamics” (SPH)
computational method. SPH originated in the late 1970’s for astrophysical problems, and
has been used since then in numerous application areas. The method is a mesh-free,
particle-based Lagrangian method, where the coordinates move with the fluid (particles).
We ask you to study, implement and test the SPH method presented in the following reference:
Louis Goffin, “Development of a didactic SPH model”, Travail de fin d’études réalisé
en vue de l’obtention du grade de Master Ingénieur Civil des Constructions, Université de
Liège, année académique 2012-2013.

The project is organized with 4 intermediate deadlines, and with students divided into 2
groups (group A and group B). For each deadline, one 8-page progress report (mandatory,
but not graded, with the first section for group A and the second for group B) is due: section
1 should detail the computer implementation; section 2 should detail the mathematical,
numerical and physical experiments:

1. Reading of L. Goffin’s master thesis (groups A and B); literature review of SPH
and particle search methods (group A); implementation of the kernel and of the
particle search using the linked list method described in Section 3.2.4 of L. Goffin’s
thesis (group B); validation and performance testing the search algorithm (group A);
design of an input file format for SPH simulations based on cubes (group A).
Deadline: February 28th (section 1 by group B; section 2 by group A).

2. Implementation of sequential 3D SPH for the Navier-Stokes formulation described in
L. Goffin’s master thesis using an Explicit Euler and a two-step Runge Kutta scheme
(group A); validation and performance testing of the code on the falling water cube
test-case from Section 5.1.1 of L. Goffin’s thesis (group B).
Deadline: March 21th (section 1 by group A; section 2 by group B).

3. Parallelization of the SPH code using OpenMP (group B); application on the dam
break problem (section 5.1.3 of L. Goffin’s thesis) (group A).
Deadline: April 4th (section 1 by group B; section 2 by group A).

4. Parallelization of the SPH using MPI (group A); application on the tsumami problem
from INFO0939 (group B).
Deadline: April 25th (section 1 by group A; section 2 by group B).

The full C/C++ code (in a single ZIP archive, directly configurable and compilable on the
NIC4 CECI cluster) should be sent for each deadline to both cgeuzaine@ulg.ac.be and
r.boman@ulg.ac.be.

The final report (single report of 60 pages that should present the method and numerical
results, the computer implementation and a detailed analysis of physical experiments on
non-trivial configurations) is due on May 19th. An oral presentation of the main project
results will be organized during the June exam session; individual theoretical and practical
questions will be asked to each member of the two student groups.


## Compilation

* Windows

A compiler (g++ MinGw) must be installed or an IDE with a build-in compiler can be used (e.g. Visual Studio). Also, the Cmake software must be installed. Compilation on Windows can be different, depending on the used IDE. Here follows some command lines for compilation:
```
mkdir build
cd ./build/
set PATH=%PATH:D:/Programy/Git/bin=%
cmake -G "MinGW Makefiles" ..
mingw32-make
```

* Linux & mac OS

The Cmake software must be installed and a g++ compiler must be available.
```
mkdir build
cd build
cmake ..
make
```

* Clusters

The code can be launched on clusters to benefit from the MPI implementation of the code and speed up the computation time of the simulations. All simulations in this project, that was composed of a large number of particles, were launched on Nic4 (a cluster hosted at the University of Liège (SEGI facility)). It features 128 computer nodes with two 8-cores Intel E5-2650 processors at 2.0 GHz and 64 GB of RAM (4GB/core). Some modules should be loaded before compiling:

```
module load slurm
module load openmpi
module load cmake
module load gcc

mkdir build
cd build
cmake ..
make
```

The following command could also be needed to set the g++ compiler by default:

```
export CXX=g++
```

## Launch an experiment 

* Create a new directory to store data


```
cd build
mkdir Results
mkdir ./Results/<myName>
```
The field "myName" should be replaced by an user chosen name for the simulation.

* Launch a new experiment (command line)


```
mpirun -np <nbProc> ./sph <pathToParameterFile> <pathToGeometryFile> <name>
```
The fields "nbProc" and "name" should be replaced by the desired number of processors and the desired name for the output files while "pathToParameterFile" and "pathToGeometryFile" should be replaced by the path to the parameter and the geometry file.


* Launch a new experiment (bash script)

An example file of a bash script is given here below

```
           #!/bin/bash
            #
            #SBATCH --job-name=wetFloor
            #SBATCH --mail-user=karl.ajumide@student.ulg.ac.be
            #SBATCH --mail-type=ALL
            #SBATCH --output=wetFloor.txt
            #
            #SBATCH --ntasks=8
            #SBATCH --cpus-per-task=4
            #SBATCH --time=100:00
            #SBATCH --mem-per-cpu=400
            
            Para="../Playgrounds/wetFloor.para"
            Geom="../Playgrounds/wetFloor.geom"
            TestName="wetFloor"
            
            
            module load openmpi/1.6.4/gcc-4.9.2 
            module load cmake/3.5.2 
            export OMP_NUM_THREADS=\$SLURM_CPUS_PER_TASK 
            mpirun sph \$Para \$Geom \$TestName 
```

The number of processors allocated to the job are fixed with "ntasks". The number of threads allocated to the job are specified with "cpus-per-task". The name of the input files as well as the output file name also have to be entered. After saving the file, it can be launched with

```
sbatch scriptTest.sh
```
            
            
## Experiment Analysis (Matlab)
A Matlab code was build in order to post-process the results. A basic "user interface" is available to ease the post-processing. This code was used to generate some graphs available in the report.
* Open a Matlab instance and go into the directory "Matlab" of the SPH project;
* Launch the GUI script in the command windows;
* Select the operation to carry out. Pay attention to the order: launching an experiment (n°2) cannot be done if the code is not compiled. Also, before analyzing an experiment (n°4), the results of this experiment has to already been computed and stored in the folder "build/Results".
* Enjoy... ;)

When the analysis is done (n°4) for a given experiment, all data are stored in the same experiment folder under a structure (.mat). This structure contains all the results computed during the analysis and can be used to plot graphs. In other words, this \matlab structure can be considered as a summary of all the output files generated by a simulation. 
