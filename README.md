Smoothed-particle hydrodynamics
MATH0471 – Spring 2017
v.1 (7/02/2017)
----------------------------------------------------------------------------------------
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

