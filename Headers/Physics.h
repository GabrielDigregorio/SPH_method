///**************************************************************************
/// HEADER: Function Used In SPH Method
///**************************************************************************

#ifndef PHYSICS_H
#define PHYSICS_H
#include "Structures.h"

// Geometry.cpp
#include <random>
void meshcube(double o[3], double L[3], double s, std::vector<double> &pos, double perturbation = 0.0, bool stack = false);
void meshcylinder(double o[3], double L[3], double s, std::vector<double> &pos, double perturbation = 0.0, bool stack = false);
void meshsphere(double o[3], double L[3], double s, std::vector<double> &pos, double perturbation = 0.0, bool stack = false);


// Neighborhood.cpp
void neighborAllPair (std::vector<double> &pos,
                         double kh,
                         std::vector<std::vector<int> > &neighborsAll,
                         std::vector<std::vector<double> > &kernelGradientsAll);
void neighborLinkedList(std::vector<double> &pos,
                          double l[3],
                          double u[3],
                          double kh,
                          std::vector<std::vector<int> > &neighborsAll,
                          std::vector<std::vector<double> > &kernelGradientsAll);
void neighborLinkedList (std::vector<double> &pos,
                         double l[3],
                         double u[3],
                         double kh,
                         std::vector<double> &values,
                         std::vector<int> &row,
                         std::vector<int> &column);
void surroundingBoxes(int box, int nBoxesX, int nBoxesY, int nBoxesZ, std::vector<int> &surrBoxes);
double distance(std::vector<double> &pos, int partA, int partB);
void findNeighbors(int particleID, std::vector<double> &pos, double kh2,
                   std::vector<std::vector<int> > &boxes,
                   std::vector<int> &surrBoxes,
                   std::vector<int> &neighbors,
                   std::vector<double> &kernelGradients);
void sortParticles(std::vector<double> &pos, double l[3], double u[3], double kh,
                   std::vector<std::vector<int> > &boxes);
void boxMesh(double l[3], double u[3], double kh,
             std::vector<std::vector<int> > &boxes,
             std::vector<std::vector<int> > &surrBoxesAll);
void boxClear(std::vector<std::vector<int> > &boxes);

// TimeIntegration.cpp
bool timeIntegration(Field* currentField, Field* nextField, Parameter* parameter, unsigned int n);


// Kernel.cpp
double Wab(std::vector<double> pos, int partA, int partB, double h, size_t choice);
double grad_Wab(std::vector<double> pos, int partA, int partB, double h, size_t choice);

// Init.cpp
void densityInit(Field* field,Parameter* parameter);
void pressureComputation(Field* field,Parameter* parameter);
void massInit(Field* field,Parameter* parameter);

// updateMovingSpeed.cpp
void updateMovingSpeed(Field* field,Parameter* parameter,double t,Mode_move mymode);

#endif
