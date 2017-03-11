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
                        std::vector<std::vector<double> > &kernelGradientsAll,
                        Kernel myKernel);
void neighborLinkedList(std::vector<double> &pos,
                        double l[3],
                        double u[3],
                        double kh,
                        std::vector<std::vector<int> > &neighborsAll,
                        std::vector<std::vector<double> > &kernelGradientsAll,
                        Kernel myKernel);
void surroundingBoxes(int box, int nBoxesX, int nBoxesY, int nBoxesZ, std::vector<int> &surrBoxes);
double distance(std::vector<double> &pos, int partA, int partB);
void findNeighbors(int particleID, std::vector<double> &pos, double kh,
                    std::vector<std::vector<int> > &boxes,
                    std::vector<int> &surrBoxes,
                    std::vector<int> &neighbors,
                    std::vector<double> &kernelGradients,
                    Kernel myKernel);
void findNeighbors(int particleID, std::vector<double> &pos, double kh,
                    std::vector<std::vector<int> > &boxes,
                    std::vector<int> &surrBoxes,
                    std::vector<int> &neighbors,
                    std::vector<double> &kernelGradients,
                    Kernel myKernel,
                    std::vector<double> &kernelGradientsSamples,
                    int resolution);
void sortParticles(std::vector<double> &pos, double l[3], double u[3], double kh,
                   std::vector<std::vector<int> > &boxes);
void boxMesh(double l[3], double u[3], double kh,
             std::vector<std::vector<int> > &boxes,
             std::vector<std::vector<int> > &surrBoxesAll);


// TimeIntegration.cpp
bool timeIntegration(Field* currentField, Field* nextField, Parameter* parameter, std::vector<std::vector<int> >& boxes,
std::vector<std::vector<int> >& surrBoxesAll, unsigned int n);


// Kernel.cpp
void kernelGradPre(Kernel myKernel, int resolution, double kh,
        std::vector<double> &kernelGradientsSamples);
int indexSamples(int resolution, double r, double kh);
double Wab(double r, double kh, Kernel choice);
double gradWab(double r, double kh, Kernel choice);

// Init.cpp
void densityInit(Field* field,Parameter* parameter);
void pressureComputation(Field* field,Parameter* parameter);
void massInit(Field* field,Parameter* parameter);

// updateMovingSpeed.cpp
void updateMovingSpeed(Field* field, Parameter* parameter, double t, MoveMod myMod);

// navierStokes.cpp
double continuity(int particleID, std::vector<int>& neighbors, std::vector<double>& kernelGradients,Field* currentField);
void momentum(int particleID, std::vector<int>& neighbors, std::vector<double>& kernelGradients,Field* currentField , Parameter* parameter,std::vector<double>& speedDerivative);

// navierStokes.cpp
double continuity(int particleID, std::vector<int>& neighbors, std::vector<double>& kernelGradients,Field* currentField);
void momentum(int particleID, std::vector<int>& neighbors, std::vector<double>& kernelGradients,Field* currentField , Parameter* parameter,std::vector<double>& speedDerivative);

// viscosityComputation.cpp
void viscosityComputation(int particleID, std::vector<int>& neighbors, Field* currentField, Parameter* parameter,std::vector<double>& viscosity, double c, double h, ViscoMod myViscoMod);

#endif
