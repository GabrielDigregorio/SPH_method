///**************************************************************************
/// HEADER: Function Used In SPH Method
///**************************************************************************

#ifndef PHYSICS_H
#define PHYSICS_H
#include "Structures.h"

// Geometry.cpp
#include <random>
void RotateVector(std::vector<double> &pos,double teta[3],int i );
void meshcube(double o[3], double L[3],double teta[3], double s, std::vector<double> &pos, int* nPart, double* volPart, double perturbation = 0.0, bool stack = false);
void meshcylinder(double o[3], double L[3], double s, std::vector<double> &pos, int* nPart, double* volPart, double perturbation = 0.0, bool stack = false);
void meshsphere(double o[3], double L[3], double s, std::vector<double> &pos, int* nPart, double* volPart, double perturbation = 0.0, bool stack = false);
Error meshBathymetry(char* batFile, int numberGroundParticles, double height0, double hFreeSurface, double s, std::vector<double> &posFree,std::vector<double> &posFixed,  int* nPartFree, int* nPartFixed, double* volPart,
     double perturbation, bool stack);



// Neighborhood.cpp
void neighborAllPair (std::vector<double> (&pos)[3],
                        double kh,
                        std::vector<std::vector<int> > &neighborsAll,
                        std::vector<std::vector<double> > &kernelGradientsAll,
                        Kernel myKernel);
void neighborLinkedList(std::vector<double> (&pos)[3],
                        double l[3],
                        double u[3],
                        double kh,
                        std::vector<std::vector<int> > &neighborsAll,
                        std::vector<std::vector<double> > &kernelGradientsAll,
                        Kernel myKernel);
void surroundingBoxes(int box, int nBoxesX, int nBoxesY, int nBoxesZ, std::vector<int> &surrBoxes);
double distance(std::vector<double> (&pos)[3], int partA, int partB);
void findNeighbors(int particleID, std::vector<double> (&pos)[3], double kh,
                    std::vector<std::vector<int> > &boxes,
                    std::vector<int> &surrBoxes,
                    std::vector<int> &neighbors,
                    std::vector<double> &kernelGradients,
                    Kernel myKernel);
void findNeighbors(int particleID, std::vector<double> (&pos)[3], double kh,
                    std::vector<std::vector<int> > &boxes,
                    std::vector<int> &surrBoxes,
                    std::vector<int> &neighbors,
                    std::vector<double> &kernelGradients,
                    Kernel myKernel,
                    std::vector<double> &kernelGradientsSamples,
                    int resolution); // USELESS
void sortParticles(std::vector<double> (&pos)[3], double l[3], double u[3], double kh,
                   std::vector<std::vector<int> > &boxes);
void sortParticles(std::vector<double> (&pos)[3], double l[3], double u[3], double kh,
                  std::vector<std::vector<int> > &boxes, bool toOptimize);
void boxMesh(double l[3], double u[3], double kh,
             std::vector<std::vector<int> > &boxes,
             std::vector<std::vector<int> > &surrBoxesAll);
double boxSizeCalc(double kh, IntegrationMethod method);

// TimeIntegration.cpp
void timeIntegration(Field* currentField, Field* nextField, Parameter* parameter, SubdomainInfo &subdomainInfo,
    std::vector<std::vector<int> >& boxes, std::vector<std::vector<int> >& surrBoxesAll,
    double t, double k, std::vector<double> &timeInfo);


// Kernel.cpp
void kernelGradPre(Kernel myKernel, int resolution, double kh,
        std::vector<double> &kernelGradientsSamples);
int indexSamples(int resolution, double r, double kh);
double Wab(double r, double kh, Kernel choice);
double gradWab(double r, double kh, Kernel choice);
double gethFromkh(Kernel kernelType, double kh);

// Init.cpp
void speedInit(Field* field,Parameter* parameter);
void densityInit(Field* field,Parameter* parameter);
void pressureInit(Field* field,Parameter* parameter);
void pressureComputation(Field* field,Parameter* parameter);
void massInit(Field* field,Parameter* parameter,std::vector<double> &vol);

// updateMovingSpeed.cpp
void updateMovingSpeed(Field* field, Parameter* parameter, double t,int IDmovingboundary,int i );

// navierStokes.cpp
double continuity(int particleID, std::vector<int>& neighbors, std::vector<double>& kernelGradients,Field* currentField);
void momentum(int particleID, std::vector<int>& neighbors, std::vector<double>& kernelGradients,Field* currentField , Parameter* parameter,std::vector<double>& speedDerivative);

// navierStokes.cpp
double continuity(int particleID, std::vector<int>& neighbors, std::vector<double>& kernelGradients,Field* currentField);
void momentum(int particleID, std::vector<int>& neighbors, std::vector<double>& kernelGradients,Field* currentField , Parameter* parameter,std::vector<double>& speedDerivative);

// viscosityComputation.cpp
void viscosityComputation(int particleID, std::vector<int>& neighbors, Field* currentField, Parameter* parameter,std::vector<double>& viscosity);

// MPI.cpp
void scatterField(Field* globalField, Field* currentField, Parameter* parameter,
    SubdomainInfo &subdomainInfo);
void gatherField(Field* globalField, Field* localField, SubdomainInfo &subdomainInfo);
void processUpdate(Field* currentField);
int getDomainNumber(double x, std::vector<double> &limits, int nTasks);
void computeDomainIndex(std::vector<double> &posX,
    std::vector<double> &limits, std::vector<int> &nbPartNode,
    std::vector< std::pair<int,int> > &index, int nTasks);
void processUpdate(Field& localField, SubdomainInfo& subdomainInfo);
void resizeField(Field& field, int nMigrate);
void computeMigrateIndex(std::vector<double>& posX,
    std::vector< std::pair<int,int> >& index, int* nMigrate,
    double Xmin, double Xmax);
void computeOverlapIndex(std::vector<double>& posX,
    std::vector< std::pair<int,int> >& index, int* nOverlap,
    double leftMinX, double leftMaxX, double rightMinX, double rightMaxX);
void sortParticles(Field& field, std::vector< std::pair<int,int> >& index);
void resizeField(Field& field, int nMigrate);
void shareRKMidpoint(Field& field, SubdomainInfo &subdomainInfo);
void shareOverlap(Field& field, SubdomainInfo &subdomainInfo);
void deleteHalos(Field &field, SubdomainInfo &subdomainInfo);
void timeStepUpdate(double &nextK, double &localProposition, SubdomainInfo &subdomainInfo);

#endif
