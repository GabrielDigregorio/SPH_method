#ifndef SPH_H
#define SPH_H

#include <string>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <map>
#include <cmath>
#include <cassert>
#include <ctime>
#include <fstream>
#include <random>
#include <sstream>
#include <float.h>
#include <iomanip>


// Memory and CPU consumption
size_t GetMemory(bool screen, bool print);
size_t GetMemoryProcess(bool screen, bool print);
size_t GetMemoryProcessPeak(bool screen, bool print);

// Geometry.cpp
void meshcube(double o[3], double L[3], double s, std::vector<double> &pos, double perturbation = 0.0, bool stack = false);
void meshcylinder(double o[3], double L[3], double s, std::vector<double> &pos, double perturbation = 0.0, bool stack = false);
void meshsphere(double o[3], double L[3], double s, std::vector<double> &pos, double perturbation = 0.0, bool stack = false);

// Playground.cpp
//void fillVector(std::vector<double> &vect, double A, double B, double C);
//Playground ReadPlayground(const char *filename);
//void GeneratePlayground( std::vector<double> &posFree, std::vector<double> &posMoving, std::vector<double> &posFixed, const char *filename);

// Neighborhood.cpp
void neighborAllPair (std::vector<double> &pos,
                         double kh,
                         std::vector<double> &values,
                         std::vector<int> &row,
                         std::vector<int> &column);
void neighborLinkedList (std::vector<double> &pos,
                         double l[3],
                         double u[3],
                         double kh,
                         std::vector<double> &values,
                         std::vector<int> &row,
                         std::vector<int> &column);
void surroundingBoxes(int box, int nBoxesX, int nBoxesY, int nBoxesZ, std::vector<int> &surrBoxes);
double distance(std::vector<double> pos, int partA, int partB);

// Kernel.cpp
double Wab(std::vector<double> pos, int partA, int partB, double h, size_t choice);
double grad_Wab(std::vector<double> pos, int partA, int partB, double h, size_t choice);

// ParaView.cpp
void paraview(std::string const &filename,
              int step,
              std::vector<double> const &pos,
              std::map<std::string, std::vector<double> *> const &scalars,
              std::map<std::string, std::vector<double> *> const &vectors);



#endif
