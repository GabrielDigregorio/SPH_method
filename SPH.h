#ifndef SPH_H
#define SPH_H

#include <string>
#include <vector>
#include <iostream>
#include <map>
#include <cmath>
#include <cassert>

// Geometry.cpp
void meshcube(double o[3], double L[3], double s, std::vector<double> &pos);
void meshcylinder(double o[3], double L, double R, double s, std::vector<double> &pos);

// Neighborhood
void neighborNaive (double p[3], double h, std::vector<double> &pos, std::vector<double> &neighbor);

// ParaView.cpp
void paraview(std::string const &filename, 
              int step,
              std::vector<double> const &pos,
              std::map<std::string, std::vector<double> *> const &scalars,
              std::map<std::string, std::vector<double> *> const &vectors);       



#endif

