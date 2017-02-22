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
#include<random>


// Geometry.cpp
void meshcube(double o[3], double L[3], double s, std::vector<double> &pos, double perturbation );
void meshcylinder(double o[3], double L, double R, double s, std::vector<double> &pos);

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


// ParaView.cpp
void paraview(std::string const &filename,
              int step,
              std::vector<double> const &pos,
              std::map<std::string, std::vector<double> *> const &scalars,
              std::map<std::string, std::vector<double> *> const &vectors);



#endif
