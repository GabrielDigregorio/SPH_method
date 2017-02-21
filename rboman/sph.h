#ifndef SPH_H
#define SPH_H

#include <string>
#include <vector>
#include <iostream>
#include <map>
#include <cmath>
#include <cassert>

void meshcube(double o[3], double L[3], double s, std::vector<double> &pos);
void paraview(std::string const &filename, 
              int step,
              std::vector<double> const &pos,
              std::map<std::string, std::vector<double> *> const &scalars,
              std::map<std::string, std::vector<double> *> const &vectors);       



#endif

