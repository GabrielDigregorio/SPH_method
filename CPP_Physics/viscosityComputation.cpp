#include "Main.h"
#include "Physics.h"

void viscosityComputation(int particleID, std::vector<int>& neighbors, Field* currentField, Parameter* parameter,std::vector<double>& viscosity)
{
  viscosity.assign(neighbors.size(),0.0); // To be changed to fit formula of Goffin
}
