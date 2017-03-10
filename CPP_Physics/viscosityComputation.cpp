#include "Main.h"
#include "Physics.h"

void viscosityComputation(int particleID, std::vector<int>& neighbors, Field* currentField, Parameter* parameter,std::vector<double>& viscosity)
{
  //viscosity.assign(neighbors.size(),0.0); // To be changed to fit formula of Goffin

  double alpha= 2;
  double beta = 0;
  double epsilon = 0.01; // See in litterature
  double c    = ????; // average speed of sound

  // For each neighbor of the particleID
  for (int i = 0; i < neighbors.size(); i++)
  {
      U_ij = currentField->speed[neighbors[i]] * currentField->speed[particleID]
             + currentField->speed[neighbors[i]] * currentField->speed[particleID]
             + currentField->speed[neighbors[i]]*currentField->speed[particleID]; // velocity
      R_ij = sqrt(distance(currentField->pos, int particleID, int neighbors[i])); // distance

      if(U_ij * R_ij < 0)
      {
        double rho  = 0.5 * (currentField->density[particleID] + currentField->density[neighbors[i]]); // average
        double nu2  = epsilon*h*h;
        double mu   = (h*U_ij*R_ij) / (R_ij*R_ij + nu2);

        viscosity[i] = ( -alpha*c*mu + beta*mu*mu ) / (rho);
      }
      else viscosity[i] = 0.0;
  }

}
