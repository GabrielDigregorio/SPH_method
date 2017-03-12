#include "Main.h"
#include "Physics.h"

void viscosityComputation(int particleID, std::vector<int>& neighbors, Field* currentField, Parameter* parameter,std::vector<double>& viscosity, double c, double h)
{
  //viscosity.assign(neighbors.size(),0.0); // To be changed to fit formula of Goffin

  double alpha   = 1.0;
  double beta    = 0.0;
  double epsilon = 0.01; // See in litterature

  switch( parameter->viscoMod){
  case 1 :
      // For each neighbor of the particleID
      for (int i = 0; i < neighbors.size(); i++)
      {
          double U_ij = currentField->speed[neighbors[i]] * currentField->speed[particleID]
                + currentField->speed[neighbors[i]] * currentField->speed[particleID]
                + currentField->speed[neighbors[i]]*currentField->speed[particleID]; // velocity
          double R_ij = sqrt(distance(currentField->pos, particleID, neighbors[i])); // distance

          if(U_ij * R_ij < 0)
          {
            double rho  = 0.5 * (currentField->density[particleID] + currentField->density[neighbors[i]]); // average
            double nu2  = epsilon*h*h;
            double mu   = (h*U_ij*R_ij) / (R_ij*R_ij + nu2);

            viscosity[i] = ( -alpha*c*mu + beta*mu*mu ) / (rho);
          }
          else viscosity[i] = 0.0;
      }

    break;

    default :
      std::cout<< "Non existing Viscosity.\n";
  }
 
}
