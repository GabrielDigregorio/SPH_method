#include "Main.h"
#include "Physics.h"

void viscosityComputation(int particleID, std::vector<int>& neighbors, Field* currentField, Parameter* parameter,std::vector<double>& viscosity)
{
  double alpha   = parameter->alpha;
  double beta    = parameter->beta;
  double epsilon = parameter->epsilon;
  double c = parameter->c;
  double h = gethFromkh(parameter->kernel ,parameter->kh);

  switch(parameter->viscosityModel){
  case violeauArtificial :

      // For each neighbor of the particleID
      /* FAUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUX !!!!!!!!!!
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
      */

    break;

    default :
      viscosity.assign(neighbors.size(),0.0);
      std::cout<< "\t \t Non existing Viscosity. Viscosity set to zero.\n";
  }

}
