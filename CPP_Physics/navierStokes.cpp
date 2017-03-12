#include "Main.h"
#include "Physics.h"

double continuity(int particleID, std::vector<int>& neighbors, std::vector<double>& kernelGradients,Field* currentField)
{
  double densityDerivative = 0.0;
  double scalarProduct = 0.0;
  for (int i = 0; i < neighbors.size(); i++){
    //Compute scalar product present in the formula
    for (int j = 0; j <= 2; j++){
      scalarProduct += (currentField->speed[3*particleID + j] - currentField->speed[3*neighbors[i] + j])
            * kernelGradients[3*i + j];
    }
    densityDerivative += currentField->mass[neighbors[i]]  * scalarProduct;
  }
  return densityDerivative;
}



void momentum(int particleID, std::vector<int>& neighbors, std::vector<double>& kernelGradients,Field* currentField , Parameter* parameter, std::vector<double>& speedDerivative)
{
  speedDerivative.assign(3,0.0);
  std::vector<double> viscosity;

  // ATTENTION !!! ON DOIT ENVOYER LA VRAIE VITESSE DU SON !!! (j'ai mis en attendant, 300m/s) OK pour vrai h
  viscosityComputation(particleID, neighbors, currentField, parameter, viscosity, 300.0, gethFromkh(parameter->kernel ,parameter->kh));

  for (int j = 0; j <= 2; j++){
    for (int i = 0; i < neighbors.size(); i++){
      speedDerivative[j] -= currentField->mass[neighbors[i]]
            * ( currentField->pressure[neighbors[i]] / ((currentField->density[neighbors[i]]*(currentField->density[neighbors[i]])))
                + currentField->pressure[particleID] / ((currentField->density[particleID]*(currentField->density[particleID])))
                + viscosity[i] )
            * kernelGradients[3*i + j];
    }
  }
  speedDerivative[2] -= parameter->g;
}
