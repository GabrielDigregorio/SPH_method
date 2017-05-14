///**************************************************************************
/// SOURCE: Function to implement continuity and momentum equations
///**************************************************************************
#include "Main.h"
#include "Physics.h"

/*
*Input:
*- particleID: ID of the particle for which equation is computed
*- neighbors: vector of the ID of the particles surrounding the particle for which the equation is computed
*- kernelGradients[i]: kernel gradients between the particle for which equation is computed and its neighbors
*- currentField: field containing information about all particles
*Output:
*- densityDerivative: derivative of the density of particleID
*/
double continuity(int particleID, std::vector<int>& neighbors, std::vector<double>& kernelGradients,Field* currentField)
{
  double densityDerivative = 0.0;
  double scalarProduct;
  for (int i = 0; i < neighbors.size(); i++)
  {
    scalarProduct = 0.0;
    //Compute scalar product present in the formula
    for (int j = 0; j <= 2; j++)
    {
      scalarProduct += (currentField->speed[j][particleID] - currentField->speed[j][neighbors[i]])
            * kernelGradients[3*i + j];
    }
    densityDerivative += currentField->mass[neighbors[i]]  * scalarProduct;
  }
  return densityDerivative;
}


/*
*Input:
*- particleID: ID of the particle for which equation is computed
*- neighbors: vector of the ID of the particles surrounding the particle for which the equation is computed
*- kernelGradients[i]: kernel gradients between the particle for which equation is computed and its neighbors
*- currentField: field containing information about all particles
*- speedDerivative: vector containing the derivative of speed in each direction for each free particle
*Decscription:
* Compute the speed derivative related to particleID and store it in the speedDerivative vector
*/
// We could make it such that gravity is a vector ?
void momentum(int particleID, std::vector<int>& neighbors, std::vector<double>& kernelGradients,Field* currentField , Parameter* parameter, std::vector<double>& speedDerivative, std::vector<double> &viscosity)
{
  viscosity.resize(neighbors.size());
  viscosityComputation(particleID, neighbors, currentField, parameter, viscosity);

  for (int j = 0; j <= 2; j++)
  {
    for (int i = 0; i < neighbors.size(); i++)
    {
      speedDerivative[3*particleID + j] -= currentField->mass[neighbors[i]]
            * ( currentField->pressure[neighbors[i]] / ((currentField->density[neighbors[i]]*(currentField->density[neighbors[i]])))
                + currentField->pressure[particleID] / ((currentField->density[particleID]*(currentField->density[particleID])))
                + viscosity[i] )
            * kernelGradients[3*i + j];
    }
  }
  //speedDerivative[3*particleID + 2] -= parameter->g; // Gravitational acceleration
    double dist2 = currentField->pos[0][particleID]*currentField->pos[0][particleID] + currentField->pos[1][particleID]*currentField->pos[1][particleID] + currentField->pos[2][particleID]*currentField->pos[2][particleID];
    double force = 10.0/dist2;
    if(force > 10.0)
        force = 10.0;

    speedDerivative[3*particleID + 0] -= currentField->pos[0][particleID] * force;
    speedDerivative[3*particleID + 1] -= currentField->pos[1][particleID] * force;
    speedDerivative[3*particleID + 2] -= currentField->pos[2][particleID] * force;

}

/*
*
*
*
*/
void xsphCorrection(int particleID, std::vector<int> &neighbors, std::vector<double>& kernelValues, Field* currentField, Parameter* parameter, std::vector<double>& positionDerivative){

    for (int j = 0; j <= 2; j++){
        double particleSpeed = currentField->speed[j][particleID];
        positionDerivative[3*particleID + j] = particleSpeed;

        for (int i = 0; i < neighbors.size(); i++){
            positionDerivative[3*particleID + j] += parameter->epsilonXSPH * (currentField->speed[j][neighbors[i]] - particleSpeed)
                * kernelValues[i] * currentField->mass[neighbors[i]] / currentField->density[neighbors[i]];
        }
    }



}
