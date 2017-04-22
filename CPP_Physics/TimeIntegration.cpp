///**************************************************************************
/// SOURCE: Function to integrate one time step.
///**************************************************************************
#include "Main.h"
#include "Physics.h"
#include "Tools.h"


/*
*Input:
*- currentField: field that contains all the information at time t
*- nextField: field in which results of step n are stored
*- parameter: pointer to the field containing the user defined parameters
*- currentDensityDerivative: vector containing derivative of density for each particle at time t
*- currentSpeedDerivative: vector containing derivative of velocity for each particle at time t
*- t: current simulation time
*- k: timestep
*Decscription:
* Knowing the field at time t (currentField) and the density and velocity derivative, computes the field at time t+k with euler integration method and store it in structure nextField
*/
void eulerUpdate(Field* currentField, Field* nextField,Parameter* parameter, SubdomainInfo &subdomainInfo, std::vector<double>& currentDensityDerivative, std::vector<double>& currentSpeedDerivative, double t, double k)
{
    // Loop on all the particles
    for(int i=subdomainInfo.startingParticle ; i<=subdomainInfo.endingParticle ; i++){
        switch (currentField->type[i]){
            // Free particles update
            case freePart:
            nextField->density[i] = currentField->density[i] + k*currentDensityDerivative[i];
            for (int j = 0; j <= 2; j++){
                nextField->speed[j][i] = currentField->speed[j][i] + k*currentSpeedDerivative[3*i + j];
                //std::cout << currentSpeedDerivative[3*i + j] << " ";
                nextField->pos[j][i] = currentField->pos[j][i] + k*currentField->speed[j][i];
            }
            // Fixed particles update
            case fixedPart:
            nextField->density[i] = currentField->density[i] + k*currentDensityDerivative[i];
            // Moving boundary particles update
            default:
            nextField->density[i] = currentField->density[i] + k*currentDensityDerivative[i];
            for (int j = 0; j <= 2; j++){
                nextField->pos[j][i] = currentField->pos[j][i] + k*currentField->speed[j][i];
            }
            //case movingPart:
            //updateMovingSpeed(nextField,parameter,t+k);//,i);
            // comment faire un switch avec n=2,3,4... indétermimé?
        }
        if(currentField->type[i]>=2)// then we have a moving boundary
        {
          int IDmovingBoundary=currentField->type[i];
          updateMovingSpeed(nextField,parameter,t+k,IDmovingBoundary,i);
        }
    }

    // To be modified !! This just entered in the switch about particles type !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //if(currentField->nMoving != 0){updateMovingSpeed(nextField,parameter,t+k);}
    // Pressure (all particles at the same time)
    pressureComputation(nextField,parameter);
}

/*
*Input:
*- currentField: field that contains all the variables
*- parameter: pointer to the field containing the user defined parameters
*- boxes: vector of vectors of int, each vector is related to a given box and contains a list
with the particles ID of the particles that are inside this box
*- surrBoxesAll: vector of vector of int, each vector is related to a given box and contains
a list with the box ID of the boxes that are adjacent to this box
*- currentDensityDerivative: vector containing derivative of density for each particle at time t
*- currentSpeedDerivative: vector containing derivative of velocity for each particle at time t
*- timeInfo: pointer to the array containing the duration of each part of the code
*Description:
* Knowing the field (currentField), computes the density and velocity derivatives and store them in vectors
*/
void derivativeComputation(Field* currentField, Parameter* parameter, SubdomainInfo &subdomainInfo, std::vector<std::vector<int> >& boxes, std::vector<std::vector<int> >& surrBoxesAll, std::vector<double>& currentDensityDerivative, std::vector<double>& currentSpeedDerivative, std::vector<double> &timeInfo, bool midPoint)
{
  // CPU time information
  std::clock_t start;
  // Neighbors vectors (declaration outside)
  std::vector<int> neighbors; // ICI
  std::vector<double> kernelGradients;

  // Sort the particles at the current time step
  start = getTime();
  if(!midPoint){sortParticles(currentField->pos, currentField->l, currentField->u, subdomainInfo.boxSize, boxes);} // At each time step, restart it (to optimize with lists?)
  timeInfo[1] += ( getTime() - start ) / (double) CLOCKS_PER_SEC;

  // Spans the boxes
  for(int box=subdomainInfo.startingBox ; box<=subdomainInfo.endingBox ; box++){
    // Spans the particles in the box
    for(unsigned int part=0 ; part<boxes[box].size() ; part++){
      start = getTime();

      // Declarations
      int particleID = boxes[box][part];
      neighbors.resize(0); // ICI
      kernelGradients.resize(0);

      // Neighbor search
      findNeighbors(particleID, currentField->pos, parameter->kh, boxes, surrBoxesAll[box], neighbors, kernelGradients, parameter->kernel);
      timeInfo[1] += ( getTime() - start ) / (double) CLOCKS_PER_SEC;

      // Continuity equation
      start = getTime();
      currentDensityDerivative[particleID] = continuity(particleID, neighbors, kernelGradients,currentField);
      timeInfo[2] += ( getTime() - start ) / (double) CLOCKS_PER_SEC;

      // Momentum equation only for free particles
      start = getTime();
      if(currentField->type[particleID] == freePart)
        momentum(particleID, neighbors, kernelGradients, currentField, parameter, currentSpeedDerivative);
      timeInfo[3] += ( getTime() - start ) / (double) CLOCKS_PER_SEC;
    }
  }

}

/*
*Input:
*- currentField: field that contains all the information about step n-1
*- nextField: field in which results of step n are stored
*- parameter: pointer to the field containing the user defined parameters
*- boxes: vector of vectors of int, each vector is related to a given box and contains a list
with the particles ID of the particles that are inside this box
*- surrBoxesAll: vector of vector of int, each vector is related to a given box and contains
a list with the box ID of the boxes that are adjacent to this box
*- n: number of the current time step
*- timeInfo: pointer to the array containing the duration of each part of the code
*Output:
*- Reboxing: flag that indicates if the box division need to be recomputed
*Description:
* Knowing the field at time t(currentField), computes the field at time t+k with euler integration method and store it in structure nextField
*/
void timeIntegration(Field* currentField, Field* nextField, Parameter* parameter,
  SubdomainInfo &subdomainInfo, std::vector<std::vector<int> >& boxes, std::vector<std::vector<int> >& surrBoxesAll,
  double t, double k, std::vector<double> &timeInfo)
  {
    std::vector<double> currentSpeedDerivative;
    std::vector<double> currentDensityDerivative;
    currentSpeedDerivative.assign(3*currentField->nTotal, 0.0);
    currentDensityDerivative.assign(currentField->nTotal, 0.0);

    // CPU time information
    std::clock_t start;
    derivativeComputation(currentField, parameter, subdomainInfo, boxes, surrBoxesAll, currentDensityDerivative, currentSpeedDerivative, timeInfo, false);

    switch (parameter->integrationMethod)
    {
      case euler:
      {
          start = getTime();
          eulerUpdate(currentField, nextField, parameter, subdomainInfo, currentDensityDerivative, currentSpeedDerivative, t, k);
          timeInfo[4] += ( getTime() - start ) / (double) CLOCKS_PER_SEC;
      }
      break;


      case RK2:
      {
          double kMid = 0.5*k/parameter->theta;
          std::vector<double> midSpeedDerivative;
          std::vector<double> midDensityDerivative;
          midSpeedDerivative.assign(3*currentField->nTotal, 0.0);
          midDensityDerivative.assign(currentField->nTotal, 0.0);

          start = getTime();
          // Storing midpoint in nextField
          eulerUpdate(currentField, nextField, parameter, subdomainInfo, currentDensityDerivative, currentSpeedDerivative, t, kMid);
          timeInfo[4] += ( getTime() - start ) / (double) CLOCKS_PER_SEC;

          // Share the mid point
          shareRKMidpoint(*nextField, subdomainInfo);

          // Compute derivatives at midPoint
          derivativeComputation(nextField, parameter, subdomainInfo, boxes, surrBoxesAll, midDensityDerivative, midSpeedDerivative, timeInfo, true);

          // Compute weighted mean derivative then update
          start = getTime();

          for(int i = 0; i < currentField->nTotal ;i++){
            if(currentField->type[i] == freePart){
              for(int j = 0; j < 3;j++){
                currentSpeedDerivative[3*i + j] =  (1-parameter->theta)*currentSpeedDerivative[3*i + j] + parameter->theta*midSpeedDerivative[3*i + j];
              }
            }
            currentDensityDerivative[i] =  (1-parameter->theta)*currentDensityDerivative[i] + parameter->theta*midDensityDerivative[i];
          }

          eulerUpdate(currentField, nextField, parameter, subdomainInfo, currentDensityDerivative, currentSpeedDerivative, t, k);

          timeInfo[4] += ( getTime() - start ) / (double) CLOCKS_PER_SEC;

      }
      break;
    }
    return;
  }
