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
void eulerUpdate(Field* currentField, Field* nextField,Parameter* parameter, std::vector<double>& currentDensityDerivative, std::vector<double>& currentSpeedDerivative, double t, double k)
{
  // Free particles
  for(int i = 0; i < currentField->nFree; i++)
  {
    nextField->density[i] = currentField->density[i] + k*currentDensityDerivative[i];
    for (int j = 0; j <= 2; j++)
    {
      nextField->speed[3*i + j] = currentField->speed[3*i + j] + k*currentSpeedDerivative[3*i + j];
      nextField->pos[3*i + j] = currentField->pos[3*i + j] + k*currentField->speed[3*i + j];
    }
  }

  // Fixed particles
  for(int i = currentField->nFree; i < (currentField->nFree + currentField->nFixed); i++)
  {
    nextField->density[i] = currentField->density[i] + k*currentDensityDerivative[i];
  }

  // Moving particles
  for(int i = (currentField->nFree + currentField->nFixed); i < currentField->nTotal; i++)
  {
    nextField->density[i] = currentField->density[i] + k*currentDensityDerivative[i];
    for (int j = 0; j <= 2; j++)
    {
      nextField->pos[3*i + j] = currentField->pos[3*i + j] + k*currentField->speed[3*i + j];
    }
  }
  if(currentField->nMoving != 0){updateMovingSpeed(nextField,parameter,t+k);}

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
*Decscription:
* Knowing the field (currentField), computes the density and velocity derivatives and store them in vectors
*/
void derivativeComputation(Field* currentField, Parameter* parameter, std::vector<std::vector<int> >& boxes, std::vector<std::vector<int> >& surrBoxesAll, std::vector<double>& currentDensityDerivative, std::vector<double>& currentSpeedDerivative, std::vector<double> &timeInfo)
{
  // CPU time information
  std::clock_t start;

  // Sort the particles at the current time step
  start = std::clock();
  boxClear(boxes); // Clear the sorting to restart it
  sortParticles(currentField->pos, currentField->l, currentField->u, parameter->kh, boxes); // At each time step (to optimize?)
  timeInfo[1] += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

  // Spans the boxes
  for(int box=0 ; box<boxes.size() ; box++){
    // Spans the particles in the box
    for(unsigned int part=0 ; part<boxes[box].size() ; part++){
      start = std::clock();

      // Declarations
      int particleID = boxes[box][part];
      std::vector<int> neighbors;
      std::vector<double> kernelGradients;

      // Neighbor search
      findNeighbors(particleID, currentField->pos, parameter->kh, boxes, surrBoxesAll[box], neighbors, kernelGradients, parameter->kernel);

      timeInfo[1] += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

      // Continuity equation
      start = std::clock();
      currentDensityDerivative[particleID] = continuity(particleID, neighbors, kernelGradients,currentField);
      timeInfo[2] += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

      // Momentum equation only for free particles
      start = std::clock();
      if(particleID < currentField->nFree)
      momentum(particleID, neighbors, kernelGradients, currentField, parameter, currentSpeedDerivative);
      timeInfo[3] += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
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
*Decscription:
* Knowing the field at time t(currentField), computes the field at time t+k with euler integration method and store it in structure nextField
*/
bool timeIntegration(Field* currentField, Field* nextField,
  Parameter* parameter, std::vector<std::vector<int> >& boxes, std::vector<std::vector<int> >& surrBoxesAll,
  double t, double k, std::vector<double> &timeInfo)
  {
    std::vector<double> currentSpeedDerivative;
    std::vector<double> currentDensityDerivative;
    currentSpeedDerivative.assign(3*currentField->nFree, 0.0);
    currentDensityDerivative.assign(currentField->nTotal, 0.0);

    // CPU time information
    std::clock_t start;

    derivativeComputation(currentField, parameter, boxes, surrBoxesAll, currentDensityDerivative, currentSpeedDerivative, timeInfo);

    switch (parameter->integrationMethod)
    {
      case euler:
      {
          start = std::clock();
          eulerUpdate(currentField, nextField, parameter, currentDensityDerivative, currentSpeedDerivative, t, k);
          timeInfo[4] += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
      }
      break;


      case RK2:
      {
          double kMid = 0.5*k/parameter->theta;
          std::vector<double> midSpeedDerivative;
          std::vector<double> midDensityDerivative;
          midSpeedDerivative.assign(3*currentField->nFree, 0.0);
          midDensityDerivative.assign(currentField->nTotal, 0.0);

          start = std::clock();
          // Storing midpoint in nextField
          eulerUpdate(currentField, nextField, parameter, currentDensityDerivative, currentSpeedDerivative, t, kMid);
          timeInfo[4] += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

          // Compute derivatives at midPoint
          derivativeComputation(nextField, parameter, boxes, surrBoxesAll, midDensityDerivative, midSpeedDerivative, timeInfo);

          // Compute weighted mean derivative then update
          start = std::clock();

          for(int i = 0; i < currentField->nFree; i++)
          {
            for(int j = 0; j < 3;j++)
            {
              currentSpeedDerivative[3*i + j] =  (1-parameter->theta)*currentSpeedDerivative[3*i + j] + parameter->theta*midSpeedDerivative[3*i + j];
            }
          }
          for(int i = 0; i < currentField->nTotal ;i++)
          {
            currentDensityDerivative[i] =  (1-parameter->theta)*currentDensityDerivative[i] + parameter->theta*midDensityDerivative[i];
          }

          eulerUpdate(currentField, nextField, parameter, currentDensityDerivative, currentSpeedDerivative, t, k);

          timeInfo[4] += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;


      }
      break;
    }

    // Reboxing criterion
    start = std::clock();
    bool reBoxing = false; // A fonction should be implemented to choose if we rebox or not
    timeInfo[4] += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

    return reBoxing;
  }
