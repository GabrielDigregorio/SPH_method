#include "Main.h"
#include "Physics.h"
#include "Tools.h"


bool timeIntegration(Field* currentField, Field* nextField, Parameter* parameter, std::vector<std::vector<int> >& boxes, std::vector<std::vector<int> >& surrBoxesAll, unsigned int n)
{


    // Time step resolution
    //double kh2 = parameter->kh*parameter->kh; // More efficient to compare distance^2
    Kernel kernelType = parameter->kernel;

    // Sort the particles at the current time step
    std::cout << "\t Sorting particles...\n" << std::endl;
    boxClear(boxes); // Clear the sorting to restart it...
    sortParticles(currentField->pos, currentField->l, currentField->u, parameter->kh, boxes); // At each time step (to optimize?)

    // Declaration
    std::vector<int> neighbors;
    std::vector<double> kernelGradients;
    int particleID;
    double densityDerivative;
    std::vector<double> speedDerivative;

    // runge kutta 2 parameters
    double k1_rho, k2_rho;
    std::vector<double> k1_U, k2_U;

    // Spans the boxes
    for(int box=0 ; box<boxes.size() ; box++)
    {
      // Spans the particles in the box
      for(unsigned int part=0 ; part<boxes[box].size() ; part++)
      {
        particleID = boxes[box][part];

        //std::cout << "\t particleID = " << particleID <<":\n";

        //std::cout << "\t \t Finding neighbors..." << std::endl;
        neighbors.resize(0);// VERY BAD
        kernelGradients.resize(0);// VERY BAD 
        findNeighbors(particleID, currentField->pos, parameter->kh, boxes, surrBoxesAll[box], neighbors, kernelGradients, kernelType);
        //std::cout << neighbors.size() << " neighbors found... \n" << std::endl;
        // Continuity equation
        densityDerivative = continuity(particleID, neighbors, kernelGradients,currentField); // also for fixed particles !
        //std::cout << "\t \t Density derivative = " << densityDerivative << "\n";
        // Momentum equation only for free particles
        if(particleID < currentField->nFree)
        {
          momentum(particleID, neighbors, kernelGradients,currentField,parameter, speedDerivative);
          //std::cout << "\t \t Speed derivative = " << speedDerivative[0] << "\n";
          //std::cout << "\t \t Speed derivative = " << speedDerivative[1] << "\n";
          //std::cout << "\t \t Speed derivative = " << speedDerivative[2] << "\n";
        }


        // Integration
        switch ( parameter->integrationMethod )
        {
          case euler: // u_n = u_(n-1) + k * du/dt
          nextField->density[particleID] = currentField->density[particleID] + parameter->k*densityDerivative;



          // Uptade speed only for Free particles
          if(particleID < currentField->nFree)
          {
            for (int i = 0; i <= 2; i++)
            {
              nextField->speed[3*particleID + i] = currentField->speed[3*particleID + i] + parameter->k*speedDerivative[i];
            }
          }

          break;
          case RK2:
          std::cout << "Runge Kutta 2 not coded.\n";
          break;

          default:
          std::cout << "Integration method not coded.\n";
          return EXIT_FAILURE;
          break;
        }
        // MAJ of different Field
        //Position ( MAJ only for non Fixed particles )
        if( (particleID < currentField->nFree) || (particleID >= currentField->nFree + currentField->nFixed) )
        {
          for (int i = 0; i <= 2; i++)
          {
            nextField->pos[3*particleID + i] = currentField->pos[3*particleID + i] + parameter->k*currentField->speed[3*particleID + i];
          }
        }

      }
    }
    //pressure
    pressureComputation(nextField,parameter);

    //Update speed of all moving particles
    updateMovingSpeed(nextField,parameter,n*parameter->k);

    bool reBoxing = false; // A fonction should be implemented to choose if we rebox or not
    return reBoxing;
}
