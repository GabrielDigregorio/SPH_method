///**************************************************************************
/// SOURCE: Function to integrate one time step.
///**************************************************************************
#include "Main.h"
#include "Physics.h"
#include "Tools.h"

/*
*Input:
*- currentField: field that contains all the information about step n-1
*- nextField: field in which results of step n are stored
*- parameter: pointer to the field containing the user defined parameters
*- boxes: vector of vector????????????????????????????????????????????????????????????????????????????
*- surrBoxesAll: vector of vector????????????????????????????????????????????????????????????????????
*- n: number of the current time step
*- timeInfo: pointer to the array containing the duration of each part of the code
*Output:
*- Reboxing: flag that indicates if the box division need to be recomputed
*Decscription:
* Knowing the field at timestep n-1 (currentField), computes the field at time n and store it in structure nextField
*/
bool timeIntegration(Field* currentField, Field* nextField,
    Parameter* parameter, std::vector<std::vector<int> >& boxes, std::vector<std::vector<int> >& surrBoxesAll,
    unsigned int n, std::vector<double> &timeInfo)
{
    // CPU time information
    std::clock_t start;

    // Sort the particles at the current time step
    start = std::clock();
    boxClear(boxes); // Clear the sorting to restart it
    sortParticles(currentField->pos, currentField->l, currentField->u, parameter->kh, boxes); // At each time step (to optimize?)
    timeInfo[1] += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

    // Runge-Kutta 2 parameters ?????????????????? What for ?
    double k1_rho, k2_rho;
    std::vector<double> k1_U, k2_U;

    // Spans the boxes
    for(int box=0 ; box<boxes.size() ; box++)
    {
        // Spans the particles in the box
        for(unsigned int part=0 ; part<boxes[box].size() ; part++)
        {
            start = std::clock();

            // Declarations
            int particleID = boxes[box][part];
            std::vector<int> neighbors;
            std::vector<double> kernelGradients;
            std::vector<double> speedDerivative;//Size is invariant and equal to 3, we could replace it by vector and make momentum return a double* ?

            // Neighbor search
            findNeighbors(particleID, currentField->pos, parameter->kh, boxes, surrBoxesAll[box], neighbors, kernelGradients, parameter->kernel);

            timeInfo[1] += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

            // Continuity equation
            start = std::clock();
            double densityDerivative = continuity(particleID, neighbors, kernelGradients,currentField);
            timeInfo[2] += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

            // Momentum equation only for free particles
            start = std::clock();
            if(particleID < currentField->nFree)
                momentum(particleID, neighbors, kernelGradients, currentField, parameter, speedDerivative);
            timeInfo[3] += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

            // Integration
            start = std::clock();
            switch(parameter->integrationMethod)
            {
                case euler: // u_n = u_(n-1) + k * du/dt
                nextField->density[particleID] = currentField->density[particleID] + parameter->k*densityDerivative;
                // Update speed only for Free particles
                if(particleID < currentField->nFree)
                {
                    for (int i = 0; i <= 2; i++)
                    {
                        nextField->speed[3*particleID + i] = currentField->speed[3*particleID + i] + parameter->k*speedDerivative[i];
                        //std::cout << speedDerivative[i] << std::endl;!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        //std::cout << currentField->speed[3*particleID + i] << std::endl;//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    }
                }
                break;

                case RK2:
                std::cout << "Integration method not coded.\n";
                break;

                default:// Default will deseaper when consistencyCheck will be coded
                std::cout << "Integration method not coded.\n";
                break;
            }
            // Position ( update only for non fixed particles )
            if( (particleID < currentField->nFree) || (particleID >= currentField->nFree + currentField->nFixed) )
            {
                for (int i = 0; i <= 2; i++)
                {
                    nextField->pos[3*particleID + i] = currentField->pos[3*particleID + i] + parameter->k*currentField->speed[3*particleID + i];
                    //std::cout << currentField->speed[3*particleID + i] << std::endl;//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                }
            }
            timeInfo[4] += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
        }
    }

    start = std::clock();
    // Pressure
    pressureComputation(nextField,parameter);

    // Speed (for all moving particles)
    if(currentField->nMoving != 0){updateMovingSpeed(nextField,parameter,n*parameter->k);}

    // Reboxing criterion
    bool reBoxing = false; // A fonction should be implemented to choose if we rebox or not

    timeInfo[4] += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

    return reBoxing;
}
