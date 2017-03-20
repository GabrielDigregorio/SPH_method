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
    for(int box=0 ; box<boxes.size() ; box++){
        // Spans the particles in the box
        for(unsigned int part=0 ; part<boxes[box].size() ; part++){
            start = std::clock();

            // Declarations
            int particleID = boxes[box][part];
            std::vector<int> neighbors;
            std::vector<double> kernelGradients;
            std::vector<double> speedDerivative;//Size is invariant and equal to 3, we could replace it by vector and make momentum return a double* ?

            // Neighbor search
            findNeighbors(particleID, currentField->pos, parameter->kh, boxes, surrBoxesAll[box], neighbors, kernelGradients, parameter->kernel);

            timeInfo[1] += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
            
            switch(parameter->integrationMethod){
                case euler:
                { 
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
                    //switch(parameter->integrationMethod){
                    //  case euler: // u_n = u_(n-1) + k * du/dt
                    nextField->density[particleID] = currentField->density[particleID] + parameter->k*densityDerivative;
                    // Update speed only for Free particles
                    if(particleID < currentField->nFree){
                        for (int i = 0; i <= 2; i++){
                            nextField->speed[3*particleID + i] = currentField->speed[3*particleID + i] + parameter->k*speedDerivative[i];
                        }
                    } 
                }
                break;

                case RK2:
                {
                    // alpha=1: heuns
                    // alpha=1/2: midpoint
                    // alpha=2/3: Ralston
                    double alpha=1.0; // to choose via struc param 
                    // see wikipedia 
                    // 1: 0  0
                    // 2:alpha 0
                    // next= 1-1/(2*alpha)   1/(2*alpha)
                    // romain Runge kutta, 2 step 

                    // // continuity
                    double k=parameter->k;
                    std::vector<double> u_vector_1={0.0,0.0,0.0};
                    std::vector<double> u_vector_2={0.0,0.0,0.0};
                    std::vector<double> speedDerivative={0.0,0.0,0.0};
                    double temp_rho_2=0.0;
                    double temp_rho_1=0.0;
                    double rho_1=0.0;
                    double rho_2=0.0;
                    // evaluating speed field for the current particle 
                
                    for(int i=0; i<3; ++i){
                       u_vector_1[i]=currentField->speed[3*particleID+i];
                    }
                    momentum(particleID, neighbors, kernelGradients,currentField,parameter, speedDerivative);
                    for(int i=0; i<3; ++i){
                        u_vector_2[i]=currentField->speed[3*particleID+i]+alpha*k*speedDerivative[i];
                    }
                    
                    // calucling intermediate value for the RK2 formula
                    temp_rho_1=continuity(particleID, neighbors, kernelGradients,currentField);
                    for(int i=0; i<3; ++i){
                        currentField->speed[3*particleID+i]  =u_vector_2[i];
                    }
                    temp_rho_2=continuity(particleID, neighbors, kernelGradients,currentField);

                    // getting back normal veocity field (it should not be affected by the calculation of rho for instance)
                    for(int i=0; i<3; ++i){
                        currentField->speed[3*particleID+i]  =u_vector_1[i];
                    }
                    // calculating new density 
                    nextField->density[particleID]= currentField->density[particleID]+ (1.0-1.0/(2.0*alpha))*temp_rho_1+1.0/(2.0*alpha)*temp_rho_2;
                     
                     
                    // // momentum
                    if(particleID < currentField->nFree)
                    {
                            rho_1=currentField->density[particleID];
                            rho_2=currentField->density[particleID]+alpha*k*continuity(particleID, neighbors, kernelGradients,currentField);
        
                            // calucling intermediate value for the RK2 formul
                            std::vector<double> temp_u_vector_1={0.0,0.0,0.0};
                            std::vector<double> temp_u_vector_2={0.0,0.0,0.0};
                            momentum(particleID, neighbors, kernelGradients,currentField,parameter, speedDerivative);
                            
                            for(int i=0; i<3; ++i){
                                temp_u_vector_1[i]=speedDerivative[i];
                            }

                            currentField->density[particleID]=rho_2;
                            momentum(particleID, neighbors, kernelGradients,currentField,parameter, speedDerivative);
                            for(int i=0; i<3; ++i){
                                temp_u_vector_2[i]=speedDerivative[i];
                            }
                            // getting back normal density (it should not be affected by the calculation of speed for instance)
                            currentField->density[particleID]=rho_1;

                            // calculating new speed
                            for(int i=0; i<3; ++i){
                                nextField->speed[3*particleID+i]  =currentField->speed[3*particleID+i]+(1.0-1.0/(2.0*alpha))*temp_u_vector_1[i]+1.0/(2.0*alpha)*temp_u_vector_2[i];
                            }
                    }       
                } 
            
            break;

            default:
                std::cout << "Integration method not coded.\n";
            return EXIT_FAILURE;
            }
            // Position ( update only for non fixed particles )
            if( (particleID < currentField->nFree) || (particleID >= currentField->nFree + currentField->nFixed) )
            {
                for (int i = 0; i <= 2; i++)
                {
                    nextField->pos[3*particleID + i] = currentField->pos[3*particleID + i] + parameter->k*currentField->speed[3*particleID + i];
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