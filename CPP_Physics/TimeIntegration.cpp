#include "Main.h"
#include "Physics.h"

bool timeIntegration(Field* currentField, Field* nextField, Parameter* parameter, unsigned int n){
  /*  // Time step resolution
    double kh2 = parameter->kh*parameter->kh; // More efficient to compare distance^2

    // Sort the particles at the current time step
    boxClear(boxes); // Clear the sorting to restart it...
    sortParticles(pos, l, u, kh, boxes); // At each time step (to optimize?)
    // Spans the boxes
    for(int box=0 ; box<boxes.size() ; box++)
    {
        // Spans the particles in the box
        for(unsigned int part=0 ; part<boxes[box].size() ; part++){
            // Declaration
            std::vector<int> neighbors;
            std::vector<double> kernelGradients;
            int particleID = boxes[box][part];
            findNeighbors(particleID, pos, kh2, boxes, surrBoxesAll[box], neighbors, kernelGradients);

            //continuity(); // also for fixed particles !
            // IF free particle
            //momentum();
        }
    }
    bool reBoxing = false; // A fonction should be implemented to choose if we rebox or not
    return reBoxing;
    */
}
