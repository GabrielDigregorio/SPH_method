#include "Main.h"
#include "Physics.h"

void timeIntegration(Field* currentField, Field* nextField, Parameter* parameter, unsigned int n){
    /*
    // Time step resolution
    double kh2 = kh*kh; // More efficient to compare distance^2

    // Creates the box mesh and describes their adjacent relations
    std::vector<std::vector<int> > boxes;
    std::vector<std::vector<int> > surrBoxesAll;
    boxMesh(l, u, kh, boxes, surrBoxesAll); // Once for the whole resolution

    // Time integration loop
    for(int iteration=0 ; iteration < 10 ; iteration++){



        // Sort the particles at the current time step
        boxClear(boxes); // Clear the sorting to restart it...
        sortParticles(pos, l, u, kh, boxes); // At each time step (to optimize?)
        // Spans the boxes
        for(int box=0 ; box<boxes.size() ; box++){
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



    }
    return;
    */
}
