#include "../Headers/SPH.hpp"

void timeIntegration(playground, ){
    // Time step resolution

    // declaration of boxes, neighbors


    // LOOP on time
        sortParticles(pos, l, u, kh, boxes);
        // LOOP on particles
            findNeighbors(boxes, surrBoxes, pos, kh, neighbors);
            continuity(); // also for fixed particles !
            // IF free particle
                momentum();
        // END LOOP on particles
    // END LOOP on time

}
