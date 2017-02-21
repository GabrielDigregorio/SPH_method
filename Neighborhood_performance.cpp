/*
 * To be finished.
 */

#include "SPH.h"
#include <fstream>
#include <vector>
#include <ctime>

int main(int argc, char *argv[])
{
    //Should be in argv[]
    double o[3] = {0.0,0.0,0.0};
    double L[3] = {10.0,10.0,10.0};
    double s = 1.0;
    //End input (end of what should be in argv

    std::vector<double> pos;

    //Generate cube
    meshcube(o,L,s,pos);

    //Record algorithm performance
    std::clock_t start;
    double duration;
    start = std::clock();

    //Apply sorting

    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

    std::cout<<"Elapsed time: "<< duration <<" [s]\n";

    return 0;
}
