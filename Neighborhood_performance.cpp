#include "SPH.h"
<<<<<<< HEAD
#include <fstream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cmath>
=======

>>>>>>> 1216497b23b0f67d0cc2c237dd404a6cf84ef544

int main(int argc, char *argv[])
{
    //Input parameters
    double s = atof(argv[1]);
    double kh = atof(argv[2]);
    double l = atof(argv[3]);
    std::cout << "\n Parameter list: " << s << ", " << kh << ", " << l << "\n";

    double o[3] = {0.0,0.0,0.0};
<<<<<<< HEAD
    double L[3] = {l,l,l};
=======
    double L[3] = {10.0,10.0,10.0};
    double s = 1.4;
    double kh = 2.2;
    //End input (end of what should be in argv)
>>>>>>> 1216497b23b0f67d0cc2c237dd404a6cf84ef544

    std::vector<double> pos;
    std::vector<double> values;
    std::vector<int> row;
    std::vector<int> column;

    //Generate cube
    meshcube(o,L,s,pos, 0);

    //Record algorithm performance
    std::clock_t start;
    double duration;

    start = std::clock();
    neighborAllPair(pos, kh, values, row, column);
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"Elapsed time AllPair: "<< duration <<" [s]\n";

    start = std::clock();
    neighborLinkedList(pos,o, L, kh, values, row, column);
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"Elapsed time Linked List: "<< duration <<" [s]\n";


    std::cout << "\n";
    return 0;
}
